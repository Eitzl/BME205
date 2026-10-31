"""Helper code for the breast cancer PRS notebook.

The notebook holds the knobs (the CFG dict). This file holds the machinery:
loading data, in-fold preprocessing, models, cross-validation, logging,
sweeps and writing submissions.

Rule followed everywhere: anything that "learns" from data (scaling,
feature selection, PCA, ancestry stats) is fit on the training fold only,
so the cross-validation AUC is an honest estimate of the leaderboard.
"""

from __future__ import annotations

import copy
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler

ID_NAMES = {"id", "sample_id", "sampleid", "sample", "iid", "eid", "row_id"}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data(cfg: dict) -> dict:
    """Read train/test CSVs and split them into PRS matrix, ancestry and target."""
    dcfg = cfg["data"]
    d = _find_data_dir(dcfg["data_dir"], dcfg["train_file"])
    cache = d / "_cache.pkl"
    if dcfg.get("use_cache") and cache.exists():
        train, test = pd.read_pickle(cache)
    else:
        train = pd.read_csv(d / dcfg["train_file"])
        test = pd.read_csv(d / dcfg["test_file"])
        if dcfg.get("use_cache"):
            pd.to_pickle((train, test), cache)

    target = dcfg["target"]
    id_col = dcfg.get("id_col") or _detect_id(train, test, target)
    anc_col = dcfg.get("ancestry_col") or _detect_ancestry(train, target, id_col)

    skip = {target, id_col, anc_col} - {None}
    prs_cols = [c for c in train.columns if c not in skip]

    rng = np.random.default_rng(cfg["seed"])
    if dcfg.get("debug_max_features"):
        k = min(dcfg["debug_max_features"], len(prs_cols))
        prs_cols = sorted(rng.choice(prs_cols, k, replace=False).tolist(),
                          key=prs_cols.index)
    if dcfg.get("debug_max_rows") and dcfg["debug_max_rows"] < len(train):
        train = train.sample(dcfg["debug_max_rows"], random_state=cfg["seed"])

    fill = dcfg.get("fillna", 0.0)
    X = train[prs_cols].to_numpy(np.float32)
    Xt = test[prs_cols].to_numpy(np.float32)
    X = np.where(np.isnan(X), fill, X)
    Xt = np.where(np.isnan(Xt), fill, Xt)

    if anc_col:
        a = train[anc_col].astype(str).to_numpy()
        at = test[anc_col].astype(str).to_numpy()
    else:
        a = np.full(len(train), "all")
        at = np.full(len(test), "all")

    y = train[target].to_numpy().astype(int)
    test_ids = test[id_col].to_numpy() if id_col else np.arange(len(test))

    print(f"train: {X.shape[0]:,} rows | test: {Xt.shape[0]:,} rows | "
          f"PRS columns: {len(prs_cols):,}")
    print(f"id column: {id_col!r} | ancestry column: {anc_col!r} | "
          f"cases: {y.sum():,} ({y.mean():.1%})")
    return dict(X=X, a=a, y=y, X_test=Xt, a_test=at, test_ids=test_ids,
                prs_cols=prs_cols, id_col=id_col, ancestry_col=anc_col)


def _find_data_dir(data_dir, train_file):
    """Use data_dir if it has the files, else try the same name one level up
    (for a repo-root data/ folder)."""
    for d in (Path(data_dir), Path("..") / data_dir):
        if (d / train_file).exists():
            return d
    raise FileNotFoundError(
        f"Couldn't find {train_file} in {Path(data_dir).resolve()} or "
        f"{(Path('..') / data_dir).resolve()}. Set CFG['data']['data_dir'].")


def _detect_id(train, test, target):
    for c in test.columns:
        if c.lower() in ID_NAMES:
            return c
    return None


def _detect_ancestry(train, target, id_col):
    cols = [c for c in train.columns if c not in (target, id_col)]
    for c in cols:
        if "ancestry" in c.lower() or c.lower() in {"pop", "population", "superpop"}:
            return c
    for c in cols:
        if not pd.api.types.is_numeric_dtype(train[c]):
            return c
    return None


# ---------------------------------------------------------------------------
# Univariate AUC (fast, column-wise)
# ---------------------------------------------------------------------------
def fast_auc(X: np.ndarray, y: np.ndarray, chunk: int = 500) -> np.ndarray:
    """AUC of every column of X on its own (Mann-Whitney U via ranks)."""
    y = np.asarray(y).astype(bool)
    n1, n0 = y.sum(), (~y).sum()
    out = np.empty(X.shape[1])
    for s in range(0, X.shape[1], chunk):
        r = rankdata(X[:, s:s + chunk], axis=0)
        out[s:s + chunk] = (r[y].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
    return out


def univariate_table(data: dict) -> pd.DataFrame:
    """AUC of each PRS column alone, overall and per ancestry group."""
    X, y, a = data["X"], data["y"], data["a"]
    df = pd.DataFrame({"feature": data["prs_cols"], "auc": fast_auc(X, y)})
    df["strength"] = (df["auc"] - 0.5).abs()
    for g in np.unique(a):
        m = a == g
        if 0 < y[m].sum() < m.sum():
            df[f"auc_{g}"] = fast_auc(X[m], y[m])
    return df.sort_values("strength", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Preprocessing (fit on the training fold, applied to val/test)
# ---------------------------------------------------------------------------
class Preprocessor:
    def __init__(self, pcfg: dict, seed: int = 0):
        self.p = pcfg
        self.seed = seed

    def fit(self, X, a, y):
        p = self.p
        X = X.astype(np.float32, copy=True)

        # 1. Re-standardize each PRS within each ancestry group
        self.groups_ = np.unique(a)
        self.global_ = (X.mean(0), _safe_std(X))
        self.stats_ = {}
        if p["restandardize_by_ancestry"]:
            for g in self.groups_:
                m = a == g
                if m.sum() >= p["min_group_size"]:
                    self.stats_[g] = (X[m].mean(0), _safe_std(X[m]))
        X = self._standardize(X, a)

        # 2. Clip extreme Z-scores
        if p["clip_z"]:
            np.clip(X, -p["clip_z"], p["clip_z"], out=X)

        # 3. Univariate strength = |AUC - 0.5| (direction ignored)
        auc = fast_auc(X, y)
        strength = np.abs(auc - 0.5)
        order = np.argsort(-strength)
        if p["min_univariate_strength"]:
            order = order[strength[order] >= p["min_univariate_strength"]]
        if p["top_k_features"]:
            order = order[:p["top_k_features"]]

        # 4. Drop near-duplicate columns (keep the stronger one)
        if p["corr_threshold"] and len(order) > 1:
            order = _corr_filter(X[:, order], order, p["corr_threshold"])
        self.keep_ = order
        self.auc_ = auc
        Xk = X[:, order]

        # 5. Optional PCA
        self.pca_ = None
        if p["pca_components"]:
            n = min(p["pca_components"], Xk.shape[1], Xk.shape[0] - 1)
            self.pca_ = PCA(n, random_state=self.seed).fit(Xk)

        # 6. Which PRS get ancestry interaction terms
        self.inter_ = np.arange(min(p["interaction_top_k"] or 0, len(order)))
        return self

    def transform(self, X, a):
        p = self.p
        X = self._standardize(X.astype(np.float32, copy=True), a)
        if p["clip_z"]:
            np.clip(X, -p["clip_z"], p["clip_z"], out=X)
        Xk = X[:, self.keep_]

        parts = []
        if self.pca_ is not None:
            pcs = self.pca_.transform(Xk).astype(np.float32)
            parts += [pcs] if p["pca_mode"] == "replace" else [Xk, pcs]
        else:
            parts.append(Xk)

        if p["ancestry_onehot"] or len(self.inter_):
            dummies = np.stack([(a == g) for g in self.groups_], 1).astype(np.float32)
            if p["ancestry_onehot"]:
                parts.append(dummies)
            if len(self.inter_):
                top = Xk[:, self.inter_]
                parts += [top * dummies[:, [j]] for j in range(dummies.shape[1])]
        return np.hstack(parts)

    def _standardize(self, X, a, cols=slice(None)):
        if not self.stats_:
            return X
        mu, sd = self.global_[0][cols], self.global_[1][cols]
        out = (X - mu) / sd if self.p["standardize_unseen_globally"] else X.copy()
        for g, (gm, gs) in self.stats_.items():
            m = a == g
            out[m] = (X[m] - gm[cols]) / gs[cols]
        return out

    def best_single(self, X, a):
        """Score = the strongest PRS on the training fold, flipped if protective."""
        j = self.keep_[0]
        z = self._standardize(X[:, [j]].astype(np.float32), a, cols=[j])[:, 0]
        return z if self.auc_[j] >= 0.5 else -z

    def n_features_out(self):
        return int(len(self.keep_))


def _safe_std(X):
    s = X.std(0)
    s[s < 1e-6] = 1.0
    return s


def _corr_filter(Xs, order, thr):
    Z = (Xs - Xs.mean(0)) / _safe_std(Xs)
    C = np.abs((Z.T @ Z) / len(Z))
    dropped = np.zeros(len(order), bool)
    kept = []
    for j in range(len(order)):
        if dropped[j]:
            continue
        kept.append(j)
        dropped |= C[j] > thr
    return order[kept]


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
def fit_predict_best_single(m, Xtr, ytr, Xva_list, seed, pre=None, raw=None):
    return [pre.best_single(X, a) for X, a in raw], {}


def fit_predict_logreg(m, Xtr, ytr, Xva_list, seed):
    scaler = StandardScaler().fit(Xtr) if m["scale"] else None
    tf = (lambda Z: scaler.transform(Z)) if scaler else (lambda Z: Z)
    clf = LogisticRegression(
        C=m["C"], l1_ratio=m["l1_ratio"], solver=m["solver"],
        max_iter=m["max_iter"], tol=m["tol"], class_weight=m["class_weight"],
        random_state=seed,
    ).fit(tf(Xtr), ytr)
    info = {"nonzero_coefs": int((clf.coef_ != 0).sum())}
    return [clf.predict_proba(tf(Z))[:, 1] for Z in Xva_list], info


def fit_predict_lgbm(m, Xtr, ytr, Xva_list, seed):
    import lightgbm as lgb

    params = dict(m["params"], objective="binary", metric="auc",
                  seed=seed, verbose=-1)
    # Early stopping uses a slice of the TRAINING fold, never the val fold.
    Xa, Xb, ya, yb = train_test_split(Xtr, ytr, test_size=m["early_stopping_frac"],
                                      stratify=ytr, random_state=seed)
    booster = lgb.train(
        params, lgb.Dataset(Xa, ya), num_boost_round=m["num_boost_round"],
        valid_sets=[lgb.Dataset(Xb, yb)],
        callbacks=[lgb.early_stopping(m["early_stopping_rounds"], verbose=False)],
    )
    it = booster.best_iteration or m["num_boost_round"]
    return [booster.predict(Z, num_iteration=it) for Z in Xva_list], {"best_iter": it}


MODEL_FNS = {"best_single": fit_predict_best_single,
             "logreg": fit_predict_logreg, "lgbm": fit_predict_lgbm}


def blend(preds: dict, ecfg: dict) -> np.ndarray:
    """Combine model predictions. Rank-mean suits AUC (only order matters)."""
    w = {k: ecfg["weights"].get(k, 1.0) for k in preds}
    w = {k: v for k, v in w.items() if v > 0}
    if not w:
        raise ValueError("All ensemble weights are 0.")
    tot = sum(w.values())
    if ecfg["method"] == "rank_mean":
        return sum(w[k] * rankdata(preds[k]) / len(preds[k]) for k in w) / tot
    return sum(w[k] * preds[k] for k in w) / tot


# ---------------------------------------------------------------------------
# Cross-validation
# ---------------------------------------------------------------------------
def run_cv(cfg: dict, data: dict, verbose: bool = True) -> dict:
    X, a, y = data["X"], data["a"], data["y"]
    Xt, at = data["X_test"], data["a_test"]
    models = [k for k, m in cfg["models"].items() if m["enabled"]]
    if not models:
        raise ValueError("Turn on at least one model in CFG['models'].")

    # Stratify on ancestry x case status so each fold looks like the whole set.
    strat = pd.Series(a).astype(str) + "_" + pd.Series(y).astype(str)
    skf = StratifiedKFold(cfg["cv"]["n_folds"], shuffle=True, random_state=cfg["seed"])

    oof = {k: np.zeros(len(y)) for k in models}
    test = {k: np.zeros(len(Xt)) for k in models}
    fold_auc = {k: [] for k in models}
    info = {k: [] for k in models}
    t0 = time.time()

    for f, (tr, va) in enumerate(skf.split(X, strat)):
        if cfg["cv"]["max_folds"] and f >= cfg["cv"]["max_folds"]:
            break
        pre = Preprocessor(cfg["preprocess"], cfg["seed"]).fit(X[tr], a[tr], y[tr])
        Ztr, Zva = pre.transform(X[tr], a[tr]), pre.transform(X[va], a[va])
        Zt = pre.transform(Xt, at)
        msg = [f"fold {f + 1}: {pre.n_features_out()} PRS kept -> {Ztr.shape[1]} cols"]
        for k in models:
            if k == "best_single":
                (pva, pt), inf = fit_predict_best_single(
                    None, None, None, None, None, pre=pre,
                    raw=[(X[va], a[va]), (Xt, at)])
            else:
                (pva, pt), inf = MODEL_FNS[k](cfg["models"][k], Ztr, y[tr],
                                              [Zva, Zt], cfg["seed"] + f)
            oof[k][va] = pva
            test[k] += rankdata(pt) / len(pt)
            fold_auc[k].append(roc_auc_score(y[va], pva))
            info[k].append(inf)
            msg.append(f"{k} {fold_auc[k][-1]:.4f}")
        if verbose:
            print(" | ".join(msg) + f"  ({time.time() - t0:.0f}s)")

    n_done = len(next(iter(fold_auc.values())))
    used = np.zeros(len(y), bool)
    for _, va in itertools.islice(skf.split(X, strat), n_done):
        used[va] = True

    oof = {k: v[used] for k, v in oof.items()}
    yu, au = y[used], a[used]
    oof["ensemble"] = blend(oof, cfg["ensemble"])
    test = {k: v / n_done for k, v in test.items()}
    test["ensemble"] = blend(test, cfg["ensemble"])

    summary = {k: roc_auc_score(yu, v) for k, v in oof.items()}
    by_anc = pd.DataFrame({
        k: {g: _auc_or_nan(yu[au == g], v[au == g]) for g in np.unique(au)}
        for k, v in oof.items()
    })
    by_anc["n"] = pd.Series(au).value_counts()
    by_anc["cases"] = pd.Series(au[yu == 1]).value_counts()

    if verbose:
        print("\nOut-of-fold AUC (higher is better):")
        for k, v in summary.items():
            spread = f" (fold sd {np.std(fold_auc[k]):.4f})" if k in fold_auc else ""
            print(f"  {k:9s} {v:.4f}{spread}")
        print(f"  time: {time.time() - t0:.0f}s")

    return dict(summary=summary, fold_auc=fold_auc, by_ancestry=by_anc,
                oof=oof, y_oof=yu, test=test, info=info, cfg=copy.deepcopy(cfg))


def _auc_or_nan(y, p):
    return roc_auc_score(y, p) if 0 < y.sum() < len(y) else np.nan


# ---------------------------------------------------------------------------
# Logging, submissions, sweeps
# ---------------------------------------------------------------------------
def log_run(res: dict, path: str = "runs_log.csv", note: str = "") -> pd.DataFrame:
    """Append this run (knobs + scores) to a CSV so you can compare runs later."""
    row = {"time": time.strftime("%Y-%m-%d %H:%M:%S"),
           "run_name": res["cfg"]["run_name"], "note": note}
    row.update({f"auc_{k}": round(v, 5) for k, v in res["summary"].items()})
    row.update(_flatten(res["cfg"]))
    df = pd.DataFrame([row])
    p = Path(path)
    if p.exists():
        df = pd.concat([pd.read_csv(p), df], ignore_index=True)
    df.to_csv(p, index=False)
    return df


def write_submission(res: dict, data: dict, cfg: dict, which: str = "ensemble") -> Path:
    out = Path(cfg["output"]["submission_dir"])
    out.mkdir(exist_ok=True)
    path = out / f"{cfg['run_name']}_{which}.csv"
    pd.DataFrame({cfg["output"]["id_col_name"]: data["test_ids"],
                  cfg["output"]["pred_col_name"]: res["test"][which]}
                 ).to_csv(path, index=False)
    print(f"wrote {path}  ({len(data['test_ids']):,} rows)")
    return path


def set_knob(cfg: dict, dotted: str, value) -> dict:
    """set_knob(cfg, 'models.logreg.C', 0.1) changes one nested knob."""
    d = cfg
    *parents, last = dotted.split(".")
    for k in parents:
        d = d[k]
    if last not in d:
        raise KeyError(f"No knob called {dotted!r}")
    d[last] = value
    return cfg


def sweep(base_cfg: dict, data: dict, grid: dict, log_path: str = "runs_log.csv",
          verbose: bool = False) -> pd.DataFrame:
    """Try every combination in grid, e.g. {'models.logreg.C': [0.01, 0.1]}."""
    keys = list(grid)
    rows = []
    combos = list(itertools.product(*grid.values()))
    for i, vals in enumerate(combos, 1):
        cfg = copy.deepcopy(base_cfg)
        for k, v in zip(keys, vals):
            set_knob(cfg, k, v)
        cfg["run_name"] = f"{base_cfg['run_name']}_sweep{i}"
        res = run_cv(cfg, data, verbose=verbose)
        log_run(res, log_path, note="sweep " + json.dumps(dict(zip(keys, vals))))
        row = dict(zip(keys, vals))
        row.update({f"auc_{k}": v for k, v in res["summary"].items()})
        rows.append(row)
        print(f"[{i}/{len(combos)}] {dict(zip(keys, vals))} -> "
              + ", ".join(f"{k} {v:.4f}" for k, v in res["summary"].items()))
    score = f"auc_{'ensemble'}"
    return pd.DataFrame(rows).sort_values(score, ascending=False).reset_index(drop=True)


def _flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, key + "."))
        else:
            out[key] = json.dumps(v) if isinstance(v, (list, tuple)) else v
    return out
