# Breast cancer prediction with PRS (Kaggle, BME 205)

Predict `breast_cancer` from 5,460 polygenic risk score (PRS) Z-scores plus an
ancestry column. Scored by AUC-ROC.

## Setup
1. Put the Kaggle files in `data/` here (or in a `data/` folder at the repo root):
   `data/train.csv`, `data/test.csv`. Both locations are git-ignored, because the
   competition data must not be redistributed.
2. Install and open:
   ```bash
   pixi install
   pixi run lab
   ```

## Files
- `prs_model.ipynb`: the notebook. All tunable settings are in the **KNOBS** cell (`CFG`).
- `prs_lib.py`: the machinery (loading, in-fold preprocessing, models, CV, sweeps, logging).
- `runs_log.csv` (created on first run): one row per run with every knob and its scores.
- `submissions/` (created on first run): CSVs to upload to Kaggle.

## Workflow
1. Run once with the defaults to get a baseline.
2. Change one or two knobs, Run All, and compare in the last cell.
3. Use the sweep cell to grid-search knobs (e.g. `"models.logreg.C": [0.01, 0.1]`).
4. Always compare against `best_single` (the strongest single PRS). A model that
   can't beat it is overfitting.

Before submitting, check that `CFG["output"]["id_col_name"]` and the row order
match Kaggle's sample submission file.
