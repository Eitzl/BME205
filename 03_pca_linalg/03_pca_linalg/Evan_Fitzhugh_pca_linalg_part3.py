import numpy as np
import pandas as pd
from sklearn.decomposition import NMF

N_COMPONENTS = 5

# 1. load 1355 dogs x 784 SNPs (allele counts 0/1/2) and each dog's clade
X = np.load("dogs_X.npy")
clades = np.load("dogs_clades.npy", allow_pickle=True)

# 2. factor X ~ W @ H; W (dogs x components) holds each dog's ancestry weights.
# fixed random_state for reproducibility; default max_iter=200 does not converge
nmf = NMF(n_components=N_COMPONENTS, random_state=0, max_iter=1000)
W = nmf.fit_transform(X)

# normalize each dog's weights to proportions that sum to 1
W_frac = W / W.sum(axis=1, keepdims=True)

# 3. mean ancestry fractions per clade; groupby sorts clades alphabetically
# ('**' names sort before letters)
summary = pd.DataFrame(W_frac).groupby(clades).mean()

# order components by prevalence across clades (each clade weighted equally)
order = summary.mean(axis=0).sort_values(ascending=False).index
summary = summary[order]
summary.columns = [f"ancestry{i}" for i in range(1, N_COMPONENTS + 1)]
summary.index.name = "clade"

summary.round(2).to_csv("dogs_ancestry_summary.tsv", sep="\t", float_format="%.2f")
