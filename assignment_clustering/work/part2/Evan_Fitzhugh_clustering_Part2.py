import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

K = 30
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def clustering_error(assignments, true_labels):
    """Sum over clusters of (cluster size - count of the majority true label)."""
    error = 0
    for c in np.unique(assignments):
        members = true_labels[assignments == c]
        error += len(members) - Counter(members).most_common(1)[0][1]
    return error


def main():
    X = np.load(os.path.join(DATA_DIR, "dogs_X.npy"), allow_pickle=True)
    clades = np.load(os.path.join(DATA_DIR, "dogs_clades.npy"), allow_pickle=True)

    Z = linkage(X, method="ward")

    # truncated dendrogram: 30 leaves, each labeled with its sample count
    plt.figure(figsize=(12, 6))
    dendrogram(Z, truncate_mode="lastp", p=K, show_leaf_counts=True)
    plt.title("Dogs SNP hierarchical clustering (Ward), truncated to 30 clusters")
    plt.xlabel("Cluster (number of samples)")
    plt.ylabel("Ward distance")
    plt.tight_layout()
    plt.savefig("Dogs_dendrogram_truncated.png", dpi=150)

    assignments = fcluster(Z, t=K, criterion="maxclust")
    print(f"k={K}, ERROR={clustering_error(assignments, clades)}")


if __name__ == "__main__":
    main()
