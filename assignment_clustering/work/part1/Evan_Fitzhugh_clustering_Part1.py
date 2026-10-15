import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def kmeans(X, k, rng, max_iter=300, tol=1e-4):
    """K-Means with random-point initialization and Euclidean distance."""
    centroids = X[rng.choice(len(X), size=k, replace=False)].astype(np.float64)
    X = X.astype(np.float64)
    x_sq = np.sum(X**2, axis=1)[:, None]

    for _ in range(max_iter):
        # squared Euclidean distances via |x|^2 - 2x.c + |c|^2
        d = x_sq - 2 * X @ centroids.T + np.sum(centroids**2, axis=1)[None, :]
        labels = np.argmin(d, axis=1)

        new = centroids.copy()
        for j in range(k):
            members = X[labels == j]
            if len(members):  # empty cluster: keep previous centroid
                new[j] = members.mean(axis=0)

        shift = np.max(np.linalg.norm(new - centroids, axis=1))
        centroids = new
        if shift < tol:
            break

    # final assignment consistent with the final centroids
    d = x_sq - 2 * X @ centroids.T + np.sum(centroids**2, axis=1)[None, :]
    return centroids, np.argmin(d, axis=1)


def clustering_error(labels, y):
    """Sum over clusters of samples not matching the cluster's majority label."""
    error = 0
    for j in np.unique(labels):
        counts = np.bincount(y[labels == j])
        error += counts.sum() - counts.max()
    return int(error)


def plot_centroids(centroids, k, path):
    fig, axes = plt.subplots(1, k, figsize=(1.5 * k, 2), squeeze=False)
    for j, ax in enumerate(axes[0]):
        ax.imshow(centroids[j].reshape(28, 28), cmap="gray", vmin=0, vmax=255)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="K-Means on MNIST subset")
    parser.add_argument("-k", type=int, required=True, help="number of clusters")
    parser.add_argument("--seed", type=int, default=None, help="optional RNG seed")
    parser.add_argument("--data-dir", default=".", help="directory with the .npy files")
    args = parser.parse_args()

    X = np.load(f"{args.data_dir}/MNIST_X_subset.npy", allow_pickle=True)
    y = np.load(f"{args.data_dir}/MNIST_y_subset.npy", allow_pickle=True).astype(int)

    rng = np.random.default_rng(args.seed)
    centroids, labels = kmeans(X, args.k, rng)
    plot_centroids(centroids, args.k, f"centroids_k{args.k}.png")
    print(f"k={args.k}, ERROR={clustering_error(labels, y)}")


if __name__ == "__main__":
    main()
