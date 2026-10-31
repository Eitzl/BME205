import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

# point in PC space chosen by eye from MNIST_PCA_2D.png, at the outer edge of
# the cluster of 1s and away from the other digits
ONE_COORD = np.array([[-1000.0, 450.0]])


def save_digit(vec, path, title):
    fig, ax = plt.subplots(figsize=(3, 3))
    # fixed 0-255 scale so reconstructions are shown on the same scale as real
    # images (reconstructed pixels can fall slightly outside that range)
    ax.imshow(vec.reshape(28, 28), cmap="gray", vmin=0, vmax=255)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# 1. load 6000 images (784 pixels each) and their digit labels
X = np.load("MNIST_X_subset.npy")
y = np.load("MNIST_y_subset.npy")

# 2. reduce 784 dimensions to 2
pca = PCA(n_components=2)
X_2d = pca.fit_transform(X)

# 3. scatter plot colored by digit label
fig, ax = plt.subplots(figsize=(8, 7))
for digit in range(10):
    mask = y == digit
    ax.scatter(X_2d[mask, 0], X_2d[mask, 1], s=4, alpha=0.6,
               color=plt.cm.tab10(digit), label=str(digit))
ax.scatter(*ONE_COORD[0], marker="x", s=150, color="black", linewidths=3)
ax.annotate("chosen '1' coord", ONE_COORD[0], xytext=(10, 10),
            textcoords="offset points", fontweight="bold")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("MNIST projected onto the first 2 principal components")
ax.legend(markerscale=3, title="Digit")
fig.tight_layout()
fig.savefig("MNIST_PCA_2D.png", dpi=150)
plt.close(fig)

# 4. reconstruct image 0 from its 2 PC coordinates
reduced_0 = pca.transform(X[[0]])
reconstructed_0 = pca.inverse_transform(reduced_0)
save_digit(X[0], "MNIST_original.png", f"Original (label {y[0]})")
save_digit(reconstructed_0, "MNIST_reconstructed_2PC.png", "Reconstructed from 2 PCs")

# 5. generate a new '1' from a chosen point in PC space
generated_1 = pca.inverse_transform(ONE_COORD)
save_digit(generated_1, "MNIST_reconstructed_1_from_coord.png",
           f"Generated from ({ONE_COORD[0, 0]:.0f}, {ONE_COORD[0, 1]:.0f})")
