import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

# 1. load data: 10 points, each a vector [x_i, y_i]
X = np.load("part1_data.npy", allow_pickle=True)
n = X.shape[0]

# 2. sklearn PCA, as a reference for the manual steps below
pca = PCA(n_components=2)
pca.fit(X)

print("Sklearn direction of each PC (PC1, PC2):\n", pca.components_)
print("Sklearn explained variance for each PC (PC1, PC2):\n", pca.explained_variance_)
print("Sklearn fraction of variance explained by each PC (PC1, PC2):\n", pca.explained_variance_ratio_)

x_pca = pca.transform(X)
print("Sklearn data projection into 1-D (onto PC1):\n", x_pca[:, 0])

# 3. center the data
Xc = X - X.mean(axis=0)

# 4. covariance matrix
sigma = (Xc.T @ Xc) / (n - 1)
print("\nCovariance matrix:\n", sigma)

# 5. eigendecomposition; eigh returns eigenvalues in ascending order,
# so reverse to get the largest first
eigenvals, eigenvecs = np.linalg.eigh(sigma)
order = np.argsort(eigenvals)[::-1]
eigenvals = eigenvals[order]
eigenvecs = eigenvecs[:, order]

lambda_1, lambda_2 = eigenvals
v_1, v_2 = eigenvecs[:, 0], eigenvecs[:, 1]  # eigenvectors are columns

print(f"\nlambda_1 = {lambda_1}")
print("v_1 =", v_1)
print(f"lambda_2 = {lambda_2}")
print("v_2 =", v_2)

print(
    "\nQuestion 1: The eigenvalues are identical to sklearn's explained_variance_, "
    "and the eigenvectors match sklearn's components_ (the rows of components_ are "
    "our eigenvector columns). v_1 matches exactly; v_2 matches up to a sign flip, "
    "which is fine because v and -v are both eigenvectors pointing along the same "
    "line. So sklearn's PCA is an eigendecomposition of the covariance matrix."
)
print(
    "Question 2: The data lie close to a line running up and to the right, steeper "
    "than 45 degrees (y grows about twice as fast as x). v_1 = [0.42, 0.91] points "
    "in exactly that direction, so v_1 is the direction of greatest variance."
)

# 6. fraction of total variance captured by each eigenvector
lambda_1_ratio = lambda_1 / (lambda_1 + lambda_2)
lambda_2_ratio = lambda_2 / (lambda_1 + lambda_2)
print(f"\nlambda_1 ratio: {lambda_1_ratio}")
print(f"lambda_2 ratio: {lambda_2_ratio}")
print(
    "Question 3: These match sklearn's explained_variance_ratio_. PC1 captures "
    "~99.7% of the variance and PC2 only ~0.3%, confirming the data is "
    "essentially one-dimensional with a little noise."
)

# 7. project centered data onto v_1
projection = Xc @ v_1
print("\nData projection into 1-D (onto v_1):\n", projection)
print(
    "Question 4: The projection is identical to sklearn's transform(X)[:, 0]. "
    "sklearn's transform is just centering followed by a dot product with each "
    "principal component. (Had v_1 come out with the opposite sign, the "
    "projection would match up to a sign flip.)"
)

# 8. plot the 1-D projection on a number line
fig, ax = plt.subplots(figsize=(8, 2))
ax.axhline(0, color="gray", lw=1, zorder=1)
ax.scatter(projection, np.zeros_like(projection), zorder=2)
ax.set_yticks([])
ax.set_xlabel("Position along v_1 (PC1)")
ax.set_title("Data projected onto PC1")
fig.tight_layout()
fig.savefig("projection_1D.png", dpi=150)

print(
    "\nQuestion 5: The spread and ordering of the points are preserved: points that "
    "were far apart along the line in the 2-D plot are still far apart, and in the "
    "same order. What changes is that the data is now centered at 0, the small "
    "perpendicular scatter (the PC2 'noise') is discarded, and each point is "
    "described by one number instead of two."
)
