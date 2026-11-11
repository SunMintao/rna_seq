"""
pca.py

Principal component analysis (PCA) implemented from scratch with NumPy.
Uses SVD on the centered data matrix.
"""

import numpy as np


class PCA:
    """
    PCA via singular value decomposition (SVD).

    Steps:
    1. Center each feature (column) by subtracting its mean.
    2. Apply SVD: X_centered = U @ diag(S) @ Vt
    3. Principal components are rows of Vt; scores = X_centered @ Vt.T

    Parameters
    ----------
    n_components : int, optional
        Number of principal components to retain. If None, keep all.
    """

    def __init__(self, n_components=None):
        self.n_components = n_components
        self.components_ = None
        self.explained_variance_ = None
        self.explained_variance_ratio_ = None
        self.singular_values_ = None
        self.mean_ = None
        self.n_samples_ = None
        self.n_features_ = None

    def fit(self, X):
        """
        Fit PCA on training data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input matrix (rows = samples, columns = features).

        Returns
        -------
        self
        """
        X = np.asarray(X, dtype=float)
        if X.ndim != 2:
            raise ValueError("X must be a 2D array")

        self.n_samples_, self.n_features_ = X.shape
        if self.n_samples_ < 2:
            raise ValueError("PCA requires at least 2 samples")

        self.mean_ = X.mean(axis=0)
        X_centered = X - self.mean_

        # X_centered = U @ diag(S) @ Vt
        _, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

        n_components = self.n_components
        if n_components is None:
            n_components = min(self.n_samples_, self.n_features_)
        n_components = min(n_components, len(S))

        self.components_ = Vt[:n_components]
        self.singular_values_ = S[:n_components]

        # Each PC variance = s^2 / (n_samples - 1)
        self.explained_variance_ = (S[:n_components] ** 2) / (self.n_samples_ - 1)
        total_var = np.sum(S ** 2) / (self.n_samples_ - 1)
        self.explained_variance_ratio_ = self.explained_variance_ / total_var

        return self

    def transform(self, X):
        """
        Project data onto principal components.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        Returns
        -------
        ndarray of shape (n_samples, n_components)
            Principal component scores.
        """
        if self.components_ is None:
            raise RuntimeError("Call fit() before transform()")

        X = np.asarray(X, dtype=float)
        X_centered = X - self.mean_
        return X_centered @ self.components_.T

    def fit_transform(self, X):
        """Fit PCA and return transformed scores."""
        self.fit(X)
        return self.transform(X)

    def inverse_transform(self, scores):
        """
        Reconstruct original space from PC scores.

        Parameters
        ----------
        scores : array-like of shape (n_samples, n_components)

        Returns
        -------
        ndarray of shape (n_samples, n_features)
        """
        if self.components_ is None:
            raise RuntimeError("Call fit() before inverse_transform()")

        scores = np.asarray(scores, dtype=float)
        return scores @ self.components_ + self.mean_


def pca_from_covariance(X, n_components=None):
    """
    Alternative PCA via covariance matrix eigendecomposition.

    Useful for understanding the math: eigenvectors of Cov(X) are PCs.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
    n_components : int, optional

    Returns
    -------
    scores : ndarray of shape (n_samples, n_components)
    components : ndarray of shape (n_components, n_features)
    explained_variance_ratio : ndarray of shape (n_components,)
    """
    X = np.asarray(X, dtype=float)
    n_samples, n_features = X.shape
    X_centered = X - X.mean(axis=0)

    # Covariance matrix (n_features x n_features)
    cov = (X_centered.T @ X_centered) / (n_samples - 1)

    # Eigen-decomposition; eigenvalues are in ascending order
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    if n_components is None:
        n_components = min(n_samples, n_features)
    n_components = min(n_components, len(eigenvalues))

    components = eigenvectors[:, :n_components].T
    explained_variance = eigenvalues[:n_components]
    explained_variance_ratio = explained_variance / eigenvalues.sum()
    scores = X_centered @ components.T

    return scores, components, explained_variance_ratio


if __name__ == "__main__":
    # Demo with synthetic 2D data
    np.random.seed(42)
    n = 100
    x1 = np.random.normal(0, 2, n)
    x2 = 0.8 * x1 + np.random.normal(0, 0.5, n)
    X = np.column_stack([x1, x2])

    pca = PCA(n_components=2)
    scores = pca.fit_transform(X)

    print("=== SVD-based PCA ===")
    print(f"Explained variance ratio: {pca.explained_variance_ratio_}")
    print(f"PC scores shape: {scores.shape}")
    print(f"Reconstruction error (MSE): {np.mean((pca.inverse_transform(scores) - X) ** 2):.6f}")

    print("\n=== Covariance-based PCA ===")
    scores2, comps2, ratio2 = pca_from_covariance(X, n_components=2)
    print(f"Explained variance ratio: {ratio2}")
    print(f"Components match SVD method: {np.allclose(pca.components_, comps2)}")
