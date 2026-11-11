#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import os
from typing import Optional, Tuple

# 设置 Matplotlib 全局字体为 Arial（必须在任何绘图之前）
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans', 'Helvetica', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False  # 正确显示负号


def cluster_samples(
    data: pd.DataFrame,
    output_dir: str = "./clustering_results",
    method: str = "both",
    max_k: int = 10,
    threshold: Optional[float] = None,
    heatmap: bool = False,
    seed: int = 42,
    scale: bool = False,
    top_genes_heatmap: int = 500
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Perform sample clustering on VST-normalized transcriptomic expression matrix.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix with genes as rows and samples as columns.
        Index: gene IDs, Columns: sample names.
    output_dir : str, default="./clustering_results"
        Directory to save output results.
    method : str, default="both"
        Clustering method: 'hierarchical', 'kmeans', or 'both'.
    max_k : int, default=10
        Maximum number of clusters for K-means to try.
    threshold : float, optional
        Height threshold to cut the hierarchical dendrogram.
        If None, automatically set to 50% of the maximum distance.
    heatmap : bool, default=False
        Whether to generate a clustered heatmap (computationally heavy).
    seed : int, default=42
        Random seed for reproducibility.
    scale : bool, default=False
        Whether to Z-score standardize data before K-means (usually not needed for VST).
    top_genes_heatmap : int, default=500
        Number of top variable genes to use in the heatmap.

    Returns
    -------
    h_labels_df : pd.DataFrame or None
        Hierarchical clustering labels.
    k_labels_df : pd.DataFrame or None
        K-means clustering labels.
    """
    np.random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    # ---------- Data Preprocessing ----------
    X = data.T  # samples × genes
    sample_names = X.index.tolist()
    print(f"Data shape: {data.shape[0]} genes × {data.shape[1]} samples")
    print(f"Samples for clustering: {X.shape[0]}, Genes: {X.shape[1]}")

    if X.isnull().any().any():
        print("Warning: Missing values detected. Removing genes with any missing values.")
        X = X.dropna(axis=1)

    h_labels_df = None
    k_labels_df = None

    # ---------- Hierarchical Clustering ----------
    if method in ["hierarchical", "both"]:
        print("\n=== Hierarchical Clustering ===")
        dist_matrix = pdist(X, metric='euclidean')
        linkage_matrix = linkage(dist_matrix, method='ward')

        # Dendrogram
        plt.figure(figsize=(12, 6))
        dendrogram(linkage_matrix, labels=sample_names, leaf_rotation=90, leaf_font_size=8)
        plt.title("Sample Hierarchical Clustering Dendrogram (Ward's Linkage, Euclidean Distance)")
        plt.xlabel("Samples")
        plt.ylabel("Distance")
        plt.tight_layout()
        dendro_file = os.path.join(output_dir, "hierarchical_dendrogram.pdf")
        plt.savefig(dendro_file, dpi=300)
        plt.close()
        print(f"Dendrogram saved to: {dendro_file}")

        # Cut tree by threshold
        max_dist = linkage_matrix[-1, 2]
        if threshold is None:
            threshold = max_dist * 0.5
            print(f"Auto-selected threshold: {threshold:.2f} (50% of max distance {max_dist:.2f})")
        else:
            if threshold > max_dist:
                print(f"Warning: threshold {threshold} > max distance {max_dist:.2f}. Auto-setting to {max_dist*0.5:.2f}")
                threshold = max_dist * 0.5

        cluster_labels_h = fcluster(linkage_matrix, t=threshold, criterion='distance')
        n_clusters_h = len(set(cluster_labels_h))
        print(f"Hierarchical clustering yielded {n_clusters_h} clusters")

        h_labels_df = pd.DataFrame({
            "Sample": sample_names,
            "Hierarchical_Cluster": cluster_labels_h
        })
        h_label_file = os.path.join(output_dir, "hierarchical_labels.csv")
        h_labels_df.to_csv(h_label_file, index=False)
        print(f"Hierarchical labels saved to: {h_label_file}")

    # ---------- K-means Clustering ----------
    if method in ["kmeans", "both"]:
        print("\n=== K-means Clustering ===")
        X_scaled = X.values
        if scale:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_scaled)
            print("Data has been Z-score standardized.")

        max_k = min(max_k, X.shape[0] - 1)
        if max_k < 2:
            print("Insufficient samples for K-means (need at least 2).")
            return h_labels_df, k_labels_df

        inertias = []
        silhouette_scores = []
        K_range = range(2, max_k + 1)

        for k in K_range:
            kmeans = KMeans(n_clusters=k, n_init=10, random_state=seed)
            kmeans.fit(X_scaled)
            inertias.append(kmeans.inertia_)
            sil = silhouette_score(X_scaled, kmeans.labels_)
            silhouette_scores.append(sil)

        # Elbow and silhouette plots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.plot(K_range, inertias, 'bo-')
        ax1.set_xlabel("K")
        ax1.set_ylabel("Inertia (SSE)")
        ax1.set_title("Elbow Method")

        ax2.plot(K_range, silhouette_scores, 'ro-')
        ax2.set_xlabel("K")
        ax2.set_ylabel("Silhouette Coefficient")
        ax2.set_title("Silhouette Coefficient")

        plt.tight_layout()
        k_plot_file = os.path.join(output_dir, "kmeans_elbow_silhouette.pdf")
        plt.savefig(k_plot_file, dpi=300)
        plt.close()
        print(f"K-means evaluation plot saved to: {k_plot_file}")

        # Select best K (max silhouette)
        best_k = K_range[np.argmax(silhouette_scores)]
        print(f"Optimal K (max silhouette): {best_k}")

        # Final K-means
        kmeans_final = KMeans(n_clusters=best_k, n_init=10, random_state=seed)
        cluster_labels_k = kmeans_final.fit_predict(X_scaled) + 1  # 1-based indexing

        k_labels_df = pd.DataFrame({
            "Sample": sample_names,
            "Kmeans_Cluster": cluster_labels_k
        })
        k_label_file = os.path.join(output_dir, "kmeans_labels.csv")
        k_labels_df.to_csv(k_label_file, index=False)
        print(f"K-means labels saved to: {k_label_file}")

        # ---------- PCA Visualization ----------
        pca = PCA(n_components=2, random_state=seed)
        pca_result = pca.fit_transform(X_scaled)
        pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
        pca_df["Sample"] = sample_names
        pca_df["Cluster"] = cluster_labels_k

        plt.figure(figsize=(10, 8))
        sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Cluster", palette="Set1", s=100)
        for i, row in pca_df.iterrows():
            plt.text(row["PC1"]+0.02, row["PC2"]+0.02, row["Sample"], fontsize=8)
        plt.title(f"PCA Scatter Plot (K={best_k})")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} Variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} Variance)")
        plt.tight_layout()
        pca_file = os.path.join(output_dir, "pca_scatter.pdf")
        plt.savefig(pca_file, dpi=300)
        plt.close()
        print(f"PCA scatter plot saved to: {pca_file}")

    # ---------- Optional Heatmap ----------
    if heatmap:
        print("\nGenerating heatmap (may take a while)...")
        gene_var = data.var(axis=1).sort_values(ascending=False)
        top_genes = gene_var.head(top_genes_heatmap).index.tolist()
        data_sub = data.loc[top_genes, :]
        data_heat = data_sub.T  # samples × genes

        g = sns.clustermap(data_heat, method='ward', metric='euclidean',
                           z_score=0,  # z-score per gene
                           col_cluster=True, row_cluster=True,
                           figsize=(12, 10), cmap='vlag',
                           dendrogram_ratio=0.15, cbar_pos=(0.02, 0.8, 0.03, 0.18))
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=6)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)
        heat_file = os.path.join(output_dir, "clustermap.pdf")
        g.savefig(heat_file, dpi=300)
        plt.close()
        print(f"Clustered heatmap saved to: {heat_file}")

    print("\nAll analyses completed successfully!")
    return h_labels_df, k_labels_df


# ============================================================
# Usage Example
# ============================================================
if __name__ == "__main__":
    # Example:
    # df = pd.read_csv("vst_counts.csv", index_col=0)
    # h_labels, k_labels = cluster_samples(
    #     data=df,
    #     output_dir="./my_results",
    #     method="both",
    #     max_k=8,
    #     heatmap=True,
    #     scale=False
    # )
    pass