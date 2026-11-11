'''
GO 富集结果聚类分析
基于基因重叠相似度的 Mean‑Shift 聚类，并按要求绘制各图为单独文件
所有图形使用 Arial 字体，网络视图中为每个点标注 GO term 文本
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.manifold import MDS
from collections import Counter
import re
import warnings

warnings.filterwarnings('ignore')

# 设置全局字体为 Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


class GOTermClusteringAnalyzer:
    """基于 GO 富集结果的聚类分析器"""

    def __init__(self, file_path, sheet_name='Sheet1'):
        self.file_path = file_path
        self.sheet_name = sheet_name
        self.df = None
        self.gene_sets = None
        self.similarity_matrix = None
        self.distance_matrix = None
        self.coords_2d = None
        self.cluster_labels = None
        self.cluster_centers = None
        self.results_df = None
        self.cluster_summary = None

    def load_data(self):
        """加载并预处理 GO 富集结果"""
        print(f"正在加载数据: {self.file_path}")
        self.df = pd.read_excel(self.file_path, sheet_name=self.sheet_name)

        def parse_genes(gene_string):
            if pd.isna(gene_string) or gene_string == '':
                return set()
            return set(gene_string.strip().split(';'))

        self.df['Gene_Set'] = self.df['All_index'].apply(parse_genes)
        self.df = self.df[self.df['Gene_Set'].apply(len) > 0].reset_index(drop=True)
        print(f"原始 GO term 数量: {len(self.df)}")
        return self

    def filter_terms(self):
        """根据 Adjusted P-value 和原始 P-value 筛选需要分析的 GO term。"""
        sig_mask = self.df['Adjusted P-value'] < 0.05
        sig_count = sig_mask.sum()
        if sig_count >= 10:
            keep_indices = self.df[sig_mask].index
            print(f"显著条目数 {sig_count} >= 10，保留所有显著条目")
        else:
            keep_df = self.df.nsmallest(10, 'P-value')
            keep_indices = keep_df.index
            print(f"显著条目数 {sig_count} < 10，按 P-value 取前 {len(keep_df)} 条")

        self.df = self.df.loc[keep_indices].reset_index(drop=True)
        print(f"筛选后保留 GO term 数量: {len(self.df)}")
        return self

    def calculate_similarity(self):
        """计算基因重叠 Jaccard 相似度矩阵"""
        print("正在计算基因重叠相似度...")
        n = len(self.df)
        self.similarity_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):
                set_i = self.df.iloc[i]['Gene_Set']
                set_j = self.df.iloc[j]['Gene_Set']
                if i == j:
                    sim = 1.0
                else:
                    intersection = len(set_i & set_j)
                    union = len(set_i | set_j)
                    sim = intersection / union if union > 0 else 0.0
                self.similarity_matrix[i, j] = sim
                self.similarity_matrix[j, i] = sim

        self.distance_matrix = 1 - self.similarity_matrix
        print("相似度计算完成")
        return self

    def perform_clustering(self, bandwidth_quantile=0.2, min_bandwidth=0.3):
        """执行 Mean‑Shift 聚类"""
        print("正在执行聚类...")
        mds = MDS(n_components=2, dissimilarity='precomputed',
                  random_state=42, max_iter=3000, eps=1e-9)
        self.coords_2d = mds.fit_transform(self.distance_matrix)

        bandwidth = estimate_bandwidth(self.coords_2d,
                                       quantile=bandwidth_quantile,
                                       n_samples=len(self.df))
        print(f"估计带宽: {bandwidth}")
        if bandwidth < min_bandwidth:
            bandwidth = min_bandwidth
            print(f"使用最小带宽: {bandwidth}")

        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=True)
        self.cluster_labels = ms.fit_predict(self.coords_2d)
        self.cluster_centers = ms.cluster_centers_

        self.results_df = self.df.copy()
        self.results_df['Cluster'] = self.cluster_labels
        self.results_df['X_coord'] = self.coords_2d[:, 0]
        self.results_df['Y_coord'] = self.coords_2d[:, 1]
        self.results_df['Bubble_Size'] = self.results_df['Count_all'] * 10

        print(f"聚类完成，共 {len(np.unique(self.cluster_labels))} 个类")
        return self

    def analyze_clusters(self):
        """生成聚类统计摘要"""
        print("正在分析聚类...")
        cluster_counts = self.results_df.groupby('Cluster').size()
        intra_sim = []
        for cl in np.unique(self.cluster_labels):
            idx = self.results_df[self.results_df['Cluster'] == cl].index.tolist()
            if len(idx) > 1:
                sub = self.similarity_matrix[np.ix_(idx, idx)]
                avg = np.mean(sub[np.triu_indices_from(sub, k=1)])
            else:
                avg = 1.0
            intra_sim.append(avg)

        self.cluster_summary = pd.DataFrame({
            'Cluster': np.unique(self.cluster_labels),
            'Size': cluster_counts.values,
            'Avg_Intra_Similarity': intra_sim
        })
        self.cluster_summary['Mean_Adjusted_P'] = self.results_df.groupby('Cluster')['Adjusted P-value'].mean().values
        print("聚类分析完成")
        return self

    # --------------------------- 绘图函数（全部独立） ---------------------------
    def plot_cluster_scatter(self, save_path='cluster_scatter.png'):
        """聚类散点图"""
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(self.coords_2d[:, 0], self.coords_2d[:, 1],
                             c=self.cluster_labels, cmap='tab20', s=100,
                             alpha=0.8, edgecolors='black')
        if self.cluster_centers is not None and len(self.cluster_centers) > 0:
            ax.scatter(self.cluster_centers[:, 0], self.cluster_centers[:, 1],
                       marker='x', s=200, linewidths=3, color='red', label='Cluster centers')
        ax.set_xlabel('MDS Dimension 1', fontname='Arial')
        ax.set_ylabel('MDS Dimension 2', fontname='Arial')
        ax.set_title('Mean‑Shift Clustering of GO Terms (based on gene overlap)', fontname='Arial')
        ax.legend(prop={'family': 'Arial'})
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"散点图已保存: {save_path}")

    def plot_cluster_distribution(self, save_path='cluster_distribution.png'):
        """聚类条目分布柱状图"""
        cluster_counts = self.results_df['Cluster'].value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(12, 6))
        bars = ax.bar(cluster_counts.index.astype(str), cluster_counts.values,
                      color=plt.cm.tab20(np.linspace(0, 1, len(cluster_counts))))
        for bar, cnt in zip(bars, cluster_counts.values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(cnt), ha='center', va='bottom', fontname='Arial')
        ax.set_xlabel('Cluster ID', fontname='Arial')
        ax.set_ylabel('Number of GO terms', fontname='Arial')
        ax.set_title('Distribution of GO terms per cluster', fontname='Arial')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"分布图已保存: {save_path}")

    def plot_similarity_heatmap(self, save_path='similarity_heatmap.png'):
        """相似度矩阵热图（按聚类排序）"""
        sorted_idx = np.argsort(self.cluster_labels)
        sorted_sim = self.similarity_matrix[sorted_idx][:, sorted_idx]
        sorted_labels = self.cluster_labels[sorted_idx]

        fig, ax = plt.subplots(figsize=(12, 10))
        im = ax.imshow(sorted_sim, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
        uniq = np.unique(sorted_labels)
        pos = 0
        for cl in uniq:
            cnt = np.sum(sorted_labels == cl)
            if cnt > 0:
                rect = plt.Rectangle((pos-0.5, pos-0.5), cnt, cnt,
                                     linewidth=2, edgecolor='blue', facecolor='none')
                ax.add_patch(rect)
                pos += cnt
        ax.set_xlabel('GO term index (sorted by cluster)', fontname='Arial')
        ax.set_ylabel('GO term index (sorted by cluster)', fontname='Arial')
        ax.set_title('Gene overlap similarity matrix (Jaccard)', fontname='Arial')
        plt.colorbar(im, ax=ax, label='Similarity')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"热图已保存: {save_path}")

    def plot_cluster_keywords(self, save_path='cluster_keywords.png'):
        """聚类关键词可视化（水平文本列表）"""
        unique_clusters = np.unique(self.cluster_labels)
        cluster_keywords = {}
        for cl in unique_clusters:
            desc_list = self.results_df[self.results_df['Cluster'] == cl]['Description'].tolist()
            words = []
            for desc in desc_list:
                tokens = re.findall(r'\b[a-zA-Z]{4,}\b', desc.lower())
                stop = {'activity', 'process', 'response', 'to', 'and', 'of', 'the', 'for', 'in', 'by', 'via', 'from', 'with'}
                words.extend([w for w in tokens if w not in stop])
            freq = Counter(words)
            cluster_keywords[cl] = freq.most_common(5)

        fig, ax = plt.subplots(figsize=(12, 0.5 * len(unique_clusters) + 2))
        ax.axis('off')
        y = 0.9
        ax.text(0.05, 1.0, 'Top keywords for each cluster:', fontsize=12, fontweight='bold', fontname='Arial', va='top')
        for cl in sorted(unique_clusters):
            kw_str = ', '.join([f'{kw}({c})' for kw, c in cluster_keywords.get(cl, [])])
            ax.text(0.05, y, f'Cluster {cl}: {kw_str}', fontsize=10, fontname='Arial', va='top')
            y -= 0.1
        ax.set_ylim(0, 1)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"关键词图已保存: {save_path}")

    def plot_bubble_chart(self, save_path='bubble_chart.png', x_cap=10, fig_width=5, height_per_term=0.4):
        """
        气泡图：确保所有气泡完整显示。
        """
        plot_df = self.results_df.copy()
        plot_df['neg_log_p'] = -np.log10(plot_df['Adjusted P-value'].clip(lower=1e-300))

        # 截断处理（超过 x_cap 的点压缩到 x_cap 并加标注）
        if x_cap is not None:
            plot_df['neg_log_p_capped'] = plot_df['neg_log_p'].clip(upper=x_cap)
            capped = plot_df['neg_log_p'] > x_cap
            if capped.any():
                print(f"  注意：有 {capped.sum()} 个点的 -log10(P) 超过 {x_cap}，将显示在 x={x_cap} 并标注实际值")
        else:
            plot_df['neg_log_p_capped'] = plot_df['neg_log_p']
            capped = pd.Series([False] * len(plot_df))

        # 截断长描述
        plot_df['Short_Desc'] = plot_df['Description'].apply(
            lambda s: s[:60] + '...' if len(s) > 63 else s
        )
        plot_df = plot_df.sort_values('Adjusted P-value')

        n = len(plot_df)
        if n == 0:
            print("没有可绘制的 GO term！")
            return

        # 计算图形高度
        fig_height = max(4, n * height_per_term + 1)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # 颜色按聚类分配
        clusters = np.unique(plot_df['Cluster'])
        colors = plt.cm.tab20(np.linspace(0, 1, len(clusters)))
        color_dict = {cl: colors[i] for i, cl in enumerate(clusters)}

        # 绘制气泡
        for i, (_, row) in enumerate(plot_df.iterrows()):
            x = row['neg_log_p_capped']
            size = row.get('Bubble_Size', 10)  # 如果不存在则默认为 10
            # 确保气泡最小可见
            size = max(size * 3, 5)
            ax.scatter(x, i,
                       s=size,
                       alpha=0.7,
                       color=color_dict[row['Cluster']],
                       edgecolors='black',
                       linewidth=0.5)
            if capped.iloc[i]:
                ax.annotate(f'({row["neg_log_p"]:.1f})', (x, i),
                            xytext=(5, 0), textcoords='offset points',
                            fontsize=7, ha='left', va='center')

        # Y轴设置
        ax.set_ylim(-0.6, n - 0.4)  # 上下留出足够边距，确保首尾标签完整
        ax.set_yticks(range(n))
        ax.set_yticklabels(plot_df['Short_Desc'])

        # X轴设置：范围 = [最小值 - 5%跨度, 最大值 + 5%跨度]
        x_vals = plot_df['neg_log_p_capped'].values
        x_min, x_max = x_vals.min(), x_vals.max()
        if x_min == x_max:
            x_min -= 0.1
            x_max += 0.1
        x_range = x_max - x_min
        left = x_min - 0.05 * x_range
        right = x_max + 0.05 * x_range
        # 避免左负值导致看不到（如果数据都 >0 则正常）
        if left < 0 and x_min >= 0:
            left = 0
        ax.set_xlim(left, right)

        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.set_title('GO term bubble chart (size = gene count)')
        ax.grid(True, axis='x', alpha=0.3)

        # 图例
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color_dict[cl], label=f'Cluster {cl}') for cl in clusters]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout(pad=0.8)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"气泡图保存至 {save_path}")
    def plot_network_view(self, similarity_threshold=0.3, save_path='network_view.png'):
        """网络视图：高相似度的 GO term 之间画连线，并为每个点标注 GO term 文本"""
        fig, ax = plt.subplots(figsize=(16, 12))
        unique_clusters = np.unique(self.cluster_labels)
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))

        # 为每个点添加文本标签（使用 GO ID 或截断描述）
        labels = []
        for i, row in self.results_df.iterrows():
            # 使用 GO ID 作为标签，如果希望用简短的描述可以改成 row['Description'][:20]
            lbl = row['Description']
            labels.append(lbl)

        # 绘制聚类点
        for cl, col in zip(unique_clusters, colors):
            idx = self.cluster_labels == cl
            ax.scatter(self.coords_2d[idx, 0], self.coords_2d[idx, 1],
                       s=120, color=col, label=f'Cluster {cl}', alpha=0.8, edgecolors='black')

        # 绘制连接线
        n = len(self.df)
        for i in range(n):
            for j in range(i+1, n):
                if self.similarity_matrix[i, j] > similarity_threshold:
                    ax.plot([self.coords_2d[i, 0], self.coords_2d[j, 0]],
                            [self.coords_2d[i, 1], self.coords_2d[j, 1]],
                            'gray', alpha=0.2, linewidth=0.5)

        # 为每个点添加文本标注，避免重叠可适当调整位置
        for i, (x, y) in enumerate(self.coords_2d):
            ax.annotate(labels[i], (x, y), fontsize=8, fontname='Arial',
                        ha='center', va='bottom', alpha=0.7,
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.6))

        ax.set_xlabel('MDS Dimension 1', fontname='Arial')
        ax.set_ylabel('MDS Dimension 2', fontname='Arial')
        ax.set_title(f'Network view (similarity > {similarity_threshold})', fontname='Arial')
        ax.legend(prop={'family': 'Arial'}, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"网络视图已保存: {save_path}")

    def save_results(self, output_path='clustering_results.xlsx'):
        """保存聚类结果到 Excel"""
        out_df = self.results_df[['GO', 'Description', 'P-value', 'Adjusted P-value',
                                  'Cluster', 'X_coord', 'Y_coord', 'Count_all',
                                  'Gene_Set', 'All_index']].copy()
        out_df.to_excel(output_path, index=False)
        print(f"聚类结果保存至: {output_path}")
        return self



if __name__ == "__main__":
    input_file = "go_enrichment_result.xlsx"  # 请根据实际文件名修改
    output_prefix = "go_clustering"

    analyzer = GOTermClusteringAnalyzer(input_file)
    analyzer.load_data()
    analyzer.filter_terms()
    analyzer.calculate_similarity()
    analyzer.perform_clustering(bandwidth_quantile=0.2, min_bandwidth=0.3)
    analyzer.analyze_clusters()

    analyzer.plot_cluster_scatter(save_path=f"{output_prefix}_scatter.png")
    analyzer.plot_cluster_distribution(save_path=f"{output_prefix}_distribution.png")
    analyzer.plot_similarity_heatmap(save_path=f"{output_prefix}_heatmap.png")
    analyzer.plot_cluster_keywords(save_path=f"{output_prefix}_keywords.png")
    analyzer.plot_bubble_chart(save_path=f"{output_prefix}_bubble.png")
    analyzer.plot_network_view(similarity_threshold=0.3, save_path=f"{output_prefix}_network.png")
    analyzer.save_results(output_path=f"{output_prefix}_results.xlsx")

    print("全部图形已生成完毕！")
