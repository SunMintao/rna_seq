import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import warnings

warnings.filterwarnings('ignore')

# 设置中文字体（如果系统有，无则默认）
plt.rcParams['font.family'] = 'Arial'


def plot_cluster_heatmaps(
        cluster_file="go_cluster_results.xlsx",
        tpm_file="gene_tpm_value.xlsx",
        output_prefix="cluster_heatmap",
        cluster_ids=None,
        z_score='row',  # 按行标准化（基因）
        figsize=(12, 10),
        cmap='viridis',
        col_cluster=True,  # 样本是否聚类
        row_cluster=True  # 基因是否聚类
):
    """
    根据 GO 聚类结果和表达量数据，为每个 cluster 绘制热图。

    参数:
        cluster_file: 包含 GO term 及其 Cluster 和 All_index 的 Excel 文件
        tpm_file: 基因表达量 Excel 文件（行=基因，列=样本）
        output_prefix: 输出图片文件名前缀
        cluster_ids: 需要绘制的 cluster 列表（默认全部）
        z_score: 标准化方式，'row' 对基因标准化，None 则不标准化
        figsize: 每张热图的尺寸
        cmap: 颜色映射
        col_cluster: 是否对样本列进行聚类
        row_cluster: 是否对基因行进行聚类
    """
    # 1. 读取聚类结果
    df_cluster = pd.read_excel(cluster_file)
    # 确保必要列存在
    if 'Cluster' not in df_cluster.columns or 'All_index' not in df_cluster.columns:
        raise ValueError("聚类文件必须包含 'Cluster' 和 'All_index' 列")

    # 解析 All_index 列（分号分隔的基因 ID）
    df_cluster['Gene_list'] = df_cluster['All_index'].fillna('').apply(
        lambda x: [g.strip() for g in x.split(';') if g.strip()]
    )

    # 获取每个 cluster 包含的所有基因（并集）
    cluster_genes = {}
    for cl, group in df_cluster.groupby('Cluster'):
        all_genes = set()
        for genes in group['Gene_list']:
            all_genes.update(genes)
        cluster_genes[cl] = list(all_genes)

    # 2. 读取表达量数据
    df_tpm = pd.read_excel(tpm_file, index_col=0)  # 第一列为基因 ID
    # 检查基因名是否匹配（全部转为字符串）
    df_tpm.index = df_tpm.index.astype(str).str.strip()

    # 3. 为每个 cluster 绘制热图
    if cluster_ids is None:
        cluster_ids = sorted(cluster_genes.keys())

    for cl in cluster_ids:
        genes = cluster_genes.get(cl, [])
        if len(genes) == 0:
            print(f"Cluster {cl} 没有基因，跳过")
            continue

        # 提取在表达矩阵中存在的基因
        exist_genes = [g for g in genes if g in df_tpm.index]
        missing = set(genes) - set(exist_genes)
        if missing:
            print(f"Cluster {cl}: 有 {len(missing)} 个基因不在表达文件中，已跳过")

        if len(exist_genes) == 0:
            print(f"Cluster {cl}: 没有可用的基因，跳过")
            continue

        # 获取子表达矩阵
        sub_df = df_tpm.loc[exist_genes].copy()

        # 可选：对基因行进行 log2 转换（避免极值影响，可根据需要开启）
        # sub_df = np.log2(sub_df + 1)

        # 跳过低表达基因（表达值全部极低），可选
        # sub_df = sub_df[sub_df.max(axis=1) > 1]   # 示例

        # 标准化（z-score）
        if z_score == 'row':
            sub_df = sub_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
            center = 0
            vmin, vmax = -2, 2  # 可调整
        else:
            center = None
            vmin, vmax = None, None

        # 绘图
        plt.figure(figsize=figsize)
        ax = sns.clustermap(
            sub_df,
            method='ward',  # 聚类方法
            metric='euclidean',  # 距离度量
            z_score=z_score,  # 已经手动标准化，这里设为 None
            col_cluster=col_cluster,
            row_cluster=row_cluster,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={'label': 'Z-score'},
            figsize=figsize
        )
        # 调整布局
        ax.ax_heatmap.set_xlabel('Samples')
        ax.ax_heatmap.set_ylabel('Genes')
        ax.ax_heatmap.set_title(f'Cluster {cl} expression pattern (n={len(exist_genes)} genes)')
        plt.tight_layout()
        # 保存图片
        output_file = f"{output_prefix}_cluster_{cl}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"已保存: {output_file} (包含 {len(exist_genes)} 个基因)")


# 运行示例
if __name__ == "__main__":
    plot_cluster_heatmaps(
        cluster_file="H-WT_vs_H-lht1_go_cluster_results.xlsx",
        tpm_file="gene_tpm_value.xlsx",
        output_prefix="cluster_heatmap",
        cluster_ids=[0, 1, 2, 3],  # 可以根据实际 cluster 数量修改，或注释掉全部绘制
        z_score='row',
        figsize=(14, 10),
        cmap='RdBu_r',
        col_cluster=True,
        row_cluster=True
    )