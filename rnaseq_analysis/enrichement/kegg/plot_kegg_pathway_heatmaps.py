import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler


def plot_pathway_heatmaps(tpm_file, gene_list_file, output_prefix="heatmap_",
                          log_transform=True, normalize='zscore',
                          vmin=None, vmax=None, clip_percentile=(2, 98),
                          cluster_rows=True, cluster_cols=False,
                          output_normalized_data=False):
    """
    为每个显著通路绘制表达热图，支持多种标准化方式，避免变化不可见。

    参数:
        tpm_file: TPM表达矩阵 (Excel, index=gene, columns=sample)
        gene_list_file: CSV 包含 'Pathway', 'Gene_list' 列
        output_prefix: 输出文件前缀
        log_transform: bool, 是否先进行 log2(TPM+1) 转换
        normalize: 标准化方式, 'zscore'(按基因), 'minmax'(按基因缩放到[0,1]), 'none'
        vmin, vmax: 颜色条的最小最大值 (若为None则自动计算)
        clip_percentile: (low, high) 用于截断异常值，如 (2,98) 表示取2%-98%分位数作为颜色范围
        cluster_rows: 是否对基因聚类
        cluster_cols: 是否对样本聚类
        output_normalized_data: 是否将每个通路的标准化数据保存为CSV
    """
    # 读取数据
    tpm_df = pd.read_excel(tpm_file, index_col=0)
    print(f"TPM矩阵维度: {tpm_df.shape}")

    genes_df = pd.read_csv(gene_list_file)
    if not all(col in genes_df.columns for col in ['Pathway', 'Gene_list']):
        raise ValueError("基因列表文件需要包含 'Pathway' 和 'Gene_list' 列")

    for idx, row in genes_df.iterrows():
        pathway = row['Pathway']
        genes = [g.strip() for g in row['Gene_list'].split(';') if g.strip()]

        # 提取可用基因
        present = [g for g in genes if g in tpm_df.index]
        if not present:
            print(f"通路 {pathway} 无可用基因，跳过")
            continue

        expr = tpm_df.loc[present].copy()

        # 1. log转换
        if log_transform:
            expr = np.log2(expr + 1)
            title_note = " (log2 TPM+1)"
        else:
            title_note = " (TPM)"

        # 2. 标准化
        if normalize == 'zscore':
            scaler = StandardScaler()
            data_plot = pd.DataFrame(scaler.fit_transform(expr.T).T,
                                     index=expr.index, columns=expr.columns)
            center = 0
            cmap = 'coolwarm'
            title_note += " | Z-score"
        elif normalize == 'minmax':
            scaler = MinMaxScaler()
            data_plot = pd.DataFrame(scaler.fit_transform(expr.T).T,
                                     index=expr.index, columns=expr.columns)
            center = None
            cmap = 'viridis'
            title_note += " | MinMax (0-1)"
        else:  # none
            data_plot = expr
            center = None
            cmap = 'viridis'
            title_note += " | Raw"

        # 3. 处理颜色范围
        if vmin is None or vmax is None:
            # 使用分位数截断
            low, high = clip_percentile
            vmin_clip = np.percentile(data_plot.values, low)
            vmax_clip = np.percentile(data_plot.values, high)
        else:
            vmin_clip, vmax_clip = vmin, vmax

        # 4. 绘图尺寸
        n_genes, n_samples = data_plot.shape
        fig_w = max(8, n_samples * 0.6)
        fig_h = max(6, n_genes * 0.4)

        # 5. 聚类热图
        g = sns.clustermap(data_plot,
                           method='average',
                           metric='euclidean',
                           cmap=cmap,
                           center=center,
                           vmin=vmin_clip,
                           vmax=vmax_clip,
                           col_cluster=cluster_cols,
                           row_cluster=cluster_rows,
                           figsize=(fig_w, fig_h),
                           dendrogram_ratio=0.1,
                           cbar_pos=(0.02, 0.85, 0.03, 0.1),
                           xticklabels=True,
                           yticklabels=True)
        plt.suptitle(f'{pathway}\n{len(present)} genes{title_note}', y=1.02)

        # 保存图片
        safe_name = pathway.replace('/', '_').replace(' ', '_')
        out_fig = f"{output_prefix}{safe_name}.png"
        plt.savefig(out_fig, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"热图已保存: {out_fig}")

        # 可选：保存标准化后的数据
        if output_normalized_data:
            out_csv = f"{output_prefix}{safe_name}_normalized_data.csv"
            data_plot.to_csv(out_csv)
            print(f"标准化数据已保存: {out_csv}")

# 使用示例
if __name__ == '__main__':
    plot_pathway_heatmaps(
        tpm_file="../2-get_gene_tpm/H-WT_vs_H-OX_tpm.xlsx",
        gene_list_file="H-WT_vs_H-OX_kegg_analysis_significant_genes.csv",
        output_prefix="heatmap_",
        log_transform=True,  # 先用log缓解高表达基因的绝对差异
        normalize='zscore',  # 按基因标准化，突出变化方向
        clip_percentile=(5, 95),  # 截断异常值，使颜色更集中
        cluster_rows=True,
        cluster_cols=False,
        output_normalized_data=True  # 输出标准化后的数值供检查
    )