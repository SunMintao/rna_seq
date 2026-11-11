import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_volcano(excel_path, sheet_name=0, log2fc_col='log2FoldChange', padj_col='padj',
                 log2fc_threshold=1.0, padj_threshold=0.05, top_n=0,
                 output_path=None, figsize=(10, 8), dpi=300,
                 print_stats=True, show_stats_on_plot=True):
    """
    根据差异表达分析结果绘制火山图，并统计上调/下调基因数量。

    Parameters
    ----------
    excel_path : str
        Excel 文件路径。
    sheet_name : str or int, default 0
        工作表名称或索引。
    log2fc_col : str, default 'log2FoldChange'
        对数倍数变化列名。
    padj_col : str, default 'padj'
        校正后 p 值列名。
    log2fc_threshold : float, default 1.0
        log2 倍数变化的阈值，绝对值大于此值视为差异表达。
    padj_threshold : float, default 0.05
        校正后 p 值的显著性阈值。
    top_n : int, default 0
        若 >0，则在图中标注 padj 最小的前 N 个基因 ID（需数据中有 'ID' 列）。
    output_path : str or None, default None
        输出图片保存路径，若 None 则显示图片。
    figsize : tuple, default (10, 8)
        图片尺寸。
    dpi : int, default 300
        图片分辨率。
    print_stats : bool, default True
        是否在控制台打印统计信息。
    show_stats_on_plot : bool, default True
        是否在图上显示统计信息。

    Returns
    -------
    fig : matplotlib.figure.Figure
        生成的图形对象。
    ax : matplotlib.axes.Axes
        图形的坐标轴对象。
    stats : dict
        包含上下调统计信息的字典，键：'up_regulated', 'down_regulated', 'not_significant'
    """
    # 读取数据
    df = pd.read_excel(excel_path, sheet_name=sheet_name)

    # 去除 padj 缺失的行
    df = df.dropna(subset=[padj_col])

    # 计算 -log10(padj)
    df['neg_log_padj'] = -np.log10(df[padj_col])

    # 标记基因调控类型
    df['regulation'] = 'Not significant'
    up_mask = (df[log2fc_col] > log2fc_threshold) & (df[padj_col] < padj_threshold)
    down_mask = (df[log2fc_col] < -log2fc_threshold) & (df[padj_col] < padj_threshold)
    df.loc[up_mask, 'regulation'] = 'Up-regulated'
    df.loc[down_mask, 'regulation'] = 'Down-regulated'

    # 统计数量
    up_count = up_mask.sum()
    down_count = down_mask.sum()
    not_sig_count = (~(up_mask | down_mask)).sum()
    stats = {
        'up_regulated': int(up_count),
        'down_regulated': int(down_count),
        'not_significant': int(not_sig_count)
    }

    # 打印统计信息
    if print_stats:
        print("=" * 40)
        print("Volcano Plot Statistics")
        print("=" * 40)
        print(f"Up-regulated genes     : {up_count}")
        print(f"Down-regulated genes   : {down_count}")
        print(f"Not significant genes  : {not_sig_count}")
        print(f"Total genes            : {len(df)}")
        print("=" * 40)

    # 设置颜色
    colors = {'Up-regulated': 'red', 'Down-regulated': 'blue', 'Not significant': 'gray'}

    # 绘图
    sns.set_style('whitegrid')
    fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(data=df, x=log2fc_col, y='neg_log_padj', hue='regulation',
                    palette=colors, alpha=0.6, s=20, edgecolor='none', ax=ax)

    # 添加阈值线
    ax.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', linewidth=0.8)
    ax.axvline(x=log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
    ax.axvline(x=-log2fc_threshold, color='black', linestyle='--', linewidth=0.8)

    # 标注 top N 显著基因（按 padj 升序）
    if top_n > 0 and 'ID' in df.columns:
        top_genes = df.nsmallest(top_n, padj_col)
        for _, row in top_genes.iterrows():
            ax.annotate(row['ID'], xy=(row[log2fc_col], row['neg_log_padj']),
                        xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.7)

    # 在图上显示统计信息
    if show_stats_on_plot:
        stats_text = f"Up: {up_count}\nDown: {down_count}"
        # 放置在右上角，稍微偏移避免覆盖图例
        ax.text(0.95, 0.95, stats_text,
                transform=ax.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                fontsize=10)

    # 坐标轴标签和标题
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title('Volcano Plot')
    ax.legend(title='Regulation')

    plt.tight_layout()

    # 保存或显示
    if output_path:
        plt.savefig(output_path, dpi=dpi)
    else:
        plt.show()

    return fig, ax, stats

# 使用示例
if __name__ == '__main__':
    # 调用函数，例如：
    # fig, ax, stats = plot_volcano('results.xlsx', sheet_name='Sheet1', top_n=10, output_path='volcano.png')
    pass