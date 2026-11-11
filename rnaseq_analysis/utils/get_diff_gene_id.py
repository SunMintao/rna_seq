import pandas as pd

def get_diff_genes(
    file_path: str,
    log2fc_threshold: float = 1.0,
    qvalue_threshold: float = 0.05,
    output_file: str = None
):
    """
    根据 log2FoldChange 和 padj 筛选差异表达基因，并可选择将结果写入 Excel 文件。

    参数
    ----------
    file_path : str
        包含差异分析结果的 Excel 文件路径（必须包含 ID、log2FoldChange、padj 列）。
    log2fc_threshold : float, optional
        log2FoldChange 的绝对值阈值，默认 1.0。
    qvalue_threshold : float, optional
        校正后 p 值（padj）的阈值，默认 0.05。
    output_file : str, optional
        如果提供，则将所有差异基因（含上下调方向）写入该 Excel 文件（.xlsx）。
        文件包含两列：GeneID 和 Regulation（Up/Down）。

    返回
    -------
    up_regulated : list
        上调基因的 ID 列表（log2FoldChange > log2fc_threshold 且 padj < qvalue_threshold）。
    down_regulated : list
        下调基因的 ID 列表（log2FoldChange < -log2fc_threshold 且 padj < qvalue_threshold）。
    """
    # 读取 Excel 文件
    df = pd.read_excel(file_path)

    # 确保必要的列存在
    required_cols = ['ID', 'log2FoldChange', 'padj']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"文件中缺少必要列: {col}")

    # 去除 padj 为缺失值的行
    df = df.dropna(subset=['padj'])

    # 筛选满足 padj 阈值条件的基因
    sig_df = df[df['padj'] < qvalue_threshold]

    # 根据 log2FoldChange 分别提取上下调
    up = sig_df[sig_df['log2FoldChange'] > log2fc_threshold]['ID'].tolist()
    down = sig_df[sig_df['log2FoldChange'] < -log2fc_threshold]['ID'].tolist()

    # 如果指定了输出文件，写入所有差异基因（含方向）到 Excel
    if output_file is not None:
        # 确保输出文件名以 .xlsx 结尾（如果没有则添加）
        if not output_file.endswith('.xlsx'):
            output_file += '.xlsx'

        # 构建结果数据框
        result_data = []
        for gene in up:
            result_data.append({'GeneID': gene, 'Regulation': 'Up'})
        for gene in down:
            result_data.append({'GeneID': gene, 'Regulation': 'Down'})
        result_df = pd.DataFrame(result_data)

        # 写入 Excel 文件
        result_df.to_excel(output_file, index=False)

    return up, down

if __name__ == '__main__':
    # 仅获取列表
    up_ids, down_ids = get_diff_genes('results.xlsx')
    # 获取列表并同时写入 Excel 文件
    up_ids, down_ids = get_diff_genes('results.xlsx', output_file='diff_genes.xlsx')