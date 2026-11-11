def filter_pydeseq2_data(counts_df, metadata, min_gene_total=10):
    """
    过滤表达矩阵和元数据：
    1. 保留 metadata['condition'] 非空的样本。
    2. 保留表达矩阵中所有样本总 counts >= min_gene_total 的基因。

    Parameters
    ----------
    counts_df : pd.DataFrame
        表达矩阵，行名为样本，列名为基因，值为整数 counts。
    metadata : pd.DataFrame
        样本元数据，行名为样本，至少包含一列 'condition'。
    min_gene_total : int, default=10
        基因保留的最小总表达量（所有样本之和）。

    Returns
    -------
    counts_df_filtered : pd.DataFrame
        过滤后的表达矩阵。
    metadata_filtered : pd.DataFrame
        过滤后的元数据（仅包含保留的样本）。
    """
    # 1. 找出 condition 列非空的样本（布尔索引）
    samples_to_keep = ~metadata['condition'].isna()

    # 2. 根据样本过滤表达矩阵和元数据
    counts_df_filtered = counts_df.loc[samples_to_keep]
    metadata_filtered = metadata.loc[samples_to_keep]

    # 3. 保留所有样本总 counts >= min_gene_total 的基因
    genes_to_keep = counts_df_filtered.columns[counts_df_filtered.sum(axis=0) >= min_gene_total]
    counts_df_filtered = counts_df_filtered[genes_to_keep]

    # 可选：打印过滤后的形状信息
    print(f"过滤后样本数: {len(metadata_filtered)}")
    print(f"过滤后基因数: {len(genes_to_keep)}")

    return counts_df_filtered, metadata_filtered

# 当该脚本作为主程序运行时（而不是被导入为模块）
if __name__ == "__main__":
    # 假设已经通过 load_pydeseq2_data 获得了 counts 和 meta
    counts, meta = load_pydeseq2_data(...)

    # 调用过滤函数
    counts_filt, meta_filt = filter_pydeseq2_data(counts, meta, min_gene_total=10)

    # 继续后续 pydeseq2 分析