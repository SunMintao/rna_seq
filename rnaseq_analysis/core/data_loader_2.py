import os
import pandas as pd

def load_pydeseq2_data(counts_file, metadata_file, transpose_counts=True, output_path=None):
    """
    读取不同格式的表达矩阵和样本分组文件，返回适用于 pydeseq2 的数据框。

    Parameters
    ----------
    counts_file : str
        表达矩阵文件路径，支持 .xlsx, .xls, .csv, .txt（自动推断分隔符）。
    metadata_file : str
        样本分组信息文件路径，格式同上。sample(样本+重复的名字) condition(处理) group(分组，每一组包含2个处理)
    transpose_counts : bool, default=True
        是否需要转置表达矩阵（pydeseq2 要求行=样本，列=基因，原始文件常为行=基因）。
    output_path : str or None
        可选，若提供则创建该目录（用于保存后续结果）。

    Returns
    -------
    counts_df : pd.DataFrame
        转置后的表达矩阵，行名为样本ID，列名为基因名，值均为整数计数。
    metadata : pd.DataFrame
        样本分组表，行名为样本ID，至少包含分组信息列（如 'condition'）。
    """
    # 若指定输出目录，则创建
    # 如果指定了输出目录路径
    if output_path:
        # 创建该目录，exist_ok=True 表示若目录已存在则不报错
        os.makedirs(output_path, exist_ok=True)

    # 定义内部函数：根据文件扩展名自动选择读取方式
    def read_table(file_path):
        # 获取文件扩展名（如 .xlsx）并转换为小写
        ext = os.path.splitext(file_path)[1].lower()
        # 如果扩展名是 .xlsx 或 .xls
        if ext in ['.xlsx', '.xls']:
            # 使用 pandas 读取 Excel 文件，第一列作为行索引
            return pd.read_excel(file_path, index_col=0)
        # 如果扩展名是 .csv
        elif ext == '.csv':
            # 读取 CSV 文件，第一列作为行索引
            return pd.read_csv(file_path, index_col=0)
        # 如果扩展名是 .txt
        elif ext == '.txt':
            # 先尝试用制表符 \t 作为分隔符读取
            df = pd.read_csv(file_path, index_col=0, sep='\t')
            # 如果读出来的列数 ≤1（说明可能分隔符不对）
            if df.shape[1] <= 1:
                # 改用逗号 , 作为分隔符重新读取
                df = pd.read_csv(file_path, index_col=0, sep=',')
            # 返回读取到的 DataFrame
            return df
        # 如果是不支持的文件格式
        else:
            # 抛出错误提示
            raise ValueError(f"不支持的文件格式: {ext}，请使用 .xlsx/.xls/.csv/.txt 文件")

    # 读取表达矩阵文件，返回原始数据（通常为基因×样本）
    counts_raw = read_table(counts_file)
    # 读取样本分组信息文件，返回元数据
    metadata = read_table(metadata_file)

    # 如果需要转置（pydeseq2 要求样本为行、基因为列）
    if transpose_counts:
        # 对表达矩阵进行转置，变为 样本×基因
        counts_df = counts_raw.T
    else:
        # 否则直接使用原始方向
        counts_df = counts_raw

    # 打印表达矩阵的形状（行数×列数），并提示为“样本×基因”
    print(f"表达矩阵形状: {counts_df.shape} (样本 x 基因)")
    # 打印元数据的形状
    print(f"元数据形状: {metadata.shape}")
    # 打印提示信息，开始检查样本名一致性
    print("样本名一致性检查：")
    # 计算表达矩阵的行索引（样本名）与元数据的行索引（样本名）的交集
    common = counts_df.index.intersection(metadata.index)
    # 打印表达矩阵中的样本数（行数）
    print(f"  表达矩阵样本数: {len(counts_df)}")
    # 打印元数据中的样本数
    print(f"  元数据样本数: {len(metadata)}")
    # 打印两者共有的样本数
    print(f"  共有样本数: {len(common)}")

    # 如果共有样本数少于任一方样本总数，说明存在不匹配的样本
    if len(common) < len(counts_df) or len(common) < len(metadata):
        # 给出警告信息
        print("警告：部分样本在表达矩阵或元数据中缺失，将自动对齐到共有样本。")
        # 将表达矩阵过滤为只保留共有样本
        counts_df = counts_df.loc[common]
        # 将元数据也过滤为只保留共有样本
        metadata = metadata.loc[common]

    # 返回处理后的表达矩阵和元数据
    return counts_df, metadata


# 当该脚本作为主程序运行时（而不是被导入为模块）
if __name__ == "__main__":
    # 定义数据文件夹的路径
    DATA_PATH_COUNTS = "D:/MyStudy/1-liulu_amino_acid/rna_seq"
    DATA_PATH_GROUPS = "D:/MyStudy/1-liulu_amino_acid/rna_seq"
    # 拼接出表达矩阵文件的完整路径（此处为 .xlsx 文件）
    counts_path = os.path.join(DATA_PATH_COUNTS, "counts.xlsx")  # 也可换成 .csv/.txt
    # 拼接出样本分组文件的完整路径
    meta_path = os.path.join(DATA_PATH_GROUPS, "group.xlsx")
    # 定义输出目录的相对路径
    OUTPUT_DIR = "../output_files/synthetic_example/"

    # 调用自定义函数加载数据，并接收返回的两个 DataFrame
    counts, meta = load_pydeseq2_data(
        counts_file=counts_path,
        metadata_file=meta_path,
        transpose_counts=True,  # 需要转置
        output_path=OUTPUT_DIR  # 创建输出目录
    )

    # 提示接下来可以用这些数据运行 pydeseq2
    print("\n表达矩阵前5行：")
    # 打印表达矩阵的前5行（预览）
    print(counts.head())
    # 打印元数据
    print("\n元数据：")
    print(meta)


