import os
import pandas as pd

def extract_group_tpm_by_names(group_names, tpm_data, output_dir='./output'):
    """
    根据指定的组名列表，从TPM表达矩阵中提取对应样本组的数据并保存为Excel文件。

    参数
    ----------
    group_names : list of str
        组名列表，每个元素格式应为 'A_vs_B'，表示比较组。
        例如 ['H-WT_vs_H-lht1', 'L-WT_vs_L-OX']
    tpm_data : str or pandas.DataFrame
        可以是TPM矩阵的Excel文件路径，或者是已读取的DataFrame。
        要求必须包含 'ID' 列，其他列名为样本标识。
    output_dir : str
        输出目录，默认为 './output'，函数会自动创建该目录。

    返回
    -------
    None
        函数无返回值，仅生成文件并打印处理信息。
    """
    # 1. 读取TPM数据（如果传入的是路径）
    if isinstance(tpm_data, str):
        if not os.path.exists(tpm_data):
            raise FileNotFoundError(f"TPM文件不存在: {tpm_data}")
        df = pd.read_excel(tpm_data)
    else:
        df = tpm_data.copy()  # 避免修改原始DataFrame

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    my_columns = df.columns.tolist()

    # 2. 解析每个组名，获得要提取的列名模式
    # 每个组名格式如 "A_vs_B"
    group_pairs = [name.split('_vs_') for name in group_names]

    # 3. 为每个比较组构建布尔索引
    indices = []
    for pair in group_pairs:
        bool_idx = []
        for col in my_columns:
            if col == 'ID':
                bool_idx.append(True)      # 保留基因ID列
            elif pair[0] in col or pair[1] in col:
                bool_idx.append(True)      # 列名包含组A或组B
            else:
                bool_idx.append(False)
        indices.append(bool_idx)

    # 4. 提取并保存
    for idx, name in zip(indices, group_names):
        group_df = df.loc[:, idx]
        out_path = os.path.join(output_dir, f"{name}_tpm.xlsx")
        group_df.to_excel(out_path, index=False)
        print(f"已保存: {out_path}")

if __name__ == "__main__":
    # 示例1：使用默认路径和组名列表
    my_groups = ['H-WT_vs_H-lht1', 'H-WT_vs_H-OX', 'L-WT_vs_L-lht1', 'L-WT_vs_L-OX']
    extract_group_tpm_by_names(
        group_names=my_groups,
        tpm_data='../2-get_gene_tpm/gene_tpm_value.xlsx',
        output_dir='../2-get_gene_tpm/'
    )

    # 示例2：先读取DataFrame再传入
    df_tpm = pd.read_excel('../2-get_gene_tpm/gene_tpm_value.xlsx')
    extract_group_tpm_by_names(my_groups, df_tpm, output_dir='./results')