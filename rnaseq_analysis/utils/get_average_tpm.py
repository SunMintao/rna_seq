#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
from pathlib import Path


def infer_group_names(columns):
    """
    根据列名自动推断组名。
    规则：去掉末尾的 '-数字' （如 -1, -2, -3...），剩余部分作为组名。
    例如：'L-WT-1' -> 'L-WT', 'H-lht1-2' -> 'H-lht1'
    """
    groups = {}
    for col in columns:
        col_stripped = str(col).strip()  # 去除可能的首尾空格
        # 匹配末尾的 '-数字'（整数）
        match = re.search(r'-\d+$', col_stripped)
        if match:
            group = col_stripped[:match.start()]
        else:
            group = col_stripped
        groups[col] = group
    return groups


def calculate_group_means(df, group_map):
    """
    对每个基因，计算每个组内所有样本的均值（行均值）。
    df : DataFrame，第一列为基因ID，其余列为样本。
    group_map : dict，列名 -> 组名
    """
    # 将样本列按组名合并
    grouped_columns = {}
    for col, group in group_map.items():
        grouped_columns.setdefault(group, []).append(col)

    # 对每个组计算均值
    result = pd.DataFrame(index=df.index)  # 临时以原始索引为索引
    result['ID'] = df.iloc[:, 0]  # 基因ID列

    for group, cols in grouped_columns.items():
        # 确保选取的列存在
        existing_cols = [c for c in cols if c in df.columns]
        if not existing_cols:
            print(f"警告：组 {group} 没有找到有效列，跳过")
            result[group] = pd.NA
            continue
        group_data = df[existing_cols]
        # 强制转换为数值类型，无法转换的变为 NaN
        group_data = group_data.apply(pd.to_numeric, errors='coerce')
        # 计算行均值，忽略 NaN
        result[group] = group_data.mean(axis=1, skipna=True)

    # 调整列顺序：ID 放在第一列
    cols = ['ID'] + [c for c in result.columns if c != 'ID']
    result = result[cols]
    return result


def compute_group_means(input_file, output_file='group_means.xlsx', sheet_name=0):
    """
    计算TPM均值（按实验组自动分组）

    参数
    ----------
    input_file : str
        输入的Excel文件路径
    output_file : str, optional
        输出的Excel文件路径（默认：'group_means.xlsx'）
    sheet_name : str or int, optional
        工作表名称或索引（默认：0，即第一个工作表）
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"文件 {input_path} 不存在！")

    # 读取 Excel，保持原始数据类型，不自动转换
    df = pd.read_excel(input_path, sheet_name=sheet_name, dtype=str, index_col=None)

    # 打印原始列名，帮助调试
    print("原始列名：", list(df.columns))

    if df.shape[1] < 2:
        raise ValueError("数据至少需要包含基因ID列和一列样本数据！")

    # 第一列为基因ID列名（保留原始名称）
    gene_col_name = df.columns[0]
    # 其余为样本列
    sample_cols = list(df.columns[1:])

    # 将基因ID列和样本数据分离，并转换样本数据为数值
    gene_ids = df[gene_col_name]
    sample_df = df[sample_cols].apply(pd.to_numeric, errors='coerce')

    # 重新拼接成完整的 DataFrame（基因ID + 数值样本）
    df_numeric = pd.concat([gene_ids, sample_df], axis=1)

    # 自动推断组名（基于样本列名）
    group_map = infer_group_names(sample_cols)
    print("自动识别的组名及对应样本列：")
    for group, cols in pd.Series(group_map).groupby(group_map).apply(lambda x: list(x.index)).items():
        print(f"  {group}: {cols}")

    # 计算每组均值
    result_df = calculate_group_means(df_numeric, group_map)

    # 保存结果
    output_path = Path(output_file)
    result_df.to_excel(output_path, index=False)
    print(f"结果已保存至：{output_path}")

# ========== 使用示例 ==========
if __name__ == '__main__':
    # 请修改为您的实际文件路径
    input_excel = "gene_tpm_value.xlsx"
    output_excel = "group_means.xlsx"

    compute_group_means(input_excel, output_excel)