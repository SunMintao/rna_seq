import os
import pandas as pd


def split_by_group_and_export(input_excel, output_base_dir="."):
    """
    读取 group.xlsx，按 group 拆分数据，为每个 group 创建文件夹，
    并在文件夹中生成包含 samples 和 condition 两列的 Excel 文件。

    文件夹命名规则：根据该 group 内 condition='a' 和 'b' 的样本前缀自动生成，
    例如：L-WT_vs_L-lht1
    """
    # 读取 Excel 文件
    df = pd.read_excel(input_excel)

    # 按 group 分组
    for group_id, group_df in df.groupby('group'):
        # 提取该 group 内的样本和 condition
        sub_df = group_df[['samples', 'condition']].copy()

        # 获取 condition='a' 和 'b' 的样本前缀（去掉末尾的 -数字）
        def get_prefix(series):
            # 例如 'L-WT-1' -> 'L-WT'
            return series.str.replace(r'-\d+$', '', regex=True).iloc[0]

        prefix_a = get_prefix(group_df[group_df['condition'] == 'a']['samples'])
        prefix_b = get_prefix(group_df[group_df['condition'] == 'b']['samples'])

        # 生成文件夹名称
        folder_name = f"{prefix_a}_vs_{prefix_b}"
        folder_path = os.path.join(output_base_dir, folder_name)
        os.makedirs(folder_path, exist_ok=True)

        # 导出 Excel 文件
        output_excel = os.path.join(folder_path, f"{folder_name}.xlsx")
        sub_df.to_excel(output_excel, index=False)
        print(f"已生成: {output_excel}")


if __name__ == "__main__":
    # 请根据实际文件路径修改
    input_file = "group.xlsx"
    output_dir = "output_groups"  # 所有子文件夹将放在这个目录下
    split_by_group_and_export(input_file, output_dir)