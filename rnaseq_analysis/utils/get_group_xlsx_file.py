import pandas as pd


def generate_group_file_from_pairs(pair_list, replicates=3, save_path=None):
    """
    根据处理对列表自动生成分组文件（三列：samples, condition, group）。

    参数:
        pair_list : list of list/tuple
            每个元素为一个组内的处理名称列表，长度不限（通常为2）。
            例如: [['L-WT', 'L-lht1'], ['L-WT', 'L-OX'], ...]
            函数会自动为每个组内的处理分配条件：第1个 → 'a'，第2个 → 'b'，第3个 → 'c'，依此类推。
            组号从1开始自动递增。
        replicates : int, default=3
            每个处理的重复次数。
        save_path : str, optional
            若提供，将生成的 DataFrame 保存为 CSV 文件。

    返回:
        pandas.DataFrame
    """
    rows = []
    condition_labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']  # 足够支持多数情况

    for group_id, treatments in enumerate(pair_list, start=1):
        # 为每个处理分配条件
        for idx, treat_name in enumerate(treatments):
            cond = condition_labels[idx] if idx < len(condition_labels) else f'cond{idx + 1}'
            for rep in range(1, replicates + 1):
                sample = f"{treat_name}-{rep}"
                rows.append([sample, cond, group_id])

    df = pd.DataFrame(rows, columns=['samples', 'condition', 'group'])

    if save_path:
        df.to_excel(save_path, index=False)
        print(f"文件已保存至: {save_path}")

    return df


if __name__ == '__main__':
    pairs = [
        ['L-WT', 'L-lht1'],
        ['L-WT', 'L-OX'],
        ['H-WT', 'H-lht1'],
        ['H-WT', 'H-OX'],
        ['L-OX', 'L-lht1'],
        ['H-OX', 'H-lht1']
    ]

    df = generate_group_file_from_pairs(pairs)
    print(df.head(12))