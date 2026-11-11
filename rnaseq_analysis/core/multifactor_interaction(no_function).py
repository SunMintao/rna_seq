"""
完整 DESeq2 交互分析 + 全局平均主效应提取 (Python + inmoose)
实验设计：叶片（或果实），3个开花时期 × 4个低温时间点
参考水平：stage = pre2d, time = 0h
"""

import pandas as pd
import numpy as np
import re
from inmoose.deseq2 import DESeqDataSet, DESeq
from inmoose.deseq2.utils import results


def parse_sample_name(sample_name: str):
    """
    解析样本名，返回 stage 和 time。
    支持格式：L-2D-0-1 (pre2d, 0h) 或 L0D-12-1 (anthesis, 12h)
    """
    if 'L-2D' in sample_name:
        stage = 'pre2d'
    elif 'L0D' in sample_name:
        stage = 'anthesis'
    elif 'L2D' in sample_name:
        stage = 'post2d'
    else:
        raise ValueError(f"Unknown stage in {sample_name}")

    # 提取时间数字
    match = re.search(r'-(\d+)-[^-]+$', sample_name)
    if match is None:
        match = re.search(r'-(\d+)-', sample_name)
    if match is None:
        raise ValueError(f"Cannot parse time from {sample_name}")
    time_num = match.group(1)
    time = f"{time_num}h"
    return stage, time


def create_metadata(sample_names: list):
    """从样本名列表生成 metadata DataFrame"""
    records = []
    for s in sample_names:
        stage, time = parse_sample_name(s)
        records.append({'sample': s, 'stage': stage, 'time': time})
    metadata = pd.DataFrame(records)
    metadata.set_index('sample', inplace=True)
    # 设置分类顺序（参考水平为第一个）
    metadata['stage'] = pd.Categorical(metadata['stage'],
                                       categories=['pre2d', 'anthesis', 'post2d'],
                                       ordered=False)
    metadata['time'] = pd.Categorical(metadata['time'],
                                      categories=['0h', '6h', '12h', '24h'],
                                      ordered=False)
    return metadata


def load_and_prepare_data(file_path: str, min_total_count: int = 10):
    """
    读取 Excel 文件，过滤低表达基因，返回 counts DataFrame 和 metadata
    """
    counts_df = pd.read_excel(file_path, index_col=0)  # 第一列为 gene_id
    # 确保列名为样本名
    sample_names = counts_df.columns.tolist()
    metadata = create_metadata(sample_names)
    # 保证列顺序与 metadata 一致
    counts_df = counts_df[metadata.index]
    # 过滤低表达 (总 counts >= min_total_count)
    keep = counts_df.sum(axis=1) >= min_total_count
    counts_df = counts_df[keep]
    print(f"过滤后保留 {counts_df.shape[0]} 个基因")
    return counts_df, metadata


def run_deseq2_analysis(counts_df, metadata, design_formula="~ stage * time"):
    """
    运行 DESeq2 分析（完整模型），返回 dds 对象
    """
    # 转置为样本行、基因列
    counts_t = counts_df.T
    # 确保索引一致
    assert all(counts_t.index == metadata.index)
    dds = DESeqDataSet(counts=counts_t, colData=metadata, design=design_formula)
    dds = DESeq(dds, test="Wald", betaPrior=False, quiet=False)
    return dds


def run_lrt(dds, reduced_formula="~ stage + time"):
    """
    运行 LRT 全局交互检验
    """
    dds_lrt = DESeq(dds, test="LRT", reduced=reduced_formula, quiet=False)
    res_lrt = results(dds_lrt)
    return res_lrt, dds_lrt


def extract_interaction(dds, stage_level, time_level, ref_stage='pre2d', ref_time='0h'):
    """
    提取特定交互项的结果 (Wald)
    例如: stage_level='post2d', time_level='24h' -> post2d × 24h
    """
    coef_name = f"stage_{stage_level}_vs_{ref_stage}:time_{time_level}_vs_{ref_time}"
    # 检查是否存在于 resultsNames
    available = dds.resultsNames()
    if coef_name not in available:
        # 尝试替换冒号为点号（inmoose可能用点）
        alt_name = f"stage_{stage_level}_vs_{ref_stage}.time_{time_level}_vs_{ref_time}"
        if alt_name in available:
            coef_name = alt_name
        else:
            raise ValueError(f"Interaction {coef_name} not found. Available: {available}")
    res = results(dds, name=coef_name)
    return res


def extract_all_main_effects(dds_obj):
    """
    从 DESeq2 对象中提取所有主效应系数（排除 Intercept 和包含 ':' 的交互项）
    返回 DataFrame: gene, 各系数_FC, _pvalue, _padj
    """
    all_coef = dds_obj.resultsNames()
    # 主效应：不包含 ':' 且不是 'Intercept'
    main_coef = [c for c in all_coef if ':' not in c and c != 'Intercept']
    print("提取的主效应:", main_coef)

    result_dfs = []
    for coef in main_coef:
        res = results(dds_obj, name=coef)
        df = pd.DataFrame({
            'gene': res.index,
            f'{coef}_FC': res.log2FoldChange,
            f'{coef}_pvalue': res.pvalue,
            f'{coef}_padj': res.padj
        })
        result_dfs.append(df)

    # 合并所有
    if result_dfs:
        merged = result_dfs[0]
        for df in result_dfs[1:]:
            merged = pd.merge(merged, df, on='gene', how='outer')
        return merged
    else:
        return pd.DataFrame()


def results_to_df(res):
    """将 DESeqResults 转为 DataFrame"""
    df = pd.DataFrame({
        'gene': res.index,
        'log2FoldChange': res.log2FoldChange,
        'pvalue': res.pvalue,
        'padj': res.padj,
        'baseMean': res.baseMean,
        'lfcSE': res.lfcSE,
        'stat': res.stat
    })
    return df


def classify_genes(final_df, lrt_padj_col='lrt_padj', time_avg_padj_col='time_24h_vs_0h_padj'):
    """
    根据 LRT 显著性和平均低温主效应显著性分类
    """
    def classify(row):
        lrt_sig = row.get(lrt_padj_col, np.nan) < 0.05 if not pd.isna(row.get(lrt_padj_col, np.nan)) else False
        time_sig = row.get(time_avg_padj_col, np.nan) < 0.05 if not pd.isna(row.get(time_avg_padj_col, np.nan)) else False
        if lrt_sig and time_sig:
            return "交互显著 + 低温主效应显著"
        elif lrt_sig and not time_sig:
            return "交互显著 + 低温主效应不显著（花期特异开关）"
        elif not lrt_sig and time_sig:
            return "交互不显著 + 仅低温主效应显著（通用冷响应）"
        else:
            return "交互不显著 + 低温主效应不显著（无响应）"

    final_df['Category'] = final_df.apply(classify, axis=1)
    return final_df


# ================== 主程序 ==================
if __name__ == "__main__":
    # 1. 加载数据
    file_path = "../../data/leaf_matrix_gene.count.xlsx"  # 请修改实际路径
    counts_df, metadata = load_and_prepare_data(file_path, min_total_count=10)

    # 2. 运行完整模型
    dds = run_deseq2_analysis(counts_df, metadata, "~ stage * time")

    # 3. 查看系数名
    print("完整模型系数名:", dds.resultsNames())

    # 4. LRT 全局交互
    res_lrt, dds_reduced = run_lrt(dds, "~ stage + time")
    df_lrt = results_to_df(res_lrt).rename(columns={
        'log2FoldChange': 'lrt_log2FC',
        'pvalue': 'lrt_pvalue',
        'padj': 'lrt_padj',
        'stat': 'lrt_stat'
    })

    # 5. 提取特定交互项 (post2d × 24h)
    res_int = extract_interaction(dds, 'post2d', '24h')
    df_int = results_to_df(res_int).rename(columns={
        'log2FoldChange': 'int_log2FC',
        'pvalue': 'int_pvalue',
        'padj': 'int_padj'
    })[['gene', 'int_log2FC', 'int_pvalue', 'int_padj']]

    # 6. 提取缩减模型的所有全局平均主效应
    df_main_avg = extract_all_main_effects(dds_reduced)

    # 7. 合并所有结果
    final_df = df_lrt[['gene', 'lrt_padj', 'lrt_stat']].merge(
        df_int, on='gene', how='outer'
    ).merge(
        df_main_avg, on='gene', how='outer'
    )

    # 8. 添加分类标签 (使用 24h 平均主效应)
    final_df = classify_genes(final_df,
                              lrt_padj_col='lrt_padj',
                              time_avg_padj_col='time_24h_vs_0h_padj')

    # 9. 保存结果
    final_df.to_csv("DESeq2_complete_results_python.csv", index=False)
    # 保存各分类子集
    for cat in final_df['Category'].unique():
        subset = final_df[final_df['Category'] == cat]
        subset.to_csv(f"Category_{cat.replace(' ', '_')}.csv", index=False)

    # 保存主效应表
    df_main_avg.to_csv("reduced_model_main_effects_all_python.csv", index=False)

    print("分析完成！结果文件已保存。")
    print("各类别基因数量:")
    print(final_df['Category'].value_counts())

    # 可选：显示前几行
    print(final_df.head())