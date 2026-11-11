import os
import pickle as pkl
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def run_deseq2(counts_df, metadata, output_path, contrast_group=None,
               design="~condition", refit_cooks=True, n_cpus=8):
    """
    运行 pydeseq2 差异表达分析，并保存结果。

    Parameters
    ----------
    counts_df : pd.DataFrame
        表达矩阵，行名=样本，列名=基因，值为整数counts。
    metadata : pd.DataFrame
        样本元数据，行名=样本，至少包含 condition 列（或设计公式中的变量）。
    output_path : str
        结果保存目录（会自动创建）。
    contrast_group : list or None, default=None
        对比组列表，例如 ["condition", "Treatment", "Control"]。
        如果为 None，则不会执行统计检验（仅运行 deseq2 并保存 dds）。
    design : str, default="~condition"
        设计公式字符串，如 "~condition" 或 "~batch+condition"。
    refit_cooks : bool, default=True
        是否重新拟合异常值（Cook's distance）。
    n_cpus : int, default=8
        并行计算的 CPU 核心数。

    Returns
    -------
    dds : DeseqDataSet
        拟合后的 DESeq2 数据集对象。
    stats : DeseqStats or None
        如果提供了 contrast_group，则返回统计结果对象；否则返回 None。
    """
    # 创建输出目录
    os.makedirs(output_path, exist_ok=True)

    # 设置并行推断
    inference = DefaultInference(n_cpus=n_cpus)

    # 创建 DESeq2 数据集并拟合
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design=design,
        refit_cooks=refit_cooks,
        inference=inference,
    )
    dds.deseq2()

    # 保存 dds 对象（转换为 anndata 格式）
    with open(os.path.join(output_path, "result_adata.pkl"), "wb") as f:
        pkl.dump(dds.to_picklable_anndata(), f)

    # 如果提供了对比组，则进行统计检验
    stats = None
    if contrast_group is not None:
        stats = DeseqStats(dds, contrast=contrast_group, inference=inference)
        stats.summary()  # 计算统计结果
        # 保存结果表格
        results_path = os.path.join(output_path, "results.xlsx")
        stats.results_df.to_excel(results_path)
        print(f"差异表达结果已保存至: {results_path}")
    else:
        print("未提供对比组，仅保存了 dds 对象。")

    return dds, stats


# 使用示例（根据原代码）
if __name__ == "__main__":
    # 假设 counts_df, metadata, OUTPUT_PATH 已定义
    dds, stats = run_deseq2(
        counts_df=counts_df,
        metadata=metadata,
        output_path=OUTPUT_PATH,
        contrast_group=["condition", "a", "b"],  # 实验组a，对照组b
        design="~condition",
        refit_cooks=True,
        n_cpus=8
    )