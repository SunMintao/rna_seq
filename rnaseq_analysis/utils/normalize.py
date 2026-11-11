"""
normalize.py

Normalization methods for RNA-seq counts data.
"""


import pandas as pd
import numpy as np
import gzip
from collections import defaultdict

def cpm(counts, normalize_by_library=True, log=False, prior_count=1):
    """
    Counts Per Million (CPM) normalization.

    Parameters
    ----------
    counts : pd.DataFrame
        Counts matrix (genes as rows, samples as columns).
    normalize_by_library : bool, default True
        If True, divide by total library size (sum of counts per sample).
        If False, uses provided scaling factors (not implemented here).
    log : bool, default False
        If True, return log2(CPM + prior_count).
    prior_count : int or float, default 1
        Pseudocount added before log transformation (only if log=True).

    Returns
    -------
    pd.DataFrame
        CPM-normalized expression matrix (or log-CPM if log=True).
    """
    # Calculate library sizes (total counts per sample)
    lib_sizes = counts.sum(axis=0)
    # CPM = (counts / lib_size) * 1e6
    cpm_values = counts.div(lib_sizes, axis=1) * 1e6

    if log:
        cpm_values = np.log2(cpm_values + prior_count)

    return cpm_values


import pandas as pd
import numpy as np
import gzip
from collections import defaultdict


import pandas as pd
import numpy as np
import gzip
import warnings
from collections import defaultdict

def tpm_from_gtf(counts, gtf_path, normalize_by_length=True, log=False, prior_count=1):
    """
    从 GTF 文件计算基因长度，并执行 TPM 标准化。

    注：若基因无 exon 记录，则依次尝试使用 transcript 或 gene 行区间作为长度。
    对于以 'TRNA' 开头的基因，如果找不到长度信息，会跳过并给出警告，不中断程序。
    对于其他基因，若找不到长度信息，则抛出 ValueError。

    Parameters
    ----------
    counts : pd.DataFrame
        计数矩阵，行为基因（索引为 gene_id），列为样本。
    gtf_path : str
        GTF 文件路径（支持 .gtf 或 .gz）。
    normalize_by_length : bool, default True
        是否按基因长度校正（标准 TPM 为 True）。
    log : bool, default False
        是否返回 log2(TPM + prior_count)。
    prior_count : int or float, default 1
        log 变换前的伪计数。

    Returns
    -------
    pd.DataFrame
        TPM 标准化后的表达矩阵（基因 × 样本）。
    """

    # ---------- 辅助函数 ----------
    def _parse_attributes(attr_str):
        attrs = {}
        for field in attr_str.strip(';').split(';'):
            field = field.strip()
            if not field:
                continue
            parts = field.split(' ', 1)
            if len(parts) == 2:
                key = parts[0]
                val = parts[1].strip('"')
                attrs[key] = val
        return attrs

    def _merge_intervals(intervals):
        """合并重叠/相邻区间，返回总长度"""
        if not intervals:
            return 0
        intervals.sort()
        merged = []
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_end + 1:  # 允许相邻（+1）合并
                cur_end = max(cur_end, e)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = s, e
        merged.append((cur_start, cur_end))
        return sum(end - start + 1 for start, end in merged)

    # ---------- 1. 从 GTF 收集各类区间 ----------
    exon_intervals = defaultdict(list)
    transcript_intervals = defaultdict(list)
    gene_intervals = defaultdict(list)

    open_func = gzip.open if gtf_path.endswith('.gz') else open
    with open_func(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attrs = _parse_attributes(parts[8])
            gene_id = attrs.get('gene_id')
            if not gene_id:
                continue

            if feat_type == 'exon':
                exon_intervals[gene_id].append((start, end))
            elif feat_type == 'transcript':
                transcript_intervals[gene_id].append((start, end))
            elif feat_type == 'gene':
                gene_intervals[gene_id].append((start, end))

    # ---------- 2. 按优先级计算基因长度 ----------
    gene_lengths_dict = {}
    all_gene_ids = set(exon_intervals) | set(transcript_intervals) | set(gene_intervals)
    for gene_id in all_gene_ids:
        intervals = exon_intervals.get(gene_id)
        if not intervals:
            intervals = transcript_intervals.get(gene_id)
        if not intervals:
            intervals = gene_intervals.get(gene_id)
        if intervals:
            gene_lengths_dict[gene_id] = _merge_intervals(intervals)

    gene_lengths = pd.Series(gene_lengths_dict)
    counts.index = counts.index.astype(str)

    # ---------- 3. 检查缺失（核心修改） ----------
    missing = set(counts.index) - set(gene_lengths.index)

    if missing:
        # 分离以 'TRNA' 开头的基因（用户指定不报错）
        trna_missing = {g for g in missing if g.startswith('TRNA')}
        other_missing = missing - trna_missing

        # 非 TRNA 基因缺失 -> 抛出错误
        if other_missing:
            raise ValueError(
                f"以下非 TRNA 基因在 GTF 文件中找不到任何可用于计算长度的特征（exon/transcript/gene）: {other_missing}"
            )

        # TRNA 基因缺失 -> 警告并跳过
        if trna_missing:
            warnings.warn(
                f"以下 {len(trna_missing)} 个 TRNA 基因在 GTF 中找不到长度信息，将被跳过（不参与 TPM 计算）: {trna_missing}"
            )
            # 只保留非 TRNA 缺失的基因（即移除 trna_missing）
            keep_genes = set(counts.index) - trna_missing
            counts = counts.loc[list(keep_genes)]
            gene_lengths = gene_lengths.reindex(keep_genes)

            # 如果过滤后没有基因，直接报错
            if counts.empty:
                raise ValueError("过滤后没有剩余基因，无法计算 TPM。")
        else:
            # 此分支理论上不会执行，因为 missing 非空且没有 trna_missing 时，other_missing 必然非空，已提前报错
            pass
    else:
        # 无缺失，正常对齐
        gene_lengths = gene_lengths.reindex(counts.index)

    # ---------- 4. 计算 TPM ----------
    if normalize_by_length:
        rpk = counts.div(gene_lengths, axis=0)
    else:
        rpk = counts

    scaling_factors = rpk.sum(axis=0) / 1e6
    tpm_values = rpk.div(scaling_factors, axis=1)

    if log:
        tpm_values = np.log2(tpm_values + prior_count)

    return tpm_values


def deseq2_size_factors(counts, pseudo_count=1, geo_mean_trim=0.2, return_df=False):
    """
    Calculate DESeq2 median-of-ratios size factors.

    This method estimates sample-specific scaling factors assuming most genes
    are not differentially expressed.

    Parameters
    ----------
    counts : pd.DataFrame
        Counts matrix (genes as rows, samples as columns).
    pseudo_count : int or float, default 1
        Pseudocount added to avoid log(0).
    geo_mean_trim : float, default 0.2
        Fraction of genes to trim from both ends when calculating geometric mean.
        A typical DESeq2 uses all genes, but trimming can remove extreme values.
    return_df : bool, default False
        If True, return normalized counts (counts divided by size factors).
        If False, return only the size factors as a Series.

    Returns
    -------
    pd.Series or pd.DataFrame
        If return_df=False: size factors (index = samples).
        If return_df=True: normalized counts matrix.
    """
    # Add pseudocount to avoid zeros
    counts_adj = counts + pseudo_count

    # Calculate geometric mean of each gene across all samples
    # (using log to avoid underflow)
    log_counts = np.log(counts_adj)
    # Geometric mean = exp(mean(log(x)))
    # Apply trimming if specified
    if geo_mean_trim > 0:
        # For each gene, compute mean of log counts, but trim extreme samples?
        # DESeq2 uses all samples for the ratio, but trimming can be applied across genes.
        # Here we compute per-gene geometric mean across samples, then trim genes.
        geo_means = np.exp(log_counts.mean(axis=1))
        # Trim genes with extreme geometric means (optional)
        if geo_mean_trim < 1:
            lower = np.percentile(geo_means, 100 * geo_mean_trim)
            upper = np.percentile(geo_means, 100 * (1 - geo_mean_trim))
            keep = (geo_means >= lower) & (geo_means <= upper)
            counts_adj = counts_adj.loc[keep]
            log_counts = log_counts.loc[keep]
        geo_means = np.exp(log_counts.mean(axis=1))
    else:
        geo_means = np.exp(log_counts.mean(axis=1))

    # For each sample, compute median of (count / geometric_mean) over genes
    size_factors = {}
    for sample in counts_adj.columns:
        ratios = counts_adj[sample] / geo_means
        # Median of ratios, ignoring NaNs/infs
        sf = np.median(ratios[~np.isnan(ratios) & np.isfinite(ratios)])
        size_factors[sample] = sf

    size_factors = pd.Series(size_factors)

    if return_df:
        norm_counts = counts.div(size_factors, axis=1)
        return norm_counts
    else:
        return size_factors


def normalize_by_size_factors(counts, size_factors):
    """
    Apply pre-calculated size factors to normalize counts.

    Parameters
    ----------
    counts : pd.DataFrame
        Counts matrix.
    size_factors : pd.Series
        Size factors (index = sample names).

    Returns
    -------
    pd.DataFrame
        Normalized counts (counts divided by size factors).
    """
    return counts.div(size_factors, axis=1)


def rle(counts, geo_mean_trim=0):
    """
    Relative Log Expression (RLE) normalization (used by edgeR/DESeq2 alternative).

    Parameters
    ----------
    counts : pd.DataFrame
        Counts matrix.
    geo_mean_trim : float, default 0
        Fraction of genes to trim when computing geometric mean (0 = no trim).

    Returns
    -------
    pd.Series
        Size factors.
    """
    # Pseudocount to avoid zeros
    counts_adj = counts + 1
    # Geometric mean per gene (across samples)
    log_counts = np.log(counts_adj)
    geo_means = np.exp(log_counts.mean(axis=1))
    # Trim if specified
    if geo_mean_trim > 0:
        lower = np.percentile(geo_means, 100 * geo_mean_trim)
        upper = np.percentile(geo_means, 100 * (1 - geo_mean_trim))
        keep = (geo_means >= lower) & (geo_means <= upper)
        counts_adj = counts_adj.loc[keep]
        geo_means = geo_means[keep]
    # For each sample, compute median of ratio to geometric mean
    size_factors = {}
    for sample in counts_adj.columns:
        ratio = counts_adj[sample] / geo_means
        size_factors[sample] = np.median(ratio)
    return pd.Series(size_factors)