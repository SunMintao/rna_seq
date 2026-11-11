import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns


class KEGGEnrichment:
    def __init__(self, pathway_file, deg_file, p_cutoff=0.05, output_prefix="kegg_"):
        """
        参数:
            pathway_file: 基因-通路映射文件路径 (Excel)
            deg_file: 差异基因列表文件路径 (Excel, 第一列为基因ID)
            p_cutoff: 显著性阈值 (校正后P值)
            output_prefix: 输出文件前缀
        """
        self.pathway_file = pathway_file
        self.deg_file = deg_file
        self.p_cutoff = p_cutoff
        self.output_prefix = output_prefix

        # 设置全局字体为 Arial
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['axes.unicode_minus'] = False

        # 存储解析后的数据
        self.pathway_to_genes = {}
        self.background_genes = set()
        self.deg_genes = []
        self.deg_in_background = []
        self.background_size = 0
        self.deg_size = 0

        # 富集结果
        self.results_df = None
        self.significant_df = None

    def load_data(self):
        """加载通路映射文件和差异基因列表"""
        df_pathway = pd.read_excel(self.pathway_file, header=0)
        for idx, row in df_pathway.iterrows():
            pathway = row['Pathway']
            gene_col = row['Gene_id']
            if pd.isna(gene_col):
                continue
            genes = [g.strip() for g in str(gene_col).split(';') if g.strip()]
            self.pathway_to_genes[pathway] = genes
            self.background_genes.update(genes)

        df_deg = pd.read_excel(self.deg_file, header=0)
        self.deg_genes = df_deg.iloc[:, 0].dropna().unique().tolist()
        self.background_size = len(self.background_genes)
        self.deg_in_background = [g for g in self.deg_genes if g in self.background_genes]
        self.deg_size = len(self.deg_in_background)

        print(f"背景基因总数: {self.background_size}")
        print(f"差异基因总数: {len(self.deg_genes)}")
        print(f"存在于背景集的差异基因数: {self.deg_size}")

    def run_enrichment(self):
        """超几何检验和多重校正，同时记录每个通路中差异基因的列表"""
        if not self.pathway_to_genes:
            raise ValueError("请先调用 load_data() 加载数据")

        results = []
        for pathway, genes in self.pathway_to_genes.items():
            pathway_set = set(genes)
            pathway_size = len(pathway_set)
            # 交集基因列表
            overlap_genes_list = [g for g in self.deg_in_background if g in pathway_set]
            overlap = len(overlap_genes_list)
            if overlap == 0:
                continue
            p_val = stats.hypergeom.sf(overlap - 1, self.background_size, pathway_size, self.deg_size)
            enrichment = (overlap / pathway_size) / (self.deg_size / self.background_size)
            results.append({
                'Pathway': pathway,
                'Pathway_size': pathway_size,
                'Overlap_genes': overlap,
                'Overlap_genes_list': overlap_genes_list,  # 存储基因列表
                'Enrichment_factor': enrichment,
                'P_value': p_val
            })

        self.results_df = pd.DataFrame(results)
        if self.results_df.empty:
            print("警告：未发现任何富集通路")
            return

        _, p_adj, _, _ = multipletests(self.results_df['P_value'], method='fdr_bh')
        self.results_df['Adjusted_P'] = p_adj
        self.significant_df = self.results_df[self.results_df['Adjusted_P'] < self.p_cutoff].copy()
        self.significant_df = self.significant_df.sort_values('Adjusted_P')
        print(f"显著富集通路数: {len(self.significant_df)}")

    def save_results(self, all_suffix="enrichment_all.csv", sig_suffix="enrichment_significant.csv"):
        """保存富集结果表格（不包含基因列表列，因为列表列不适合直接存为常规CSV）"""
        if self.results_df is not None:
            # 保存时去除列表列，避免CSV格式问题
            all_df_to_save = self.results_df.drop(
                columns=['Overlap_genes_list']) if 'Overlap_genes_list' in self.results_df.columns else self.results_df
            all_path = self.output_prefix + all_suffix
            all_df_to_save.to_csv(all_path, index=False)
            print(f"全部结果已保存至 {all_path}")
        if self.significant_df is not None and not self.significant_df.empty:
            sig_df_to_save = self.significant_df.drop(columns=[
                'Overlap_genes_list']) if 'Overlap_genes_list' in self.significant_df.columns else self.significant_df
            sig_path = self.output_prefix + sig_suffix
            sig_df_to_save.to_csv(sig_path, index=False)
            print(f"显著结果已保存至 {sig_path}")

    def save_significant_genes(self, suffix="significant_genes.csv"):
        """保存显著通路对应的差异基因列表文件"""
        if self.significant_df is None or self.significant_df.empty:
            print("没有显著通路，无法生成基因列表文件")
            return

        # 复制显著通路数据，将基因列表转换为分号分隔的字符串
        genes_df = self.significant_df[['Pathway', 'Overlap_genes', 'Overlap_genes_list']].copy()
        genes_df['Gene_list'] = genes_df['Overlap_genes_list'].apply(lambda x: ';'.join(x))
        genes_df = genes_df.drop(columns=['Overlap_genes_list'])
        # 重新排列列顺序
        genes_df = genes_df[['Pathway', 'Overlap_genes', 'Gene_list']]

        out_path = self.output_prefix + suffix
        genes_df.to_csv(out_path, index=False)
        print(f"显著通路差异基因列表已保存至 {out_path}")

    def plot_bubble(self, top_n=20, suffix="bubble.png"):
        """气泡图，字体 Arial"""
        if self.results_df is None or self.results_df.empty:
            print("无数据，无法绘制气泡图")
            return

        plt.rcParams['font.family'] = 'Arial'
        plot_df = self.results_df.sort_values('Adjusted_P').head(top_n).copy()
        plot_df['neg_log10_padj'] = -np.log10(plot_df['Adjusted_P'])

        plt.figure(figsize=(10, 8))
        scatter = plt.scatter(
            x=plot_df['Enrichment_factor'],
            y=plot_df['Pathway'],
            s=plot_df['Overlap_genes'] * 10,
            c=plot_df['neg_log10_padj'],
            cmap='Reds',
            alpha=0.8,
            edgecolors='black'
        )
        plt.colorbar(scatter, label='-log10(Adjusted P)')
        plt.xlabel('Enrichment Factor', fontfamily='Arial')
        plt.ylabel('KEGG Pathway', fontfamily='Arial')
        plt.title('KEGG Enrichment Bubble Plot', fontfamily='Arial')
        plt.tight_layout()
        out_path = self.output_prefix + suffix
        plt.savefig(out_path, dpi=300)
        plt.show()
        print(f"气泡图已保存至 {out_path}")

    def plot_bar(self, top_n=20, pvalue_type='raw', suffix="bar.png"):
        """柱形图，字体 Arial"""
        if self.results_df is None or self.results_df.empty:
            print("无数据，无法绘制柱形图")
            return

        plt.rcParams['font.family'] = 'Arial'

        if pvalue_type == 'raw':
            sort_col = 'P_value'
            xlabel = '-log10(P-value)'
            title = 'KEGG Enrichment Bar Plot (Raw P-value)'
        else:
            sort_col = 'Adjusted_P'
            xlabel = '-log10(Adjusted P-value)'
            title = 'KEGG Enrichment Bar Plot (Adjusted P-value)'

        plot_df = self.results_df.sort_values(sort_col).head(top_n).copy()
        plot_df['neg_log10_p'] = -np.log10(plot_df[sort_col])

        plt.figure(figsize=(10, 8))
        plt.barh(plot_df['Pathway'], plot_df['neg_log10_p'], color='steelblue')
        plt.xlabel(xlabel, fontfamily='Arial')
        plt.ylabel('KEGG Pathway', fontfamily='Arial')
        plt.title(title, fontfamily='Arial')
        plt.tight_layout()
        out_path = self.output_prefix + suffix
        plt.savefig(out_path, dpi=300)
        plt.show()
        print(f"柱形图已保存至 {out_path}")

    def run_all(self, save_results=True, save_genes=True, plot_bubble=True, plot_bar=True, bar_pvalue='adjusted',
                all_suffix="enrichment_all.csv", sig_suffix="enrichment_significant.csv",
                gene_suffix="significant_genes.csv", bubble_suffix="bubble.png", bar_suffix="bar.png"):
        """完整流程"""
        self.load_data()
        self.run_enrichment()
        if save_results:
            self.save_results(all_suffix=all_suffix, sig_suffix=sig_suffix)
        if save_genes:
            self.save_significant_genes(suffix=gene_suffix)
        if plot_bubble:
            self.plot_bubble(suffix=bubble_suffix)
        if plot_bar:
            self.plot_bar(pvalue_type=bar_pvalue, suffix=bar_suffix)


# 使用示例
if __name__ == "__main__":
    kegg = KEGGEnrichment(
        pathway_file="gene.pathway.xls",
        deg_file="H-WT_vs_H-OX_diff_genes.xlsx",
        output_prefix="my_kegg_"
    )
    kegg.run_all()