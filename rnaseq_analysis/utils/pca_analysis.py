import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
import warnings
warnings.filterwarnings('ignore')

class PCAAnalyzer:
    """
    对基因表达数据进行 PCA 分析，生成可视化图形和输出文件。
    支持按样本名前缀自动分组，并在二维 PCA 图中绘制 95% 置信椭圆。
    输入 Excel 文件：行 = 基因，列 = 样本（第一列为基因 ID）。
    """

    def __init__(self, file_path, sheet_name=0, index_col=0):
        self.file_path = file_path
        self.sheet_name = sheet_name
        self.index_col = index_col
        self.data = None
        self.X = None
        self.X_raw = None
        self.sample_names = None
        self.gene_names = None
        self.pca = None
        self.scores = None
        self.loadings = None
        self.explained_variance_ratio = None
        self.components = None
        self.groups = None

    def load_data(self):
        df = pd.read_excel(self.file_path, sheet_name=self.sheet_name, index_col=self.index_col)
        self.data = df.copy()
        df_T = df.T
        self.sample_names = df_T.index.tolist()
        self.gene_names = df_T.columns.tolist()
        self.X_raw = df_T.values.astype(float)
        self._extract_groups()
        print(f"数据加载完成：{self.X_raw.shape[0]} 个样本，{self.X_raw.shape[1]} 个基因")
        print(f"发现分组：{set(self.groups)}")

    def _extract_groups(self):
        groups = []
        for name in self.sample_names:
            parts = name.split('-')
            if len(parts) > 1 and parts[-1].isdigit():
                group = '-'.join(parts[:-1])
            else:
                group = name
            groups.append(group)
        self.groups = groups

    def standardize(self):
        scaler = StandardScaler()
        self.X = scaler.fit_transform(self.X_raw)
        print("数据标准化完成")

    def perform_pca(self, n_components=None):
        if self.X is None:
            if self.X_raw is not None:
                print("数据未标准化，直接使用原始数据进行 PCA（注意：量纲差异可能影响结果）。")
                self.X = self.X_raw.copy()
            else:
                raise ValueError("请先调用 load_data() 加载数据")
        if n_components is None:
            n_components = min(self.X.shape)
        self.pca = PCA(n_components=n_components)
        self.scores = self.pca.fit_transform(self.X)
        self.components = self.pca.components_
        self.explained_variance_ratio = self.pca.explained_variance_ratio_
        self.loadings = pd.DataFrame(
            self.components.T,
            index=self.gene_names,
            columns=[f"PC{i+1}" for i in range(n_components)]
        )
        print(f"PCA 完成，保留 {n_components} 个主成分，累计方差解释率: {self.explained_variance_ratio.sum():.2%}")

    def _confidence_ellipse(self, x, y, ax, n_std=2.447, **kwargs):
        cov = np.cov(x, y)
        mean_x, mean_y = np.mean(x), np.mean(y)
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
        width, height = 2 * n_std * np.sqrt(vals)
        ellipse = Ellipse(xy=(mean_x, mean_y), width=width, height=height,
                          angle=theta, **kwargs)
        ax.add_patch(ellipse)

    def plot_2d(self, save_path=None, figsize=(10, 8), conf_level=0.95, group_order=None):
        """
        绘制二维 PCA 散点图（PC1 vs PC2），按分组着色并绘制置信椭圆

        参数:
            group_order (list): 自定义分组显示顺序（图例顺序），如 ['Control', 'Treatment']
        """
        if self.scores is None:
            raise RuntimeError("请先执行 perform_pca()")
        if self.scores.shape[1] < 2:
            print("主成分不足2个，无法绘制二维图")
            return

        plt.figure(figsize=figsize)
        ax = plt.gca()

        # 获取所有分组并处理排序
        unique_groups = set(self.groups)
        if group_order is None:
            ordered_groups = sorted(unique_groups)  # 默认按字母排序，更稳定
        else:
            # 检查用户提供的顺序是否完整
            ordered_groups = []
            for g in group_order:
                if g in unique_groups:
                    ordered_groups.append(g)
                else:
                    print(f"警告：分组 '{g}' 不在数据中，已忽略")
            # 补充未列出的分组（放在末尾）
            for g in sorted(unique_groups):
                if g not in ordered_groups:
                    ordered_groups.append(g)
                    print(f"注意：分组 '{g}' 未在 group_order 中，已追加到末尾")

        # 颜色生成（按顺序分配）
        n_groups = len(ordered_groups)
        if n_groups <= 10:
            colors = [plt.cm.tab10(i) for i in range(n_groups)]
        elif n_groups <= 20:
            colors = [plt.cm.tab20(i) for i in range(n_groups)]
        else:
            colors = [plt.cm.hsv(i / n_groups) for i in range(n_groups)]

        from scipy.stats import chi2
        chi2_val = chi2.ppf(conf_level, df=2)
        n_std = np.sqrt(chi2_val)

        for group, color in zip(ordered_groups, colors):
            idx = [i for i, g in enumerate(self.groups) if g == group]
            if not idx:
                continue
            pc1_vals = self.scores[idx, 0]
            pc2_vals = self.scores[idx, 1]
            plt.scatter(pc1_vals, pc2_vals, label=group, color=color, s=60, alpha=0.8)
            if len(idx) >= 3:
                self._confidence_ellipse(pc1_vals, pc2_vals, ax, n_std=n_std,
                                         edgecolor=color, facecolor='none', linewidth=1.5, linestyle='--')
            plt.scatter(np.mean(pc1_vals), np.mean(pc2_vals), marker='X', s=100, color=color, edgecolors='black')

        plt.xlabel(f'PC1 ({self.explained_variance_ratio[0]:.2%})')
        plt.ylabel(f'PC2 ({self.explained_variance_ratio[1]:.2%})')
        plt.title('PCA 2D Projection with 95% Confidence Ellipses')
        plt.legend(title='Group')
        plt.grid(alpha=0.3)

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"二维图已保存至 {save_path}")
        else:
            plt.show()

    def plot_3d(self, save_path=None, figsize=(10, 8), group_order=None):
        """
        绘制三维 PCA 散点图（前三个主成分），按分组着色

        参数:
            group_order (list): 自定义分组显示顺序（图例顺序）
        """
        if self.scores is None:
            raise RuntimeError("请先执行 perform_pca()")
        if self.scores.shape[1] < 3:
            print("警告：主成分数不足3，无法绘制三维图")
            return

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        unique_groups = set(self.groups)
        if group_order is None:
            ordered_groups = sorted(unique_groups)
        else:
            ordered_groups = []
            for g in group_order:
                if g in unique_groups:
                    ordered_groups.append(g)
                else:
                    print(f"警告：分组 '{g}' 不在数据中，已忽略")
            for g in sorted(unique_groups):
                if g not in ordered_groups:
                    ordered_groups.append(g)
                    print(f"注意：分组 '{g}' 未在 group_order 中，已追加到末尾")

        n_groups = len(ordered_groups)
        if n_groups <= 10:
            colors = [plt.cm.tab10(i) for i in range(n_groups)]
        elif n_groups <= 20:
            colors = [plt.cm.tab20(i) for i in range(n_groups)]
        else:
            colors = [plt.cm.hsv(i / n_groups) for i in range(n_groups)]

        for group, color in zip(ordered_groups, colors):
            idx = [i for i, g in enumerate(self.groups) if g == group]
            if not idx:
                continue
            xs = self.scores[idx, 0]
            ys = self.scores[idx, 1]
            zs = self.scores[idx, 2]
            ax.scatter(xs, ys, zs, label=group, color=color, s=50, alpha=0.8)
            ax.scatter(np.mean(xs), np.mean(ys), np.mean(zs), marker='X', s=100, color=color, edgecolors='black')

        ax.set_xlabel(f'PC1 ({self.explained_variance_ratio[0]:.2%})')
        ax.set_ylabel(f'PC2 ({self.explained_variance_ratio[1]:.2%})')
        ax.set_zlabel(f'PC3 ({self.explained_variance_ratio[2]:.2%})')
        ax.set_title('PCA 3D Projection')
        ax.legend(title='Group')

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"三维图已保存至 {save_path}")
        else:
            plt.show()

    def save_loadings(self, output_path="loadings.csv"):
        if self.loadings is None:
            raise RuntimeError("请先执行 perform_pca()")
        self.loadings.to_csv(output_path, index=True)
        print(f"载荷矩阵已保存至 {output_path}")

    def save_scores(self, output_path="factor_scores.csv"):
        if self.scores is None:
            raise RuntimeError("请先执行 perform_pca()")
        scores_df = pd.DataFrame(
            self.scores,
            index=self.sample_names,
            columns=[f"PC{i+1}" for i in range(self.scores.shape[1])]
        )
        scores_df.to_csv(output_path, index=True)
        print(f"样本主成分得分表已保存至 {output_path}")

    def save_factor_analysis_result(self, output_path="factor_analysis_result.txt"):
        if self.pca is None:
            raise RuntimeError("请先执行 perform_pca()")
        eigenvalues = self.pca.explained_variance_
        ev_ratio = self.explained_variance_ratio
        cum_ratio = np.cumsum(ev_ratio)
        n_comp = len(eigenvalues)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("========== PCA 因子分析结果 ==========\n\n")
            f.write("特征值 (Eigenvalues):\n")
            for i, val in enumerate(eigenvalues, 1):
                f.write(f"PC{i}: {val:.4f}\n")
            f.write("\n方差解释比例 (Explained Variance Ratio):\n")
            for i, val in enumerate(ev_ratio, 1):
                f.write(f"PC{i}: {val:.4f} ({val:.2%})\n")
            f.write("\n累计方差解释比例:\n")
            for i, val in enumerate(cum_ratio, 1):
                f.write(f"前{i}个主成分: {val:.4f} ({val:.2%})\n")
            f.write("\n载荷矩阵 (Loadings，前5个主成分，仅展示前10个基因):\n")
            f.write(self.loadings.iloc[:10, :min(5, n_comp)].to_string())
            f.write("\n\n完整载荷矩阵请参见 loadings.csv 文件")
        print(f"因子分析结果已保存至 {output_path}")

    def run_full_analysis(self, n_components=None, output_dir=".", standardize=False, group_order=None):
        """
        执行完整的分析流程，生成所有输出

        参数:
            n_components (int): 保留的主成分个数
            output_dir (str): 输出文件保存目录
            standardize (bool): 是否对数据进行标准化（默认 False）
            group_order (list): 自定义分组显示顺序（图例顺序）
        """
        import os
        self.load_data()
        if standardize:
            self.standardize()
        self.perform_pca(n_components=n_components)
        self.plot_2d(save_path=os.path.join(output_dir, "pca_2d.png"), group_order=group_order)
        if self.scores.shape[1] >= 3:
            self.plot_3d(save_path=os.path.join(output_dir, "pca_3d.png"), group_order=group_order)
        self.save_loadings(os.path.join(output_dir, "loadings.csv"))
        self.save_scores(os.path.join(output_dir, "factor_scores.csv"))
        self.save_factor_analysis_result(os.path.join(output_dir, "factor_analysis_result.txt"))
        print("\n所有分析完成！")


if __name__ == "__main__":
    # 使用示例：自定义图例顺序（假设分组名称为 Control, Treatment1, Treatment2）
    analyzer = PCAAnalyzer("1.xlsx", sheet_name="Sheet1", index_col=0)
    analyzer.run_full_analysis(
        output_dir=".",
        standardize=False,               # 不标准化
        group_order=['Control', 'Treatment1', 'Treatment2']  # 自定义顺序
    )
    # 若所有分组均列出，则图例完全按此顺序；若有遗漏，会自动补在末尾并给出提示。