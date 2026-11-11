"""
多组 Venn 图绘制工具（面向对象版本，已修复颜色参数为 None 的错误）
支持 2-6 组基因 ID 的 Venn 图绘制。
依赖库：matplotlib, matplotlib-venn (2-3组), venn (4-6组)
"""

from typing import Dict, List, Set, Union, Optional, Tuple
import matplotlib.pyplot as plt


class MultiVenn:
    """多组 Venn 图绘制器，支持 2-6 组"""

    def __init__(
        self,
        groups: Dict[str, Union[List[str], Set[str]]],
        colors: Optional[List[str]] = None,
        title: Optional[str] = None,
        figsize: Tuple[int, int] = (8, 8),
    ):
        self.groups = groups
        self.num_groups = len(groups)
        if self.num_groups < 2:
            raise ValueError("至少需要提供2组数据")
        if self.num_groups > 6:
            raise ValueError("本类最多支持6组 Venn 图")

        self.sets = {name: set(ids) for name, ids in groups.items()}
        self.names = list(groups.keys())
        self.colors = colors
        self.title = title
        self.figsize = figsize

        self._import_libs()
        self.figure = None
        self.axes = None

    def _import_libs(self):
        if self.num_groups <= 3:
            try:
                from matplotlib_venn import venn2, venn3
                self.venn2 = venn2
                self.venn3 = venn3
                self.venn_lib = None
            except ImportError:
                raise ImportError("2-3组 Venn 图需要 matplotlib-venn，请安装：pip install matplotlib-venn")
        else:
            try:
                import venn as venn_lib
                self.venn_lib = venn_lib
                self.venn2 = self.venn3 = None
            except ImportError:
                raise ImportError("4-6组 Venn 图需要 venn 库，请安装：pip install venn")

    def draw(self):
        plt.figure(figsize=self.figsize)

        if self.num_groups == 2:
            if self.venn2 is None:
                raise RuntimeError("matplotlib_venn 未正确导入")
            set_a, set_b = self.sets[self.names[0]], self.sets[self.names[1]]
            # 动态构建参数，避免传递 None
            kwargs = {"set_labels": self.names}
            if self.colors is not None:
                kwargs["set_colors"] = self.colors
            self.venn2([set_a, set_b], **kwargs)

        elif self.num_groups == 3:
            if self.venn3 is None:
                raise RuntimeError("matplotlib_venn 未正确导入")
            set_a, set_b, set_c = self.sets[self.names[0]], self.sets[self.names[1]], self.sets[self.names[2]]
            kwargs = {"set_labels": self.names}
            if self.colors is not None:
                kwargs["set_colors"] = self.colors
            self.venn3([set_a, set_b, set_c], **kwargs)

        else:  # 4-6 组
            if self.venn_lib is None:
                raise RuntimeError("venn 库未正确导入")
            set_list = [self.sets[name] for name in self.names]
            labels = self.venn_lib.get_labels(set_list, fill=["number"])
            # venn 库的绘图函数接受 colors 参数，同样需要处理 None
            kwargs = {"names": self.names}
            if self.colors is not None:
                kwargs["colors"] = self.colors
            if self.num_groups == 4:
                fig, ax = self.venn_lib.venn4(labels, **kwargs)
            elif self.num_groups == 5:
                fig, ax = self.venn_lib.venn5(labels, **kwargs)
            elif self.num_groups == 6:
                fig, ax = self.venn_lib.venn6(labels, **kwargs)
            else:
                raise ValueError(f"不支持的组数: {self.num_groups}")
            # 重新关联当前 figure
            plt.figure(fig.number)

        if self.title:
            plt.title(self.title)

        self.figure = plt.gcf()
        self.axes = plt.gca()

    def save(self, output_path: str, dpi: int = 300):
        if self.figure is None:
            self.draw()
        self.figure.savefig(output_path, bbox_inches="tight", dpi=dpi)
        print(f"Venn 图已保存至: {output_path}")

    def show(self):
        if self.figure is None:
            self.draw()
        plt.show()

    def get_figure(self):
        if self.figure is None:
            self.draw()
        return self.figure


# 使用示例
if __name__ == "__main__":
    # 3组数据
    data_3 = {
        "Group1": ["geneA", "geneB", "geneC", "geneD"],
        "Group2": ["geneC", "geneD", "geneE", "geneF"],
        "Group3": ["geneD", "geneE", "geneG", "geneH"],
    }
    venn3 = MultiVenn(data_3, title="3组 Venn 图")
    venn3.save("venn3.png")

    # 4组数据（需要安装 venn 库）
    data_4 = {
        "A": ["gene1", "gene2", "gene3", "gene4"],
        "B": ["gene2", "gene3", "gene5", "gene6"],
        "C": ["gene3", "gene4", "gene5", "gene7"],
        "D": ["gene4", "gene5", "gene6", "gene8"],
    }
    venn4 = MultiVenn(data_4, title="4组 Venn 图")
    venn4.save("venn4.png")