#!/usr/bin/env python3
"""
create_rnaseq_structure.py

Create empty directory structure and empty Python files for the rnaseq_analysis package.
Run this script in the desired parent directory.
"""

import os

PROJECT_NAME = "rnaseq_analysis"

# 需要创建的目录
DIRECTORIES = [
    PROJECT_NAME,
    os.path.join(PROJECT_NAME, "core"),
    os.path.join(PROJECT_NAME, "visualization"),
    os.path.join(PROJECT_NAME, "utils"),
    "tests",
    "examples",
    os.path.join("docs", "source"),
]

# 需要创建的空文件（相对路径）
EMPTY_FILES = [
    # 根目录
    "README.md",
    "requirements.txt",
    "setup.py",
    # 主包
    os.path.join(PROJECT_NAME, "__init__.py"),
    # core 模块
    os.path.join(PROJECT_NAME, "core", "__init__.py"),
    os.path.join(PROJECT_NAME, "core", "data_loader_2.py"),
    os.path.join(PROJECT_NAME, "core", "preprocess_3.py"),
    os.path.join(PROJECT_NAME, "core", "normalize.py"),
    os.path.join(PROJECT_NAME, "core", "diff_expr_4.py"),
    os.path.join(PROJECT_NAME, "core", "enrichment.py"),
    # visualization 模块
    os.path.join(PROJECT_NAME, "visualization", "__init__.py"),
    os.path.join(PROJECT_NAME, "visualization", "pca.py"),
    os.path.join(PROJECT_NAME, "visualization", "volcano.py"),
    os.path.join(PROJECT_NAME, "visualization", "heatmap.py"),
    os.path.join(PROJECT_NAME, "visualization", "expression.py"),
    # utils 模块
    os.path.join(PROJECT_NAME, "utils", "__init__.py"),
    os.path.join(PROJECT_NAME, "utils", "io.py"),
    os.path.join(PROJECT_NAME, "utils", "pca_analysis.py"),
    os.path.join(PROJECT_NAME, "utils", "plot_helpers.py"),
    # pipeline
    os.path.join(PROJECT_NAME, "pipeline.py"),
    # tests
    os.path.join("tests", "__init__.py"),
    os.path.join("tests", "test_data_loader.py"),
    # examples (可以放示例数据，但用户要求空文件，只创建占位文件)
    os.path.join("examples", "counts_example.csv"),
    os.path.join("examples", "metadata_example.csv"),
    os.path.join("examples", "run_analysis.ipynb"),
]

def create_structure():
    base = os.getcwd()
    print(f"Creating project structure in: {base}")

    # 创建目录
    for d in DIRECTORIES:
        path = os.path.join(base, d)
        os.makedirs(path, exist_ok=True)
        print(f"Created directory: {d}")

    # 创建空文件
    for f in EMPTY_FILES:
        path = os.path.join(base, f)
        # 确保父目录存在（理论上已创建）
        parent = os.path.dirname(path)
        os.makedirs(parent, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fp:
            # 如果是 Python 文件，可以添加简单的注释头（可选），但用户要求空文件，就留空
            # 如果希望完全空文件，直接 pass
            pass
        print(f"Created empty file: {f}")

    print("\n" + "=" * 60)
    print(f"Project '{PROJECT_NAME}' structure created successfully!")
    print("Next steps:")
    print(f"  cd {PROJECT_NAME}")
    print("  pip install -e .")
    print("=" * 60)

if __name__ == "__main__":
    create_structure()