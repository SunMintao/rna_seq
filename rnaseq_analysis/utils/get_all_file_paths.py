from pathlib import Path

def get_all_file_paths(folders, recursive=True):
    """
    获取指定文件夹下所有文件的路径
    :param folders: 文件夹路径列表
    :param recursive: 是否递归子文件夹，默认 True
    :return: 文件路径列表
    """
    all_files = []
    for folder in folders:
        path = Path(folder)
        if not path.exists():
            print(f"警告：{folder} 不存在")
            continue
        if not path.is_dir():
            print(f"警告：{folder} 不是文件夹")
            continue

        if recursive:
            # 递归遍历所有子文件夹
            for item in path.rglob('*'):
                if item.is_file():
                    all_files.append(str(item.resolve()))  # resolve() 获取绝对路径
        else:
            # 只遍历当前文件夹（不递归）
            for item in path.iterdir():
                if item.is_file():
                    all_files.append(str(item.resolve()))
    return all_files

if __name__ == "__main__":
    # 你的四个文件夹路径（请替换为实际路径）
    folders = [
        r"C:\your\folder1",
        r"C:\your\folder2",
        r"C:\your\folder3",
        r"C:\your\folder4"
    ]

    files = get_all_file_paths(folders, recursive=True)  # 递归获取所有文件
    # 或者只获取直接文件：get_all_file_paths(folders, recursive=False)

    # 打印结果
    for f in files:
        print(f)

    # 可选：保存到文本文件
    with open("all_files.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(files))