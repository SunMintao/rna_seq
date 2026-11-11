#!/usr/bin/env python3
"""
获取所有 KEGG 参考通路（ko 编号）的 ENTRY、NAME、DESCRIPTION
输出为 Excel 文件
依赖：requests, pandas, openpyxl
"""

import requests
import time
import sys
import os
import re
from typing import List, Dict, Tuple
import pandas as pd

# ---------- 配置 ----------
KEGG_LIST_URL = "https://rest.kegg.jp/list/pathway/ko"
KEGG_GET_URL = "https://rest.kegg.jp/get"
BATCH_SIZE = 10
REQUEST_DELAY = 0.5


# ---------- 清理函数 ----------
def clean_excel_string(text: str) -> str:
    """
    移除 Excel 不支持的字符（控制字符），保留换行、制表、回车
    """
    if not isinstance(text, str):
        return str(text) if text is not None else ""

    # 移除 ASCII 控制字符（0-31），但保留 \n (10), \r (13), \t (9)
    # 使用正则：保留 \n \r \t，移除其他控制字符
    # 方法：将除 \n \r \t 以外的控制字符替换为空格
    def is_allowed_char(ch):
        code = ord(ch)
        if code == 9 or code == 10 or code == 13:  # \t \n \r
            return True
        if code < 32:  # 其他控制字符
            return False
        return True

    cleaned = ''.join(ch for ch in text if is_allowed_char(ch))
    return cleaned


# ---------- 1. 获取所有 ko 通路 ID ----------
def fetch_ko_pathway_ids() -> List[str]:
    """获取所有 ko 通路 ID（如 ko00010）"""
    print("正在获取 ko 通路列表...")
    resp = requests.get(KEGG_LIST_URL, timeout=30)
    resp.raise_for_status()

    ko_ids = []
    for line in resp.text.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 1:
            continue
        raw_id = parts[0].strip()
        if raw_id.startswith('pathway:ko'):
            ko_id = raw_id[8:]
        elif raw_id.startswith('ko'):
            ko_id = raw_id
        else:
            continue
        ko_ids.append(ko_id)

    print(f"共找到 {len(ko_ids)} 个 ko 通路")
    return ko_ids


# ---------- 2. 批量获取通路详细信息 ----------
def fetch_pathway_details(ko_ids: List[str]) -> Dict[str, str]:
    """批量获取通路详情，返回字典 {ko_id: full_entry_text}"""
    details = {}
    total = len(ko_ids)

    for i in range(0, total, BATCH_SIZE):
        batch = ko_ids[i:i + BATCH_SIZE]
        ids_str = '+'.join(batch)
        url = f"{KEGG_GET_URL}/{ids_str}"
        print(f"正在获取第 {i // BATCH_SIZE + 1} 批: {ids_str}")

        resp = requests.get(url, timeout=30)
        resp.raise_for_status()

        entries = resp.text.strip().split('///')
        for entry in entries:
            if not entry.strip():
                continue
            lines = entry.strip().split('\n')
            if not lines:
                continue
            first_line = lines[0]
            if first_line.startswith('ENTRY'):
                parts = first_line.split()
                if len(parts) >= 2:
                    ko_id = parts[1]
                    details[ko_id] = entry.strip()
        time.sleep(REQUEST_DELAY)

    return details


# ---------- 3. 解析 NAME 和 DESCRIPTION ----------
def parse_kegg_entry(entry_text: str) -> Tuple[str, str]:
    """从 KEGG 条目文本中提取 NAME 和 DESCRIPTION"""
    lines = entry_text.split('\n')
    name = ""
    description = ""
    in_description = False
    desc_lines = []

    for line in lines:
        line_stripped = line.strip()
        if line_stripped.startswith('NAME'):
            parts = line_stripped.split(None, 1)
            if len(parts) > 1:
                name = parts[1].strip()
        elif line_stripped.startswith('DESCRIPTION'):
            parts = line_stripped.split(None, 1)
            if len(parts) > 1:
                desc_lines.append(parts[1].strip())
            in_description = True
        elif in_description:
            if line_stripped and not line_stripped[0].isupper():
                desc_lines.append(line_stripped)
            else:
                in_description = False

    description = ' '.join(desc_lines) if desc_lines else ""
    return name, description


# ---------- 4. 主执行函数 ----------
def run(output_file: str = "ko_pathway_descriptions.xlsx"):
    """获取所有 ko 通路信息并保存为 Excel 文件"""
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"已创建输出目录: {output_dir}")

    print("=== KEGG ko 通路 ENTRY, NAME, DESCRIPTION 获取工具 ===")
    print(f"输出文件: {output_file}\n")

    ko_ids = fetch_ko_pathway_ids()
    if not ko_ids:
        print("未找到任何 ko 通路，请检查网络或 API")
        sys.exit(1)

    details = fetch_pathway_details(ko_ids)
    print(f"成功获取 {len(details)} 个通路的详细信息")

    data = []
    for ko_id in ko_ids:
        entry_text = details.get(ko_id)
        if not entry_text:
            continue
        name, desc = parse_kegg_entry(entry_text)
        # 清理 Excel 不支持的控制字符
        name = clean_excel_string(name)
        desc = clean_excel_string(desc)
        data.append({
            "Pathway_id": ko_id,
            "Pathway": name,
            "DESCRIPTION": desc
        })

    df = pd.DataFrame(data)

    try:
        df.to_excel(output_file, index=False, engine='openpyxl')
        print(f"\n结果已保存至: {output_file}")
        print(f"共写入 {len(df)} 条通路记录")
    except Exception as e:
        print(f"保存 Excel 文件失败: {e}")
        sys.exit(1)


# ---------- 5. 直接运行（在此修改输出路径） ----------
if __name__ == "__main__":
    # 修改为你需要的路径，例如：
    # run("/home/user/data/my_ko.xlsx")
    run()  # 默认当前目录