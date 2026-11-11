import pandas as pd
import re

'''
注意：传入的 两个文件 必须 与 结果文件的列名一样（gene.pathway_template.xlsx）
'''
def parse_ko_string(ko_str):
    """
    从 KEGG_ORTHOLOGY 列中提取 KO 编号（去掉 'ko:' 前缀）。
    输入如 "ko:K00033" 或 "ko:K00033;ko:K00001" 或 "_" 或空。
    返回列表。
    """
    if pd.isna(ko_str) or ko_str == '_' or ko_str.strip() == '':
        return []
    # 用分号分隔
    parts = ko_str.split(';')
    kos = []
    for p in parts:
        p = p.strip()
        if p.startswith('ko:'):
            # 去掉 ko: 前缀
            kos.append(p[3:])
        else:
            # 如果没有前缀，可能直接是 K 号？但按照样本是带 ko: 的
            kos.append(p)
    return kos

def parse_pathway_string(pw_str):
    """
    从 Pathway_id 列中提取通路编号（去掉 'path:' 前缀，取 ko 号）。
    输入如 "path:ko00030;path:ko00480;..." 或空。
    返回列表。
    """
    if pd.isna(pw_str) or pw_str.strip() == '':
        return []
    parts = pw_str.split(';')
    pws = []
    for p in parts:
        p = p.strip()
        if p.startswith('path:'):
            # 提取 ko 号，如 "ko00030"
            pws.append(p[5:])
        else:
            pws.append(p)
    return pws

def generate_gene_pathway(ann_file, map_file, output_file):
    """
    生成 gene.pathway.xlsx
    """
    # 1. 读取注释文件
    df_ann = pd.read_excel(ann_file)
    # 确保列名存在
    required_cols = ['Gene_id', 'KEGG_ORTHOLOGY', 'Pathway_id']
    for col in required_cols:
        if col not in df_ann.columns:
            raise ValueError(f"注释文件缺少列: {col}")

    # 构建 gene -> list of KOs, gene -> list of pathways
    gene_to_kos = {}
    gene_to_pathways = {}

    for _, row in df_ann.iterrows():
        gene = row['Gene_id']
        if pd.isna(gene) or gene == '':
            continue
        kos = parse_ko_string(row['KEGG_ORTHOLOGY'])
        pws = parse_pathway_string(row['Pathway_id'])

        if kos:
            gene_to_kos.setdefault(gene, []).extend(kos)
        if pws:
            gene_to_pathways.setdefault(gene, []).extend(pws)

    # 2. 读取通路描述文件
    df_map = pd.read_excel(map_file)
    pathway_to_name = {}
    for _, row in df_map.iterrows():
        pid = row['Pathway_id']
        name = row['Pathway']
        if pd.notna(pid) and pd.notna(name):
            pathway_to_name[pid] = name

    # 3. 构建 pathway -> set of genes, pathway -> set of KOs
    pathway_genes = {}
    pathway_kos = {}

    for gene, pws in gene_to_pathways.items():
        for pw in pws:
            if pw not in pathway_genes:
                pathway_genes[pw] = set()
                pathway_kos[pw] = set()
            pathway_genes[pw].add(gene)
            # 添加该基因的所有KO（如果该基因有KO）
            if gene in gene_to_kos:
                pathway_kos[pw].update(gene_to_kos[gene])

    # 4. 生成输出 DataFrame
    rows = []
    for pw, genes in pathway_genes.items():
        name = pathway_to_name.get(pw, pw)  # 若找不到名称，用ID代替
        gene_list = sorted(genes)
        gene_str = ';'.join(gene_list)
        ko_list = sorted(pathway_kos.get(pw, set()))
        ko_str = '+'.join(ko_list)
        rows.append({
            'Pathway': name,
            'Pathway_id': pw,
            'Gene_id': gene_str,
            'KEGG_ORTHOLOGY': ko_str
        })

    df_out = pd.DataFrame(rows)
    # 按通路名称排序
    df_out = df_out.sort_values('Pathway').reset_index(drop=True)

    # 5. 保存为 Excel
    df_out.to_excel(output_file, index=False)
    print(f"已生成 {output_file}，共 {len(df_out)} 个通路。")

if __name__ == '__main__':
    generate_gene_pathway(
        ann_file='kegg_annotation.xlsx',
        map_file='ko_pathway_descriptions.xlsx',
        output_file='gene.pathway.xlsx'
    )