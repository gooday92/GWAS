import pandas as pd
import os
from Bio import AlignIO
from icecream import ic


def align_df(path):
    # 输入alignment文件，输出df
    alignment = AlignIO.read(path, "fasta")
    df_dict = {"ID": [], "Sequence": []}
    for record in alignment:
        df_dict["ID"].append(record.id.split(";")[0])
        df_dict["Sequence"].append("".join(record.seq))
    df = pd.DataFrame(df_dict)
    return df

# 读取菌株名
strain_list = []
with open("strains_list.txt") as file_strain_list:
    strain_list = file_strain_list.read().splitlines()

# 构建结果数据
dict_result = {"gene": []}
for strain in strain_list:
    dict_result[strain] = []
# print(dict_result)
# 读取每个gene_alignment的矩阵
path_alns = r"D:\pythonProject\GWAS_selfmade\表葡\panaroo\aligned_gene_sequences"
for align in os.listdir(path_alns):
    gene_name = align.split(".")[0]
    dict_result["gene"].append(gene_name)
    # print(gene_name)
    file_aln = os.path.join(path_alns, align)
    df_aln = align_df(file_aln)
    # print(df_aln.groupby)
    grouped = df_aln.groupby(["Sequence"])
    n = 1
    for seqs in grouped:
        for strain_n in seqs[1].ID:
            dict_result[strain_n].append(n)
        n += 1
    # 每次都检验结果中是不是每一个菌株都加入了结果，是不是同样长度
    for key, item in dict_result.items():
        # 如果改基因组内没有这个基因，则加0
        if len(item) < len(dict_result["gene"]):
            dict_result[key].append(0)
        # 如果该基因组内存在双拷贝，则减一
        if len(item) > len(dict_result["gene"]):
            del(dict_result[key][-1])

for k, i in dict_result.items():
    if len(i) != 249:
        print(k, print(len(i)))
df_result = pd.DataFrame(dict_result)
df_result.to_csv("gene_type_matrix.csv")

