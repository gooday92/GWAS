import pandas as pd
from scipy import stats
import time

# 先合并表型和基因的表格
epi_file = r"D:\pythonProject\GWAS_selfmade\表葡\统计\原始文件\耐药表型.csv"
gene_pre_ab = r"D:\pythonProject\GWAS_selfmade\表葡\统计\原始文件\gene_presence_absence.csv"
gene_type_matrix = r"D:\pythonProject\GWAS_selfmade\表葡\统计\原始文件\gene_type_matrix.csv"
df_epi = pd.read_csv(epi_file)
epi_list = df_epi.columns.to_list()[1:]

# 存在与否与耐药表型的关系
df_pre_ab = pd.read_csv(gene_pre_ab)
gene_pre_list = df_pre_ab.columns.to_list()[1:]
df_mer_pre_ab = pd.merge(df_pre_ab, df_epi, how="left", on="strain")
dic_result = dict()
for drug in epi_list:
    print(drug)
    dic_result[drug] = dict()
    for gene in gene_pre_list:
        contingency_table = pd.crosstab(df_mer_pre_ab[gene], df_mer_pre_ab[drug], margins=True)
        # row = contingency_table.shape[0] - 1
        # col = contingency_table.shape[1] - 1
        stat, p, dof, expected = stats.chi2_contingency(contingency_table)
        dic_result[drug][gene] = p
        if drug == "OXA-RS" and gene == "mecA2":
            print(contingency_table)
            print(stat, p, dof, expected)
df_result = pd.DataFrame(dic_result)
df_result.to_csv("epi_pre_ab_v2.csv")

# 基因型别与耐药表型的关系
# df_gene_type = pd.read_csv(gene_type_matrix)
# gene_pre_list = df_gene_type.columns.to_list()[1:]
# df_mer_pre_ab = pd.merge(df_gene_type, df_epi, how="left", on="strain")
# dic_result = dict()
# for drug in epi_list:
#     print(drug)
#     dic_result[drug] = dict()
#     for gene in gene_pre_list:
#         contingency_table = pd.crosstab(df_mer_pre_ab[gene], df_mer_pre_ab[drug], margins=True)
#         # row = contingency_table.shape[0] - 1
#         # col = contingency_table.shape[1] - 1
#         stat, p, dof, expected = stats.chi2_contingency(contingency_table)
#         dic_result[drug][gene] = p
# df_result = pd.DataFrame(dic_result)
# df_result.to_csv("epi_types.csv")
