from Bio import SeqIO
import os

# file_gene_list = r"D:\pythonProject\GWAS_selfmade\表葡\3.结果分析\gene_list_oxa.txt"
file_gene_reference = r"D:\pythonProject\GWAS_selfmade\表葡\3.结果分析\pan_genome_reference.fa"
# out_put_path = r"D:\pythonProject\GWAS_selfmade\表葡\3.结果分析\new_seq\candicate_genes_oxa.fa"

files_list_path = r"D:\pythonProject\GWAS_selfmade\表葡\3.结果分析\candidate_gene_lists"
file_output_path = r"D:\pythonProject\GWAS_selfmade\表葡\3.结果分析\new_seq"
for f in os.listdir(files_list_path):
    file_input = os.path.join(files_list_path, f)
    print(file_input)
    drug = file_input.split("\\")[-1].split(".")[0].replace("gene_list_", "")
    file_out = f"candidate_genes_{drug}.fa"
    file_output = os.path.join(file_output_path, file_out)
    # 输入基因列表，从pan-genome-reference中抽取序列构成新的fasta文件
    # 读取基因列表
    with open(file_input) as file:
        gene_list = file.read().splitlines()

    # 读取gene_reference
    genes_records = list(SeqIO.parse(file_gene_reference, "fasta"))

    # 筛选后保存
    new_records = []
    for record in genes_records:
        if record.id in gene_list:
            new_records.append(record)

    SeqIO.write(new_records, file_output, "fasta")
