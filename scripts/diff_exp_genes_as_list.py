import csv
import json

def diff_exp_genes_as_list(inpath: str, outpath: str):
    genes = set()
    with open(inpath, 'r') as handler:
        reader = csv.DictReader(handler)
        for row in reader:
            genes.update(row["Proteins in dataset"].split(','))
    with open(outpath, 'w') as handler:
        handler.write("\n".join(genes))




def diff_gene_list_to_entrez(inpath, outpath, translation_path):
    with open(translation_path, 'r') as handler:
        tran_dict = json.load(handler)
    with open(inpath, 'r') as handler:
        translations = [tran_dict[g.strip()] for g in handler.readlines()]
    with open(outpath, 'w') as handler:
        handler.write("\n".join(translations))






if __name__ == "__main__":
    diff_gene_list_to_entrez(r"D:\differential_gene_expression\https_doi.org_10_1038_s41586_020_2332_7\upregulated_genes_24h.txt",
                             r"D:\differential_gene_expression\https_doi.org_10_1038_s41586_020_2332_7\upregulated_genes_24h_entrez.txt",
                             r"D:\complete_translation.json")