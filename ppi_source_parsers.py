import re

from utils.queue_managers import load_json, dump_json
import biomart
import pandas as pd
from typing import Union


def translate_table_to_entrez(table_path: str,
                              cols_to_copy: Union[str, list[str]], cols_to_translate: Union[str, list[str]],
                              trans_dict_path: str, output_table_path: str, use_index: bool=False):
    if isinstance(cols_to_copy, str):
        cols_to_copy = [cols_to_copy]
    if isinstance(cols_to_translate, str):
        cols_to_translate = [cols_to_translate]

    table = pd.read_csv(table_path, usecols=cols_to_copy + cols_to_translate)
    trans_dict = load_json(trans_dict_path)
    for col in cols_to_translate:
        table[col] = table[col].apply(lambda prot_name: trans_dict.get(prot_name.upper(), prot_name))

    table.to_csv(output_table_path, index=use_index)


def query_ensembl():
    server = biomart.BiomartServer("http://uswest.ensembl.org/biomart")
    db = server.databases["ENSEMBL_MART_ENSEMBL"]
    ds = db.datasets["hsapiens_gene_ensembl"]
    attributes = ("ensembl_gene_id", "protein_id", "entrezgene_id")
    res = ds.search({"attributes": attributes})
    with open("ens_res.txt", "w") as handler:
        handler.write(res.text) # "ensembl_gene_id\tentrezgene_id\n" +
    # df = pd.read_csv("ens_res.txt")
    # df.to_csv("ens_res_with_protein_id.csv", index=False, sep=",")
    # translate_table_to_entrez("ens_res.csv", "ensembl_gene_id", "entrezgene_id", r"D:\data\other\symbol_to_entrezgene.json",
    #                           "ens_res_translated.csv")








# Returns a list of tuples of the form (virus_protein_name, human_protein_entrez)
def extract_human_uniprot_ids(mitab_file_path):
    mitab_df = pd.read_csv(mitab_file_path)
    species_A_col = "Taxid interactor A"
    prot_A_col = "# ID(s) interactor A"
    prot_B_col = "ID(s) interactor B"

    human_uniprots = []
    for i, row in mitab_df.iterrows():
        species_A = re.search("\d+", row[species_A_col]).group()
        species_A = "human" if species_A == "9606" else "other"
        human_prot_id_col = prot_A_col if species_A == "human" else prot_B_col
        human_prot_id = row[human_prot_id_col].split(":")[-1]
        human_uniprots.append(human_prot_id)

    return human_uniprots






# Returns a list of tuples of the form (virus_protein_name, human_protein_entrez)
def parse_mitab_interactions(mitab_file_path, trans_dict_parth):
    mitab_df = pd.read_csv(mitab_file_path)
    trans_dict = load_json(trans_dict_parth)
    species_A_col = "Taxid interactor A"
    # species_B_col = "Taxid interactor B"
    prot_A_col = "# ID(s) interactor A"
    prot_B_col = "ID(s) interactor B"

    ppis = []
    for i, row in mitab_df.iterrows():
        species_A = re.search("\d+", row[species_A_col]).group()
        species_A = "human" if species_A == "9606" else "other"
        human_prot_id_col = prot_A_col if species_A == "human" else prot_B_col
        other_prot_id_col = prot_A_col if human_prot_id_col == prot_B_col else prot_B_col
        # human_prot_id = re.search("ensembl:.+?\|", row[human_prot_id_col]).group().split(":")[-1][:-1]
        # translated_human_prot_id = trans_dict.get(human_prot_id, None)
        # other_prot_id = re.search("ensembl:.+?\|", row[human_prot_id_col]).group().split(":")[-1][:-1]
        human_prot_id = row[human_prot_id_col].split(":")[-1]
        human_prot_id_translated = trans_dict.get(human_prot_id, None)
        other_prot_id = row[other_prot_id_col].split(":")[-1]
        ppis.append((other_prot_id, human_prot_id_translated))

        # ppis.append((other_prot_id, translated_human_prot_id))

    return ppis


def add_uniprot_translations(uniprot_to_symbol: str, trans_dict_file: str, new_trans_dict_file: str):
    df = pd.read_csv(uniprot_to_symbol)
    trans_dict = load_json(trans_dict_file)

    for _, row in df.iterrows():
        uniprot = row["From"]
        symbol = row["To"]
        if symbol not in trans_dict:
            print(symbol)
        else:
            trans_dict[uniprot] = trans_dict[symbol]

    dump_json(trans_dict, new_trans_dict_file)

def main():
    input_tuples = [
        (r"D:\data\networks\virus_to_human\pr81av_host_protein_interactome.csv", "Viral Gene", "Gene ID",
         r"D:\data\other\symbol_to_entrezgene.json", r"D:\data\networks\virus_to_human\pr81av_host_protein_interactome_translated.csv")
    ]
    df = pd.read_csv(r"D:\data\networks\virus_to_human\intact_IM-21137.txt", delimiter="\t")
    x = 7
    for tup in input_tuples:
        translate_table_to_entrez(*tup)




if __name__ == "__main__":
    # add_uniprot_translations(r"D:\data\other\uniprot_to_symbol\uniprot-compressed_true_download_true_format_xlsx-2022.06.30-13.00.35.71.csv",
    #                          r"D:\data\other\symbol_to_entrezgene.json", r"D:\data\other\symbol_to_entrezgene.json")
    uniprots = parse_mitab_interactions(r"D:\data\networks\virus_to_human\intact_IM-21137.csv",
                            r"D:\data\other\symbol_to_entrezgene.json")
    viral_gene_col = [u[0] for u in uniprots]
    human_entrez_col = [u[1] for u in uniprots]
    df = pd.DataFrame({"Viral Gene": viral_gene_col, "Gene ID": human_entrez_col})
    df.to_csv(r"D:\data\networks\virus_to_human\hep_c_host_protein_interactome_translated.csv", index=False)

    # with open("uniprots_to_query.txt", "w") as handler:
    #     handler.write("\n".join(uniprots))
    # query_ensembl()
    # with open("ens_res.txt", 'r') as handler:
    #     lines = handler.readlines()
    # # txt = "".join([l.replace("\t\n", "\t\t\n") for l in lines])
    # # with open("ens_res_fixed.txt", 'w') as handler:
    # #     handler.write(txt)
    # for i, l in enumerate(lines):
    #     if len(l.split("\t")) != 2:
    #         lines[i] = l[:-1] + "\t\n"
    #     lines[i] = l[:-1]
    #
    # for i, l in enumerate(lines):
    #     if len(l.split("\t")) != 2:
    #         lines[i] = l[:-1] + "\t\n"
    # ens_col = [l.split("\t")[0] for l in lines]
    # entrez_col = [l.split("\t")[1] for l in lines]
    # df = pd.DataFrame({"ens": ens_col, "entrez": entrez_col})
    # df = df.set_index("ens")
    # df = df[df["entrez"] != ""]
    # df = df[df["entrez"] != "\n"]
    # s = df["entrez"]
    #
    # # df = pd.read_csv(r"ens_res_fixed.txt", delimiter="\t", index_col="ens")
    # d = s.to_dict()
    # existing_trans = load_json(r"D:\data\other\symbol_to_entrezgene.json")
    # for k, v in d.items():
    #     k = k.upper()
    #     if k not in existing_trans:
    #         existing_trans[k] = v
    # dump_json(existing_trans, r"D:\data\other\symbol_to_entrezgene.json")





    # with open(r"D:\downloads\2022-06-29-13-00.txt", 'r') as handler:
    #     lines = handler.readlines()
    # x = 7
    # df = pd.read_csv(r"D:\data\networks\virus_to_human\biomart_res.csv", usecols=["entrezgene_id", "uniprot_gn_symbol"],
    #                  index_col="uniprot_gn_symbol")
    # s = df["entrezgene_id"]
    # s = s[s != "none"].drop(index="none")
    # d = s.to_dict()
    # dump_json(d, "new_uniprot_to_entrez.json")
    # biomart = load_json(r"D:\data\other\biomart_uniprot_to_entrez.json")
    # existing_trans = load_json(r"D:\data\other\symbol_to_entrezgene.json")
    # for k, v in biomart.items():
    #     k = k.upper()
    #     if k not in existing_trans:
    #         existing_trans[k] = v
    # dump_json(existing_trans, r"D:\data\other\symbol_to_entrezgene.json")

    # x = 7
    # with open(r"D:\data\networks\virus_to_human\biomart_res.txt", 'r') as handler:
    #     lines = handler.readlines()
    # as_str = "\n".join([",".join(line.replace("\t\t", " none ").replace("\t\n", " none").split()) for line in lines])
    # with open(r"D:\data\networks\virus_to_human\biomart_res.csv", 'w') as handler:
    #     handler.write(as_str)
    # main()
