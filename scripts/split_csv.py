import pandas as pd
from utils.utils import load_json, dump_json

def split_csv(file_path: str, split_col: str):
    df = pd.read_csv(file_path, index_col=split_col)
    df = df[df.index.notnull()]
    return pd.Series({ind: df.loc[ind].reset_index() for ind in set(df.index)})


def to_entrez(csv_path, outpath, trans_dict):
    trans = load_json(trans_dict)
    trans_func = lambda x: trans.get(x, "etay")

    df = pd.read_csv(csv_path)
    df["Entrez"] = df["Gene name"].apply(trans_func)
    failed_translations = df[df["Entrez"] == "etay"]
    failed = list(failed_translations["Gene name"])
    print(f"failed to translate the following genes: {failed}")
    df.to_csv(outpath)


def extract_uniprot(csv_path, outpath):
    df = pd.read_csv(csv_path, index_col="Entrez")
    problematics = set(df.loc["etay"]["Gene name"].values)
    with open(outpath, 'w') as handler:
        for p in problematics:
            handler.write(f"{p}\n")


def ncbi_query_strings(input_file, output_file):
    with open(input_file, 'r') as handler:
        gene_names = handler.readlines()
    query_strings = [f"({gene_name.strip()}[gene]) AND (Homo sapiens[orgn])" for gene_name in gene_names]
    with open(output_file, 'w') as handler:
        for qs in query_strings:
            handler.write(f"{qs}\n")


def back_to_entrez(query_strings, naama_exe, naama_output):
    # call naama exe with to_translate and naama_output
    with open(query_strings, 'r') as handler:
        raw_queries = handler.readlines()   # Format: (<what I care about>[gene]) bla bla bla
    query_gene_names = [s.split('[')[0][1:] for s in raw_queries]
    with open(naama_output, 'r') as handler:
        raw_trans = handler.readlines() # format: bla bla: <what I care about>, bla bla

    just_ids = [s.split(",")[0].split()[-1] for s in raw_trans]
    if len(query_gene_names) != len(just_ids):
        raise Exception("Not all queries have been succesfully completed")
    return {query_gene_names[i]: just_ids[i] for i in range(len(query_gene_names)) if just_ids[i] != "found"}


def update_translation_dict(translation_file, new_translations, new_translation_file):
    transaltion_dict = load_json(translation_file)
    existing_keys = set(transaltion_dict.keys()) & set(new_translations.keys())
    if existing_keys:
        if diffs := [k for k in existing_keys if transaltion_dict[k] != new_translations[k]]:
            error_str = "\n".join(diffs)
            raise Exception(f"The following keys have different translation in existing and new dicts: {error_str}")
    transaltion_dict.update(new_translations)
    dump_json(transaltion_dict, new_translation_file)
