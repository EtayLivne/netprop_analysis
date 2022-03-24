import pandas as pd
"""
load_metadata_rna: to be used with "D:\\data\\networks\\metastudy\\rna_interactome_translated.csv"
load_metadata_protein_interactions: to be used with "D:\\data\\networks\\metastudy\\protein_interactome_translated.csv"
"""

def load_metadata_rna(file_path, approved_cell_lines):
    # if not self.rna_interactions_path:
    #     return []
    id_col = "entrez"
    df = pd.read_csv(file_path,
                     usecols=["Assay cell line", id_col])
    df[id_col].astype("int32")
    df.set_index("Assay cell line", inplace=True)
    interactions = []
    for cell_line in set(df.index):
        if cell_line not in approved_cell_lines:
            continue
        interactions.extend([("covid_rna", row[id_col]) for _, row in df.loc[cell_line].iterrows()])

    return interactions


def load_metadata_protein_interactions(file_path, approved_cell_lines, merge=False):
    # if not self.protein_interactions_path:
    #     return []

    id_col = "entrez"
    df = pd.read_csv(file_path,
                     usecols=["Assay cell line", "Bait", "Reference", id_col], index_col="Assay cell line")
    df.dropna(inplace=True)
    df[id_col] = df[id_col].astype("int32")
    interactions = []

    validation_counter = dict()
    for cell_line in set(df.index):
        if cell_line not in approved_cell_lines:
            continue
        for _, row in df.loc[cell_line].iterrows():
            key = (row["Bait"], row[id_col])
            validation_counter[key] = validation_counter.get(key, 0) + 1

        interactions.extend([k for k in validation_counter if validation_counter[k] > 1])

    if merge:
        interactions = [("covid_proteins", x[1]) for x in interactions]

    return interactions