from pathlib import Path
from itertools import repeat
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd


def _get_whitelist(df: pd.DataFrame, rank_by: str=None, threshold: int=1000) -> list:    # Returns the indexes of the values not to be discarded
    if rank_by is None:
        print("rankin by nothin")
        df = df.sort_index()
    else:
        try:
            print(f"ranking by {rank_by}")
            df = df.sort_values(rank_by, ascending=True)    #ascending = true because sorting by rank
        except KeyError:
            raise ValueError(f"cannot filter by rank of {rank_by} because df does not have this column")
    return list(df.iloc[:threshold].index)


def _column_concat_average_dist_from_reference(dfs: list[pd.DataFrame], reference_df: pd.DataFrame, protein: str) -> pd.DataFrame:
    avg = "__average__"
    abs_avg = "__abs_average__"
    nodes = set(reference_df.index)
    print(f"*** {protein} with reference df {reference_df.name}***")

    # for df in dfs:
    #     filtered_df_idx = set(df.index)
    #     print(f"{df.name}: {len(filtered_df_idx - nodes)} / {len(filtered_df_idx - nodes)}")
    # indexed_reference = reference_df.set_index("nodes")
    # indexes = set(indexed_reference.index)
    # print(f"*** {protein} with reference df {reference_df.name}***")
    # for df in dfs:
    #     nodes = set(df["nodes"])
    #     print(f"{df.name}: {len(indexes - nodes)} / {len(nodes - indexes)} / {len(df) == 1000} / {1 in nodes}")

    try:
        for df in dfs:
            try:
                """
                Each column in the df are the ranks of the proteins in a specific subset.
                Each row is the ranks of a specific protein across all subsets.
                to calculate avg value of a row, find the rank of the protein it represents in the original df (that is 
                filtered to only have values for the propagation from the covid protein of relevance), subtract it from
                all elements in this row, and calculate the average of that. 
                What you get is the average *offset* of an item in this set from its rank in the original propagation.
                This number might be interesting 
                For each
                """
                df[avg] = df.apply(lambda row: int((row - reference_df.loc[row["nodes"]][protein]).mean()), axis=1)
                # df[abs_avg] = df.apply(lambda row: int((row - reference_df.loc[row["nodes"]][protein]).abs().mean()), axis=1)
            except Exception as e:
                print(df.name)
                print(df.loc)
                print(reference_df.loc[1])
                raise e
        res_dict = {"nodes": dfs[0]["nodes"]}
        res_dict.update({df.name: df[avg] for df in dfs})
        res_df = pd.DataFrame(res_dict)
    finally:
        for df in dfs:
            if avg in df.columns:
                df.drop(avg, axis=1)
    return res_df


def _sort_by_average_diff(df: pd.DataFrame) -> pd.DataFrame:
    AVG = "__avg__"
    df[AVG] = df.apply(lambda row: row[[c for c in row.index if c != "nodes"]].sum(), axis=1)
    df = df.sort_values(AVG)
    df.drop(AVG, axis=1, inplace=True)

    return df


def single_protein_sorted_df(prop_df: pd.DataFrame, knockout_dfs: list[pd.DataFrame], protein: str,
                             top_propagated_threshold: int) -> pd.DataFrame:
    prop_df = prop_df.set_index("nodes",)
    whitelist = _get_whitelist(prop_df, rank_by=protein, threshold=top_propagated_threshold)
    # whitelist_set = set(whitelist)    # For debug only

    filtered_dfs = []   # dfs where only proteins that rank high in original propagation are retained
    for df in knockout_dfs:

        indexed_df = df.set_index("nodes", drop=False)
        filtered_df = indexed_df.loc[whitelist]
        filtered_df.name = df.name

        filtered_dfs.append(filtered_df)
    # TODO IS'nt all this just re-discovering the whitelist? I', like 90% sure that it is [NO IT ISNT BECAUSE THE INDEX IS NOT THE SAME AS THE PROPAGATION RANKING
    node_mask = list(filtered_dfs[0].index)
    idxs = []
    for tup in prop_df.itertuples():
        if tup.Index in node_mask:
            idxs.append(tup.Index)  # the index
    masked_prop_df = prop_df.loc[idxs]
    masked_prop_df.name = protein

    # counter = 0
    # for df in filtered_dfs:
    #     filtered_df_idx = set(df.index)
    #     masked_nodes = set(masked_prop_df.index)
    #     print(f"{df.name}: {len(filtered_df_idx - whitelist_set)} / {len(filtered_df_idx - whitelist_set)} / {len(filtered_df_idx - masked_nodes)} / {len(masked_nodes - filtered_df_idx)}")

    averaged_df = _column_concat_average_dist_from_reference(filtered_dfs, masked_prop_df, protein)
    sorted_df = _sort_by_average_diff(averaged_df)
    return sorted_df


def _worker(tup: tuple[str, Path]) -> pd.DataFrame:
    top_propagated_threshold, prop_res_file, knockout_res_folder = tup
    prop_df = pd.read_csv(prop_res_file)
    prop_df.name = "full_prop"
    knockout_dfs = []
    for file in knockout_res_folder.glob("*"):
        df = pd.read_csv(str(file))
        df.name = file.stem
        knockout_dfs.append(df)
    protein = knockout_res_folder.name

    if not knockout_dfs:
        print(f"no dfs for {protein}")
        res_df = pd.DataFrame()
    else:
        prop_df = prop_df[["nodes", protein]]   #.sort_values(protein, ascending=True)
        descending_prop_df = prop_df.sort_values(protein, ascending=False, ignore_index=True)
        descending_prop_df[protein] = pd.Series(np.arange(len(prop_df)))
        res_df = single_protein_sorted_df(descending_prop_df, knockout_dfs, protein, top_propagated_threshold)
    res_df._metadata.append(protein)   # setting res_df.name = protein doesn't work because of deserialization in return from subprocess. See: https://stackoverflow.com/questions/50372509/why-are-attributes-lost-after-copying-a-pandas-dataframe
    return res_df


def sort_per_protein(top_propagated_threshold,
                     prop_res_file: str, knockout_res_folders: list[str]) -> list[pd.DataFrame]:
    prop_df = pd.read_csv(prop_res_file)
    proteins = [c for c in prop_df.columns if c not in ["all", "nodes"] and not c.isnumeric()]
    folder_paths = [Path(p) for p in knockout_res_folders]
    folder_names = [p.name for p in folder_paths]
    unmatched_proteins = (set(proteins) | set(folder_names)) - (set(proteins) & set(folder_names))
    if unmatched_proteins:
        # raise ValueError(f"The following proteins are missing from folders or df: {unmatched_proteins}")
        print(f"no worries! {proteins} are here!")
    with Pool(min(len(proteins), 2*cpu_count())) as p:
        result_dfs = p.map(_worker, zip(repeat(top_propagated_threshold), repeat(prop_res_file), folder_paths))
    for result_df in result_dfs:
        result_df.name = result_df._metadata[0]
    return result_dfs


def resiliency_scores(top_propagated_threshold: int,
                      prop_res_file: str, knockout_res_folders: list[str], output_root: str) -> None:
    result_dfs = sort_per_protein(top_propagated_threshold, prop_res_file, knockout_res_folders)
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    for result_df in result_dfs:
        filename = result_df.name + ".csv"
        result_file_path = root / filename
        result_df.to_csv(result_file_path, index=False)


def test_resiliency_scores(top_propagated_threshold: int,
                           prop_res_file: str,  knockout_res_folders_root: str, output_root: str):
    # prop_res_file = r"/data/propagations/krogan_interactors/individual_interactors/all.csv"
    # output_root = r"/code/pipelines/infection_resiliency_pipeline/pipeline_objs/temp_resiliency"
    # knockout_res_folders = [r"/code/netprop_analysis/pipelines/infection_resiliency_pipeline/pipeline_objs/temp_outputs/e"]
    # prop_res_file = r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv"
    # output_root = r"C:\studies\thesis\code\netprop_analysis\pipelines\infection_resiliency_pipeline\pipeline_objs\temp_resiliency"
    # knockout_res_folders_root =  r"C:\studies\thesis\code\netprop_analysis\pipelines\infection_resiliency_pipeline\pipeline_objs\temp_outputs2"
    knockout_res_folders = [str(f) for f in Path(knockout_res_folders_root).glob("*")]
    resiliency_scores(top_propagated_threshold, prop_res_file, knockout_res_folders, output_root)


if __name__ == "__main__":
    test_resiliency_scores()