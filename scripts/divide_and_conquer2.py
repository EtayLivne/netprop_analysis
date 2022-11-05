from typing import Union
import pandas as pd
import numpy as np
import json
from random import shuffle
from math import log
from utils.queue_managers import load_json, dump_json
from multiprocessing import Process, Queue
from metrics.roc_auc import SinglePropROC
from itertools import product
from crispr import get_huh_crispr


class PropResNewFormat:
    def __init__(self, path_to_csv: str, path_to_metadata: str):
        self.prop_res = pd.read_csv(path_to_csv)
        with open(path_to_metadata, 'r') as handler:
            self.protein_to_interactors = json.load(handler)




def random_cross_proteins_splits(res_path: str, metadata_path: str, proteins: list[str]=None):
    res = pd.read_csv(res_path)
    prot_to_interators = load_json(metadata_path)
    if proteins is None:
        proteins = [c for c in res.columns if c != "nodes" and not c.isnumeric()]

    actual_proteins = []
    for i in range(len(proteins) - 1):
        p1 = proteins[i]
        for j in range(i + 1, len(proteins)):
            p2 = proteins[j]
            key = f"{p1},{p2}"
            prot_to_interators[key] = list(set(prot_to_interators[p1]) | set(prot_to_interators[p2]))
            actual_proteins.append(key)
    return {
        protein: [_split_a_list(prot_to_interators[protein])
                  for i in range(_calc_num_splits(len(prot_to_interators[protein])))] for protein in actual_proteins if
        len(prot_to_interators[protein]) >= 10
    }


def _calc_single_intersection(df: pd.DataFrame, col_set1: list[str], col_set2: list[str], threshold: int) -> list[str]:
    df1 = df[col_set1]
    df2 = df[col_set2]
    combined_set1 = df1.apply(lambda row: row.sum(), axis=1).sort_values(ascending=False)
    combined_set2 = df2.apply(lambda row: row.sum(), axis=1).sort_values(ascending=False)
    discovered = list(set(combined_set1.iloc[:threshold].index) & set(combined_set2.iloc[:threshold].index))
    return discovered


def _get_protein_intersections_dict(df: pd.DataFrame, protein_splits: list[list[list[str]]], threshold: int) -> dict[str, int]:
    intersection_sets = []
    for split in protein_splits:
        intersection_sets.append(_calc_single_intersection(df, split[0], split[1], threshold))

    all_sets = set().union(*intersection_sets)
    return {node: len([s for s in intersection_sets if node in s])/len(protein_splits) for node in all_sets}

def _get_protein_intersections_dict_worker(df: pd.DataFrame, protein: str, protein_splits: list[list[list[str]]], threshold: int, queue: Queue) -> dict[str, int]:
    print(f"started on {protein}")
    intersection_sets = []
    counter = 0
    for split in protein_splits:
        intersection_sets.append(_calc_single_intersection(df, split[0], split[1], threshold))
        counter +=1
        if counter % 10 == 0:
            print(f"did {counter} for {protein} with threshold {threshold}")
    all_sets = set().union(*intersection_sets)
    queue.put({
        protein: {node: len([s for s in intersection_sets if node in s])/len(protein_splits) for node in all_sets}
    })

    print(f"done with {protein} with threshold {threshold}")


def create_intersections_data(res_file: str, splits_file: str, output_prefix: str, thresholds: list[int]) -> None:
    df = pd.read_csv(res_file)
    splits = load_json(splits_file)
    q = Queue()
    proteins_to_check = list(splits.keys())[:15] # TEMP DELETE
    for t in thresholds:
    # input_tuples = product([list(splits.keys()), thresholds])
        workers = [Process(target=_get_protein_intersections_dict_worker, args=(df, p, splits[p], t, q)) for p in proteins_to_check]
        for worker in workers:
            worker.start()
        for worker in workers:
            worker.join()
        overall_res = dict()
        for i in range(len(workers)):
            overall_res.update(q.get())
        # split_res = {protein: _get_protein_intersections_dict(df, splits[protein], threshold) for protein in splits}
        output_file = output_prefix + f"_{t}.json"
        print(f"dumping to {output_file}")
        dump_json(overall_res, output_file)



def kinda_correlate_interesection_to_prop(res_file: str, splits_res_file: str, top_res_threshold: int) -> dict[str, float]:
    df = pd.read_csv(res_file)
    splits_res = load_json(splits_res_file)
    top_prop_nodes = {protein:  list(df[protein].sort_values(ascending=False).index)[:top_res_threshold] for protein in splits_res}
    top_intersection_nodes = {protein: [int(i) for i in intersections if intersections[i] > 0.5]
                              for protein, intersections in splits_res.items()}

    kinda_correlation = {protein: len([n for n in top_intersected if n in top_prop_nodes[protein]])
                         for protein, top_intersected in top_intersection_nodes.items()}

    return kinda_correlation


# def intersect_top_propagated(res_file: str, splits_res_file: str, top_res_threshold: int) -> dict[str, float]:
#     df = pd.read_csv(res_file)
#     proteins = [c for c in df.columns if c.isnumeric()]
#     df = df[proteins]
#     top_prop_nodes = pd.DataFrame({column:  list(df[column].sort_values().index)[:top_res_threshold] for column in df.columns})
#     intersection_df = pd.DataFrame({protein: [list(set(df[protein]) & set(other_protein)) for other_protein in proteins] for protein in proteins}, index=proteins)
#     # np.fill_diagonal(intersection_df.values, None) do away with intersections with themselves   
#
#     top_intersection_nodes = {protein: [int(i) for i in intersections if intersections[i] > 0.5]
#                               for protein, intersections in splits_res.items()}
#
#     kinda_correlation = {protein: len([n for n in top_intersected if n in top_prop_nodes[protein]])
#                          for protein, top_intersected in top_intersection_nodes.items()}
#
#     return kinda_correlation


def intersect_top_propagated_with_cripsr_kos(res_file: str, crispr_path: str, top_prop_threshold: int) -> dict[str, list[str]]:
    df = pd.read_csv(res_file, index_col="nodes")
    proteins = [c for c in df.columns if not c.isnumeric()]
    df = df[proteins]
    for i in range(1, 3):
        df[f"L{i}"] = df[proteins].apply(lambda row: sum(row**i)**(1./i), axis=1)
    df["L_inf"] = df[proteins].apply(lambda row: row.max(), axis=1)
    df = df.reset_index()
    df["nodes"] = df["nodes"].apply(str)
    # top_prop_nodes = pd.DataFrame({column:  list(df[column].sort_values(ascending=False).index)[:top_prop_threshold] for column in df.columns})
    # top_prop_nodes = top_prop_nodes.reset_index()
    # top_prop_nodes.to_csv(r"D:\data\propagations\krogan_interactors\individual_interactors\proteins.csv")
    df.to_csv(r"D:\data\propagations\krogan_interactors\individual_interactors\proteins.csv", index=False)
    metric = SinglePropROC(r"D:\data\propagations\krogan_interactors\individual_interactors\proteins.csv",
                           hit_set_load_method=get_huh_crispr)
    metric.load_hit_set(crispr_path)

    metric.show_results(metric.calc_metric())



    # crispr_kos = [float(ko) for ko in get_huh_crispr(crispr_path)]
    #
    # return {protein: [target for target in top_prop_nodes[protein] if target in crispr_kos] for protein in top_prop_nodes.columns}

def _get_top_intersections(intersections_file: str, threshold: float=0.5) -> dict[str, list[str]]:
    intersections = load_json(intersections_file)
    return {
        protein: [node for node, percent_intersections in intersections[protein].items() if percent_intersections > threshold]
        for protein in intersections
    }


def crispr_targets_in_top_intersections(intersections_file: str, crispr_file: str) -> dict[str, list[str]]:
    intersections = _get_top_intersections(intersections_file)
    crispr = get_huh_crispr(crispr_file)
    return {protein: [n for n in protein if n in crispr] for protein in intersections}


def split_by_size(res_path: str, splits_sizes: list[int]):
    res = pd.read_csv(res_path)
    interactors_columns = [c for c in res.columns if c.isnumeric()]
    splits = dict()
    for split_size in splits_sizes:
        shuffle(interactors_columns)
        interactors_in_split = interactors_columns[:split_size]
        splits[split_size] = [_split_a_list(interactors_in_split) for i in range(_calc_num_splits(split_size))]

    return splits