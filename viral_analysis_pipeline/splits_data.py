from pathlib import Path
from utils.utils import load_json, dump_json
from functools import partial
import pandas as pd
from multiprocessing import Queue, Process
from math import log
from random import shuffle
from typing import Union


# calculated so that there's more than a 95% chance for every pair of proteins to be in opposiong halves at least once.
def _calc_num_splits(set_size: int) -> int:
    log_base = (set_size - 2)/(set_size - 1)
    try:
        return int(log(0.1, log_base)) + 1
    except:
        raise ValueError(f"log_base: {log_base}")


def _split_a_list(to_split: list) -> tuple[list, list]:
    shuffled = list(to_split)
    shuffle(shuffled)
    return shuffled[:len(shuffled)//2], shuffled[len(shuffled)//2:]


def _randomly_split(prot_to_interactors: dict, proteins: Union[str, list[str]]=None) -> dict[list[list[str]]]:
    prot_map = load_json(prot_to_interactors)

    if proteins is None:
        proteins = list(prot_map.keys())
    elif not isinstance(proteins, list):
        proteins = [proteins]

    prot_map["all"] = list(set().union(*[v for v in prot_map.values()]))
    return {
        protein: [_split_a_list(prot_map[protein]) for
                  i in range(_calc_num_splits(len(prot_map[protein])))]
        for protein in proteins if len(prot_map[protein]) >= 10
    }

def create_splits_files(viral_prot_to_interactor: str, splits_file_path: str):
    prot_to_interactor = load_json(viral_prot_to_interactor)
    splits_dict = _randomly_split(prot_to_interactor)
    dump_json(splits_dict, splits_file_path)



def _calc_single_intersection(df: pd.DataFrame, col_set1: list[str], col_set2: list[str], threshold: int) -> list[str]:
    # print(df.columns)
    # print(col_set1)
    # raise Exception
    df1 = df[col_set1]
    df2 = df[col_set2]
    combined_set1 = df1.apply(lambda row: row.sum(), axis=1).sort_values(ascending=False)
    combined_set2 = df2.apply(lambda row: row.sum(), axis=1).sort_values(ascending=False)
    discovered = list(set(combined_set1.iloc[:threshold].index) & set(combined_set2.iloc[:threshold].index))
    return discovered


def _get_protein_intersections_dict_worker(df: pd.DataFrame, protein: str, protein_splits: list[list[list[str]]], threshold: int, queue: Queue) -> dict[str, int]:
    print(f"started on {protein}")
    intersection_sets = []
    counter = 0
    for split in protein_splits:
        intersection_sets.append(_calc_single_intersection(df, split[0], split[1], threshold))
        counter += 1
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
    proteins_to_check = list(splits.keys())
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


def split_by_size(res_path: str, splits_sizes: list[int]):
    res = pd.read_csv(res_path)
    interactors_columns = [c for c in res.columns if c.isnumeric()]
    splits = dict()
    for split_size in splits_sizes:
        shuffle(interactors_columns)
        interactors_in_split = interactors_columns[:split_size]
        splits[split_size] = [_split_a_list(interactors_in_split) for i in range(_calc_num_splits(split_size))]

    return splits

def generate_intersection_data(metadata_file: str, res_file: str, output_root: str, split_repetitions: int, intersection_thresholds: list[int]):
    root_path = Path(output_root)
    inter_output_dir = root_path / "inter"
    cross_output_dir = root_path / "cross"
    by_size_output_dir = root_path / "by_size"
    all_interactors_output_dir = root_path / "all_interactors"
    intersections_output_dir = root_path / "intersection_results"

    root_path.mkdir(exist_ok=True)
    inter_output_dir.mkdir(exist_ok=True)
    cross_output_dir.mkdir(exist_ok=True)
    by_size_output_dir.mkdir(exist_ok=True)
    all_interactors_output_dir.mkdir(exist_ok=True)
    intersections_output_dir.mkdir(exist_ok=True)

    split_types = [
        (inter_output_dir, partial(_randomly_split, metadata_file, None), "inter"),
        (cross_output_dir, partial(random_cross_proteins_splits, res_file, metadata_file), "cross"),
        (by_size_output_dir, partial(split_by_size, res_file, [10, 15, 20, 30, 40]), "by_size"),
        (all_interactors_output_dir, partial(_randomly_split, metadata_file, ["all"]), "all")
    ]

    for split_type in split_types:
        output_dir = split_type[0]
        splitting_func = split_type[1]
        type_name = split_type[2]

        print(f"   ***********   WORKING ON SPLIT TYPE {type_name}  ***********")
        for i in range(split_repetitions):
            splits = splitting_func()
            dump_json(splits, str(output_dir / f"splits_{i}.json"))

        split_files = list(output_dir.glob("*"))
        for i in range(len(split_files)):
            print(f"   ***********   WORKING ON INTERSECTIONS {i}   ***********")
            specific_file_output_dir = intersections_output_dir / type_name / f"intersections_from_file_{i}"
            specific_file_output_dir.mkdir(exist_ok=True, parents=True)
            create_intersections_data(res_file,
                                      str(split_files[i]),
                                      str(specific_file_output_dir / f"intersection_res"),
                                      intersection_thresholds)