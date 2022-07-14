from viral_analysis_pipeline.pipeline_objs.pipeline_objs_consts import *
from pathlib import Path
from utils.utils import load_json, dump_json
import pandas as pd
from viral_analysis_pipeline.analysis_helpers import _split_a_list, _calc_num_splits, _analyze_intersection_res_dir
from random import shuffle
from pipelines import NonRepeatingPipeline

def _prep_env(attrs: dict):
    # metadata_file = attrs[METADATA_FILE_PATH_ATTR]
    # res_file = attrs[REAL_PROP_RES_PATH_ATTR]
    virus_res_root = attrs[VIRUS_RES_ROOT_ATTR]

    output_root = Path(virus_res_root) / "intersection_data"

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

    attrs[INTER_OUTPUT_DIR_ATTR] = inter_output_dir
    attrs[CROSS_OUTPUT_DIR_ATTR] = cross_output_dir
    attrs[BY_SIZE_OUTPUT_DIR_ATTR] = by_size_output_dir
    attrs[ALL_OUTPUT_DIR_ATTR] = all_interactors_output_dir
    attrs[INTERSECTION_OUTPUT_DIR_ATTR] = intersections_output_dir


def _randomly_split(attrs: dict) -> dict[list[list[str]]]:
    metadata_path = attrs[METADATA_FILE_PATH_ATTR]
    prot_map = load_json(metadata_path)
    proteins = attrs.get("proteins", None)

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


def _random_cross_proteins_splits(attrs: dict):
    metadata_path = attrs[METADATA_FILE_PATH_ATTR]
    res_path = attrs[REAL_PROP_RES_PATH_ATTR]
    proteins = attrs.get("proteins", None)

    res = pd.read_csv(res_path)
    prot_to_interactors = load_json(metadata_path)
    if proteins is None:
        proteins = [c for c in res.columns if c != "nodes" and not c.isnumeric()]

    actual_proteins = []
    for i in range(len(proteins) - 1):
        p1 = proteins[i]
        for j in range(i + 1, len(proteins)):
            p2 = proteins[j]
            key = f"{p1},{p2}"
            prot_to_interactors[key] = list(set(prot_to_interactors[p1]) | set(prot_to_interactors[p2]))
            actual_proteins.append(key)
    return {
        protein: [_split_a_list(prot_to_interactors[protein])
                  for i in range(_calc_num_splits(len(prot_to_interactors[protein])))] for protein in actual_proteins if
        len(prot_to_interactors[protein]) >= 10
    }


def _split_by_size(attrs: dict):
    res_path = attrs[REAL_PROP_RES_PATH_ATTR]
    res = pd.read_csv(res_path)
    interactors_columns = [c for c in res.columns if c.isnumeric()]
    splits = dict()
    for split_size in [10, 15, 20, 30, 40]:
        shuffle(interactors_columns)
        interactors_in_split = interactors_columns[:split_size]
        splits[split_size] = [_split_a_list(interactors_in_split) for _ in range(_calc_num_splits(split_size))]

    return splits

def create_inter_splits_data(attrs: dict):
    inter_output_dir = attrs[INTER_OUTPUT_DIR_ATTR]
    print(f"   ***********   WORKING ON SPLIT TYPE INTER  ***********")
    for i in range(10):
        splits = _randomly_split(attrs)
        dump_json(splits, str(inter_output_dir / f"splits_{i}.json"))
    _analyze_intersection_res_dir(inter_output_dir, attrs[INTERSECTION_OUTPUT_DIR_ATTR], "inter", attrs[REAL_PROP_RES_PATH_ATTR])

def create_cross_splits_data(attrs: dict):
    cross_output_dir = attrs[CROSS_OUTPUT_DIR_ATTR]
    print(f"   ***********   WORKING ON SPLIT TYPE CROSS  ***********")
    for i in range(10):
        splits = _random_cross_proteins_splits(attrs)
        dump_json(splits, str(cross_output_dir / f"splits_{i}.json"))
    _analyze_intersection_res_dir(cross_output_dir, attrs[INTERSECTION_OUTPUT_DIR_ATTR], "cross",
                                  attrs[REAL_PROP_RES_PATH_ATTR])

def create_by_size_splits_data(attrs: dict):
    by_size_output_dir = attrs[BY_SIZE_OUTPUT_DIR_ATTR]
    print(f"   ***********   WORKING ON SPLIT TYPE CROSS  ***********")
    for i in range(10):
        splits = _split_by_size(attrs)
        dump_json(splits, str(by_size_output_dir / f"splits_{i}.json"))

    _analyze_intersection_res_dir(by_size_output_dir, attrs[INTERSECTION_OUTPUT_DIR_ATTR], "by_size",
                                  attrs[REAL_PROP_RES_PATH_ATTR])
# def create_all_splits_data(attrs:dict):
#     by_size_output_dir = attrs[BY_SIZE_OUTPUT_DIR_ATTR]
#     print(f"   ***********   WORKING ON SPLIT TYPE CROSS  ***********")
#     for i in range(10):
#         splits = _split_by_size(attrs)
#         dump_json(splits, str(by_size_output_dir / f"splits_{i}.json"))


def get_pipeline(virus_name: str, reset_state: bool=False) -> NonRepeatingPipeline:
    pipe = NonRepeatingPipeline(state_suffix=virus_name + "intersection_analysis_state", reset_state=reset_state)
    pipe.add_steps(
        [
            _prep_env,
            create_inter_splits_data,
            create_by_size_splits_data,
            create_cross_splits_data
        ]
    )
    return pipe