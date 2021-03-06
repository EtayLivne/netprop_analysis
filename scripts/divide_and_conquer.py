import json
from math import log
from pathlib import Path
from random import shuffle
from typing import Union, Callable
from crispr import get_huh_crispr
import networkx as nx
import pandas as pd
from random import sample

from netprop.models import PropagationResultModel, ConfigModel




def split_to_randomized_component(multi_prior_set_conf: Union[str, PropagationResultModel],
                                  output_dir: str, max_splits_per_component: int, override_conf_out_path: str=None) -> None:
    if type(multi_prior_set_conf) is str:
        multi_prior_set_conf = ConfigModel.parse_file(multi_prior_set_conf)

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    if override_conf_out_path:
        multi_prior_set_conf.output_dir_path = override_conf_out_path
    conf_out_dir_path = Path(Path(multi_prior_set_conf.output_dir_path))
    new_conf = multi_prior_set_conf.copy(deep=True)
    for prior_set in multi_prior_set_conf.prior_set:
        if len(prior_set.nodes) < 7:
            continue
        prior_set_id = prior_set.id
        new_conf.output_dir_path = conf_out_dir_path / prior_set_id
        new_conf = multi_prior_set_conf.copy(deep=True)
        new_conf.prior_set = []
        for i in range(max_splits_per_component):
            shuffled_prior_nodes = list(prior_set.nodes)
            shuffle(shuffled_prior_nodes)
            new_conf.prior_set.extend([
                {
                    "id": f"split_{i}_chunk_1",
                    "nodes": shuffled_prior_nodes[:len(shuffled_prior_nodes)//2],
                    "confidence": 0.7
                 },
                {
                    "id": f"split_{i}_chunk_2",
                    "nodes": shuffled_prior_nodes[len(shuffled_prior_nodes)//2:],
                    "confidence": 0.7}
            ])

        Path(new_conf.output_dir_path).mkdir(parents=True, exist_ok=True)
        filename = f"{prior_set_id}_split_conf.json"
        with open(str(output_dir/filename), 'w') as handler:
            json.dump(new_conf.dict(exclude_unset=True), handler, indent=4)
#
#
#
#
#

#
#
def merge_results_as_different_liquids(res_files: list[str], output_file: str) -> None:
    score_dfs = []
    merged_node_data = {}
    for f in res_files:
        file_prefix = Path(f).stem
        res = PropagationResultModel.parse_file(f)
        for node, node_data in res.nodes.items():
            if node not in merged_node_data:
                merged_node_data[node] = node_data
            merged = merged_node_data[node]
            if node_data.source_of:
                merged.source_of.union(node_data.source_of)
            for l in list(node_data.liquids.keys()):
                merged.liquids[file_prefix + f"_{l}"] = node_data.liquids[l]

    for node_data in merged_node_data.values():
        node_data.source_of = list(node_data.source_of)
    output = PropagationResultModel.parse_file(res_files[0])
    output.nodes = merged_node_data

    with open(output_file, 'w') as handler:
        json.dump(output.dict(), handler, indent=4)


def _intersect_two_chunks(chunk_1: PropagationResultModel, chunk_2: PropagationResultModel, k: int) -> set[str]:
    chunk_1_top_propagated = sorted([n for n in chunk_1.nodes.keys()],
                                    key=lambda n: chunk_1.nodes[n].liquids["info"], reverse=True)[:k]
    chunk_2_top_propagated = sorted([n for n in chunk_2.nodes.keys()],
                                    key=lambda n: chunk_2.nodes[n].liquids["info"], reverse=True)[:k]

    print(f"found {len(set(chunk_1_top_propagated) & set(chunk_2_top_propagated))} intersections")
    return set(chunk_1_top_propagated) & set(chunk_2_top_propagated)


def _count_intersection_appearances_per_node(intersection_sets: list[set[str]]) -> dict[str, int]:
    all = set().union(*intersection_sets)
    return {node: len([s for s in intersection_sets if node in s]) for node in all}


def get_nodes_in_intersection(results: list[str], k: int=50, output_file: str=None):
    sorted_results = sorted(results)    #alphabatical sorting would place the files as minor-major all the way
    even_chunks = sorted_results[::2]
    odd_chunks = sorted_results[1::2]
    # names = [res.split("\\")[-1].split("_chunk_1")[0] for res in even_chunks]
    intersection_sets = []
    for i in range(len(even_chunks)):
        intersection_sets.append(
            _intersect_two_chunks(PropagationResultModel.parse_file(even_chunks[i]),
                                  PropagationResultModel.parse_file(odd_chunks[i]), k)
        )

    return _count_intersection_appearances_per_node(intersection_sets)


def intersect_crispr_ko_with_node_intersections(path_to_crispr_ko: str, intersection_appearance_dicts: list[dict[str, int]]) -> float:
    crispr_kos = set(get_huh_crispr(path_to_crispr_ko))
    # crispr_kos = set(crispr_kos.to_dict().keys())
    all_intersection = set().union(*[set(d.keys()) for d in intersection_appearance_dicts])

    interesting = crispr_kos & all_intersection

    # for node in interesting:
    #     scores = []
    #     for i, d in enumerate(intersection_appearance_dicts):
    #         if node in d:
    #             scores.append((i, d[node]))
    #     print(f"{node}: {scores}")
    import pandas as pd
    res = PropagationResultModel.parse_file(r"D:\data\propagations\krogan_interactors\merged\all.json").prop_scores_as_series()
    res = res.sort_values("all.info", ascending=False)
    res_node_index = pd.Index(res["nodes"])
    for node in interesting:
        print(f"[{node}: {res_node_index.get_loc(node)}")

    return len(crispr_kos & all_intersection) / len(crispr_kos)















"""

def cov_interactors_prior_set(nw: nx.Graph, base_conf_path: str, new_conf_path: str,
                              prop_output_dir: str=None, confidence=0.7, method="iterative",
                              merge_sources: bool=False):
    conf_to_patch = ConfigModel.parse_file(base_conf_path)
    if merge_sources:
        conf_to_patch.prior_set = [
            {
                "nodes": [{"id": n, "source_of": ["info"]} for n in get_covid_interactors(nw)],
                "id": "cov_interactors",
                "confidence": confidence
            }
        ]

    else:
        conf_to_patch.prior_set = [
            {
                "nodes": [{"id": n, "source_of": ["info"]} for n in interactors],
                "id": f"{cov_protein}_interactors",
                "confidence": confidence
            }
            for cov_protein, interactors in gen_covid_interactors_by_protein(nw)
        ]

    if prop_output_dir:
        conf_to_patch.output_dir_path = prop_output_dir
    conf_to_patch.method = method
    with open(new_conf_path, 'w') as handler:
        json.dump(conf_to_patch.dict(exclude_unset=True), handler, indent=4)



def partial_cov_interactors_prior_set(confidence=0.7):


    cov_interactors = get_covid_interactors(nw)

    boring_old_conf = ConfigModel.parse_file(r"D:\configurations\stefan1\stefan_conf_backup.json")
    all_nodes = [{"id": n, "source_of": ["info"]} for n in cov_interactors]
    shuffle(all_nodes)
    chunks = [list(all_nodes[ceil(len(all_nodes) * i / 5):ceil(len(all_nodes) * (i + 1) / 5)]) for i in range(5)]
    chunks = [[n["id"] for n in chunk] for chunk in chunks]
    boring_old_conf.prior_set = []
    for i, chunk in enumerate(chunks):
        boring_old_conf.prior_set.append({
            "nodes": [n for n in all_nodes if n["id"] not in chunk],
            "id": f"chunk_{i}_majority",
            "confidence": confidence
        })
        boring_old_conf.prior_set.append({
            "nodes": [n for n in all_nodes if n["id"] in chunk],
            "id": f"chunk_{i}_minority",
            "confidence": 0.7
        })

    boring_old_conf.output_dir_path = f"D:\\data\\propagations\\stefan_cross_validation"
    with open(f"D:\\configurations\\stefan1\\stefan_conf.json", 'w') as handler:
        json.dump(boring_old_conf.dict(exclude_unset=True), handler, indent=4)


def specific_protein_partial_cov_interactors_prior_set(confidence=0.7):
    nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
                                         r"D:\data\networks\krogan_cov_to_human.cx",
                                         r"D:\data\other\symbol_to_entrezgene.json",
                                         merge_covid=False).load()

    for covid_protein, interactors in gen_covid_interactors_by_protein(nw):
        boring_old_conf = ConfigModel.parse_file(r"D:\configurations\stefan1\stefan_conf_backup.json")
        all_nodes = [{"id": n, "source_of": ["info"]} for n in interactors]
        shuffle(all_nodes)
        chunks = [list(all_nodes[ceil(len(all_nodes) * i / 5):ceil(len(all_nodes) * (i + 1) / 5)]) for i in range(5)]
        chunks = [[n["id"] for n in chunk] for chunk in chunks]
        boring_old_conf.prior_set = []
        for i, chunk in enumerate(chunks):
            boring_old_conf.prior_set.append({
                "nodes": [n for n in all_nodes if n["id"] not in chunk],
                "id": f"chunk_{i}_majority",
                "confidence": confidence
            })
            boring_old_conf.prior_set.append({
                "nodes": [n for n in all_nodes if n["id"] in chunk],
                "id": f"chunk_{i}_minority",
                "confidence": 0.7
            })

        boring_old_conf.output_dir_path = f"D:\\data\\propagations\\stefan_cross_validation\\{covid_protein}"
        Path.mkdir(Path(boring_old_conf.output_dir_path), exist_ok=True)
        with open(f"D:\\configurations\\stefan1\\stefan_conf_{covid_protein}.json", 'w') as handler:
            json.dump(boring_old_conf.dict(exclude_unset=True), handler, indent=4)


def random_subgroups_partial_cov_interactors_prior_set(group_size: float, num_groups: int, confidence=0.7, ):
    nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
                                         r"D:\data\networks\krogan_cov_to_human.cx",
                                         r"D:\data\other\symbol_to_entrezgene.json",
                                         merge_covid=False).load()

    cov_interactors = get_covid_interactors(nw)
    boring_old_conf = ConfigModel.parse_file(r"D:\configurations\stefan1\stefan_conf_backup.json")
    all_nodes = [{"id": n, "source_of": ["info"]} for n in cov_interactors]

    subgroups = [sample(all_nodes, group_size) for _ in range(num_groups)]
    Path.mkdir(Path(f"D:\\configurations\\stefan1\\random_subgroups\\{group_size}"), exist_ok=True, parents=True)
    for subgroups_num, subgroup in enumerate(subgroups):
        chunks = [list(subgroup[ceil(len(subgroup) * i / 5):ceil(len(subgroup) * (i + 1) / 5)]) for i in range(5)]
        chunks = [[n["id"] for n in chunk] for chunk in chunks]
        boring_old_conf.prior_set = []
        for i, chunk in enumerate(chunks):
            boring_old_conf.prior_set.append({
                "nodes": [n for n in subgroup if n["id"] not in chunk],
                "id": f"chunk_{i}_majority",
                "confidence": confidence
            })
            boring_old_conf.prior_set.append({
                "nodes": [n for n in subgroup if n["id"] in chunk],
                "id": f"chunk_{i}_minority",
                "confidence": 0.7
            })

        boring_old_conf.output_dir_path = f"D:\\data\\propagations\\stefan_cross_validation\\random_subgroups\\{group_size}\\{subgroups_num}"
        Path.mkdir(Path(boring_old_conf.output_dir_path), exist_ok=True, parents=True)
        with open(f"D:\\configurations\\stefan1\\random_subgroups\\{group_size}\\stefan_rand_conf_{subgroups_num}.json",'w') as handler:
            json.dump(boring_old_conf.dict(exclude_unset=True), handler, indent=4)



# construct a prior set where each covid interactor propagates a unique liquid
def individual_cov_interactors_prior_set(confidence=0.7):
    nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
                                         r"D:\data\networks\krogan_cov_to_human.cx",
                                         r"D:\data\other\symbol_to_entrezgene.json").load()

    cov_interactors = get_covid_interactors(nw)

    boring_old_conf = ConfigModel.parse_file(r"D:\configurations\stefan1\stefan_conf_backup.json")
    boring_old_conf.prior_set = [
        {
            "nodes": [{"id": n, "source_of": [n]} for n in cov_interactors],
            "id": "covid_interactors",
            "confidence": confidence
        }
    ]
    boring_old_conf.output_dir_path = f"D:\\data\\propagations\\stefan_unique_liquids"
    with open(f"D:\\configurations\\stefan1\\stefan_conf.json", 'w') as handler:
        json.dump(boring_old_conf.dict(exclude_unset=True), handler, indent=4)


"""


#
#
#
# def merge_prior_sets(conf: str, new_conf: str, confidence: float=0.7) -> None:
#     conf_to_patch = ConfigModel.parse_file(conf)
#     all_nodes = []
#     for p in conf_to_patch.prior_set:
#         all_nodes += p.nodes
#
#     conf_to_patch.prior_set = [
#         {
#             "nodes": all_nodes,
#             "id": "all",
#             "confidence": 0.7
#         }
#     ]
#
#     with open(new_conf, 'w') as handler:
#         json.dump(conf_to_patch.dict(exclude_unset=True), handler, indent=4)