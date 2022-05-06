from utils.new_file_loaders import NetpropCovToHumanLoader
from score_from_netprop import *
from crispr import get_crispr_rankings
from scipy.stats import spearmanr
from pathlib import Path
from metrics.roc_auc import HuhCrisprROC
from netprop.networks.loaders import NetpropNetwork, CombinedHumanCovidNetworkLoader
from netprop.models import ConfigModel
from netprop.propagation import propagate_from_config
import networkx as nx
from random import shuffle, sample
from math import ceil
import json
import numpy as np
from scripts.draw_network import draw_prop_network, draw_network
from multiprocessing import Process
from scripts.file import calc_protein_p_value, get_intersection_quality_metrics, intersection_bar_plot


def merge_prior_sets(conf: str, new_conf: str, confidence: float=0.7) -> None:
    conf_to_patch = ConfigModel.parse_file(conf)
    all_nodes = []
    for p in conf_to_patch.prior_set:
        all_nodes += p.nodes

    conf_to_patch.prior_set = [
        {
            "nodes": all_nodes,
            "id": "all",
            "confidence": 0.7
        }
    ]

    with open(new_conf, 'w') as handler:
        json.dump(conf_to_patch.dict(exclude_unset=True), handler, indent=4)

# Returns all nodes in the graph that are not covid nodes themselves but have a covid neighbor
def get_covid_interactors(network: nx.Graph) -> set[str]:
    cov_interactors = set()
    for _, interactors in gen_covid_interactors_by_protein(network):
        cov_interactors |= interactors
    return cov_interactors


def gen_covid_interactors_by_protein(network: nx.Graph) -> set[str]:
    cov_nodes = [n for n, data in network.nodes(data=True) if data["species_id"] == "sars-cov-2"]
    for node in cov_nodes:
        yield node, set(network.neighbors(node))


def cov_interactors_prior_set(nw: nx.Graph, base_conf_path: str, new_conf_path: str,
                              prop_output_dir: str=None, confidence=0.7, method="iterative",
                              merge_sources: bool=False):
    # nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
    #                                      r"D:\data\networks\krogan_cov_to_human.cx",
    #                                      r"D:\data\other\symbol_to_entrezgene.json").load()

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
    nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
                                         r"D:\data\networks\krogan_cov_to_human.cx",
                                         r"D:\data\other\symbol_to_entrezgene.json").load()

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

def propagate_stefan():
    propagate_from_config(r"D:\configurations\stefan1\stefan_conf.json", ordering={"prior_set": 100})


def stefan_propagations(conf_files: list[str]):
    for file in conf_files:
        print(f"propagating {file}")
        propagate_from_config(str(file), ordering={"prior_set": 100}, max_processes=7   )


def intersect_top_propagated(results: list[str], k: int=50, output_file: str=None):
    sorted_results = sorted(results)    #alphabatical sorting would place the files as minor-major all the way
    majorities = sorted_results[::2]
    minorities = sorted_results[1::2]
    names = [res.split("\\")[-1].split("_majority")[0] for res in majorities]
    res_dict = dict()
    for i in range(len(majorities)):
        major_res = PropagationResultModel.parse_file(majorities[i])
        major_res_top_propagated = sorted([n for n in major_res.nodes.keys()], key=lambda n: major_res.nodes[n].liquids["info"], reverse=True)[:k]
        minor_res = PropagationResultModel.parse_file(minorities[i])
        minor_res_top_propagated = sorted([n for n in minor_res.nodes.keys()], key=lambda n: minor_res.nodes[n].liquids["info"], reverse=True)[:k]
        intersection = set(major_res_top_propagated) & set(minor_res_top_propagated)
        if len(intersection)/k > 0:
            print(f"{names[i]}: {len(intersection)/k}")
        res_dict[names[i]] = len(intersection)/k

    if output_file:
        with open(output_file, 'w') as handler:
            json.dump(res_dict, handler, indent=4)


def average_intersection(result_folders: list[str], k: int=50):
    scores = []
    for result_folder in result_folders:
        files = sorted([str(f) for f in Path(result_folder).glob("*.json")]) #alphabatical sorting would place the files as minor-major all the way
        majorities = files[::2]
        minorities = files[1::2]
        names = [res.split("\\")[-1].split("_majority")[0] for res in majorities]
        for i in range(len(majorities)):
            major_res = PropagationResultModel.parse_file(majorities[i])
            major_res_top_propagated = sorted([n for n in major_res.nodes.keys()], key=lambda n: major_res.nodes[n].liquids["info"], reverse=True)[:k]
            minor_res = PropagationResultModel.parse_file(minorities[i])
            minor_res_top_propagated = sorted([n for n in minor_res.nodes.keys()], key=lambda n: minor_res.nodes[n].liquids["info"], reverse=True)[:k]
            intersection = set(major_res_top_propagated) & set(minor_res_top_propagated)
            scores.append(len(intersection)/k)

    return scores


def top_prop_by_source(res_file: str, k: int=50):
    nodes = PropagationResultModel.parse_file(res_file).nodes
    liquids = nodes['1'].liquids.keys()
    return {
        liquid: sorted([n for n in nodes.keys()], key=lambda n: nodes[n].liquids[liquid], reverse=True)[:k]
        for liquid in liquids
    }


def most_intersected(dict_by_liquid: dict):
    nodes = {}
    for liquid, top_nodes in dict_by_liquid.items():
        for node in top_nodes:
            if node not in nodes:
                nodes[node] = []
            nodes[node].append(liquid)
    return nodes


def compare_results(res_1: str, res_2: str) -> bool:
    res_1_nodes = PropagationResultModel.parse_file(res_1).nodes
    res_2_nodes = PropagationResultModel.parse_file(res_2).nodes
    from math import fabs
    diffs = []
    counters = {"normal": 0, "weird": 0}
    for node, data in res_1_nodes.items():
        score_1, score_2 = data.liquids["info"], res_2_nodes[node].liquids["info"]

        if min(score_1, score_2) == 0:
            deviation = max(score_1, score_2)
            counters["weird"] += 1
        else:
            deviation = max(score_1, score_2) / min(score_1, score_2)
            counters["normal"] += 1

        #diffs.append(fabs(score_1- score_2)/max(score_1, score_2))
        diffs.append(deviation)
    print(f"max: {max(diffs)}\navg: {sum(diffs)/len(diffs)}")
    print(len([d for d in diffs if d < 1]))


# def validate_statistical_significance()


def more_detailed_intersection_data(results: list[str], k: int=50, output_file: str=None):
    sorted_results = sorted(results)    #alphabatical sorting would place the files as minor-major all the way
    majorities = sorted_results[::2]
    minorities = sorted_results[1::2]
    names = [res.split("\\")[-1].split("_majority")[0] for res in majorities]
    res_dict = {"chunks": dict(), "proteins": dict()}
    all_intersections = set()
    for i in range(len(majorities)):
        major_res = PropagationResultModel.parse_file(majorities[i])
        major_res_top_propagated = sorted([n for n in major_res.nodes.keys()], key=lambda n: major_res.nodes[n].liquids["info"], reverse=True)[:k]
        minor_res = PropagationResultModel.parse_file(minorities[i])
        minor_res_top_propagated = sorted([n for n in minor_res.nodes.keys()], key=lambda n: minor_res.nodes[n].liquids["info"], reverse=True)[:k]
        intersection = set(major_res_top_propagated) & set(minor_res_top_propagated)
        if len(intersection)/k > 0:
            print(f"{names[i]}: {len(intersection)/k}")
        res_dict["chunks"][names[i]] = {"size": len(intersection)/k, "intersecting_proteins": list(intersection)}
        all_intersections = all_intersections.union(intersection)

    res_dict["proteins"] = {
        p: len([1 for chunk in res_dict["chunks"].values() if p in chunk["intersecting_proteins"]]) for p in all_intersections
    }

    if output_file:
        with open(output_file, 'w') as handler:
            json.dump(res_dict, handler, indent=4)



if __name__ == "__main__":
    # g = nx.cycle_graph(24)
    # np.random.seed(100)
    # g = nx.connected_watts_strogatz_graph(20, 5, 0.5, seed=100)
    # labels = {}
    # scores = np.random.randint(0, 255, 20)
    # node_scores = {f"{i}": scores[i] for i in range(20)}
    # draw_network(g, labels, node_scores)
    # df = pd.DataFrame({"node": ["CHD6", "MBPTS2", "ACE2", "VAC14", "EXOC2"], "score": [0.86, 0.85, 0.63, 0.57, 0.52]})
    # print(df)

    huh_7_nw = NetpropNetwork(r"D:\data\networks\metastudy\complete_networks\huh7_only_fixed.json").load()
    base_conf = r"C:\studies\thesis\code\cicd\data_cache\metastudy_conf.json"
    new_conf = r"C:\studies\thesis\code\cicd\data_cache\metastudy_all_interactors_conf.json"
    prop_output_dir = "/output/props/all_interactors"
    cov_interactors_prior_set(huh_7_nw, base_conf, new_conf, prop_output_dir=prop_output_dir)

    # df = pd.read_csv(r"D:\data\networks\metastudy\protein_interactome_translated.csv")
    # df["Bait"] = df["Bait"].apply(lambda name: name.upper())
    # df.to_csv(r"D:\data\networks\metastudy\protein_interactome_translated_upper_case.csv")




















































    #roc = HuhCrisprROC("prop_file_path", "hits_file_path")
    #roc_df = roc.calc_metric()
    #partial_cov_interactors_prior_set()
    # cov_interactors_prior_set(method="iterative")
    # propagate_stefan()
    # compare_results(r"D:\data\propagations\stefan_cross_validation\self test\cov_interactors_analytical.json",
    #                 "D:\data\propagations\stefan_cross_validation\self test\cov_interactors_iterative.json")

    # nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
    #                                      r"D:\data\networks\krogan_cov_to_human.cx",
    #                                      r"D:\data\other\symbol_to_entrezgene.json",
    #                                      merge_covid=False).load()
    #
    # nw.remove_nodes_from([n for n, data in nw.nodes(data=True) if data["species_id"] == "sars-cov-2"])
    # NetpropNetwork.record_network(nw, r"D:\data\networks\no_covid.json")
    #
    #
    # with open(r"D:\configurations\stefan1\stefan_conf_backup.json", 'r') as handler:
    #     d = json.load(handler)
    # for prior_set in d["prior_set"]:
    #     for node in prior_set["nodes"]:
    #         node["source_of"] = ["info"]
    # with open(r"D:\configurations\stefan1\stefan_conf.json", 'w') as handler:
    #     json.dump(d, handler, indent=4)
    #propagate_stefan()

    # for conf_score in [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
    #     print(f"***  confidence {conf_score}  ***")
    #     partial_cov_interactors_prior_set(confidence=conf_score)
    #     propagate_stefan()
    #     files = [str(p) for p in Path(r"D:\data\propagations\stefan_cross_validation").glob("*.json")]
    #     intersect_top_propagated(files)
    #     print("\n\n\n")

    # for conf_score in [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
    #     print(f"***  confidence {conf_score}  ***")
    #     individual_cov_interactors_prior_set(confidence=conf_score)
    #     propagate_stefan()
    #     tops = top_prop_by_source(r"D:\data\propagations\stefan_unique_liquids\covid_interactors.json")
    #     top_intersected = most_intersected(tops)
    #     print(f"total liquids: {len(tops)}")
    #     histogram = {}
    #     for node, intersection_list in top_intersected.items():
    #         num = len(intersection_list)
    #         if num not in histogram:
    #             histogram[num] = 0
    #         histogram[len(intersection_list)] += 1
    #
    #     for num in sorted(list(histogram.keys())):
    #         print(f"{num}: {histogram[num]}")
    #
    #     print("\n\n\n")

    # specific_protein_partial_cov_interactors_prior_set()
    # stefan_propagations([str(file) for file in Path(r"D:\configurations\stefan1\specific_proteins").glob("*.json")])
    #sizes = set()
    #for file in [str(f) for f in Path(f"D:\configurations\stefan1\specific_proteins").glob("*.json")]:
        #with open(file, 'r') as handler:
        #    d = json.load(handler)
        #size = len(d["prior_set"][0]["nodes"]) + len(d["prior_set"][1]["nodes"])
        #print(f"{file}: {size}")
        #sizes.add(size)

    #nsp4: 8, nsp7: 32, nsp8: 24
    # for size in [9, 11, 15, 20, 26, 30, 32, 40, 47]:
    #     random_subgroups_partial_cov_interactors_prior_set(size, 100)

    #
    # files = []
    # processes = []
    # for triplet in [(40, 47)]:
    #     for folder in [Path(f"D:\\configurations\\stefan1\\random_subgroups\{num}") for num in triplet]:
    #         processes.append(Process(target=stefan_propagations, args=([str(f) for f in Path(folder).glob("*.json")],)))
    #     #     for file in [str(f) for f in Path(folder).glob("*.json")]:
    #     #         files.append(file)
    #     # stefan_propagations(files)
    #     for p in processes:
    #         p.start()
    #     for p in processes:
    #         p.join()



    # for folder in Path(r"D:\data\propagations\stefan_cross_validation\individual_proteins").glob("*"):
    #     files = [str(p) for p in folder.glob("*chunk*")]
    #     print(str(folder))
    #     more_detailed_intersection_data(files, output_file=str(folder / "detailed_intersection_data.json"))


    # for folder in Path(r"D:\data\propagations\stefan_cross_validation\random_subgroups").glob("*"):
    #     for subfolder in folder.glob("*"):
    #         files = [str(p) for p in subfolder.glob("*chunk*")]
    #         print(str(subfolder))
    #         more_detailed_intersection_data(files, output_file=str(subfolder / "detailed_intersection_data.json"))

    # p_value_dict = calc_protein_p_value(r"D:\configurations\stefan1\specific_proteins",
    #                                     r"D:\data\propagations\stefan_cross_validation\individual_proteins")
    # for k, v in p_value_dict.items():
    #     print(f"{k}: {v}")


    # intersection_dict = get_intersection_quality_metrics(r"D:\data\propagations\stefan_cross_validation\individual_proteins")
    # intersection_bar_plot(intersection_dict)

    # folders = [r"D:\data\propagations\stefan_cross_validation\random_subgroups\nsp4",
    #            r"D:\data\propagations\stefan_cross_validation\nsp7",
    #            r"D:\data\propagations\stefan_cross_validation\nsp8"]
    # for folder in folders:
    #     protein = folder.split("\\")[-1]
    #     files = [str(p) for p in Path(folder).glob("*.json")]
    #
    #     intersect_top_propagated(files,
    #                              output_file=f"D:\\data\\propagations\\stefan_cross_validation\\analysis_{protein}.json")
    # root_folders = [folder_path for folder_path in Path(r"D:\data\propagations\stefan_cross_validation\random_subgroups").glob("*") if folder_path.is_dir()]
    # folder_lists = [
    #     [str(subfolder) for subfolder in folder.glob("*")] for folder in root_folders
    # ]
    # for i in range(len(folder_lists)):
    #     print(f"{root_folders[i]}: {average_intersection(folder_lists[i])}")
    # average_intersection()
    # for folder in Path(r"D:\data\propagations\stefan_cross_validation\random_subgroups").glob("*"):
    #     files = [str(p) for p in folder.glob("*.json")]
    #     intersect_top_propagated(files, output_file=r"D:\data\propagations\stefan_cross_validation\random_subgroups\analysis8_16_24.json")

    # from utils.new_file_loaders import CovToHumanMeta
    # fixed_huh7_nw = CovToHumanMeta(r"C:\studies\thesis\code\cicd\data_cache\H_sapiens_aug_2020.net",
    #                                r"D:\data\networks\metastudy\protein_interactome_translated.csv",
    #                                r"D:\data\networks\metastudy\rna_interactome_translated.csv",
    #                                rna_cell_lines=["HUH7", "HUH7.5"], protein_cell_lines="all", min_num_studies=2).load()
    # NetpropNetwork.record_network(fixed_huh7_nw,  r"D:\data\networks\metastudy\huh7_only_fixed.json")


