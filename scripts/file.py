from netprop.networks.randomization import randomize_network
from netprop.networks.loaders import NetpropNetwork, MetaCovidHumanLoader, CombinedHumanCovidNetworkLoader
from utils.new_file_loaders import CovToHumanMeta
from crispr import get_huh_crispr, get_crispr_rankings
from netprop.models import PropagationResultModel
import json
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def gen_networks(human_ppi_path: str):
    # generate 101 (original + 100 randomized) of networks with different combinations of PPIs
    krogan_nw = CombinedHumanCovidNetworkLoader(human_ppi_path,
                                                "/data/ppi/cov_to_human_ppi.cx", "/data/symbol_to_entrezgene",
                                                merge_covid=False).load()
    full_metadata_nw = CovToHumanMeta(human_ppi_path,
                                      protein_interactions_path="/data/ppi/protein_interactome_translated.csv",
                                      protein_cell_lines="all", merged_covid=False).load()

    huh_crispr_set = get_huh_crispr('/data/huh7_crispr_translated.csv')
    krogan_crispr_scores = get_crispr_rankings('/data/crispr_screen.csv', "/data/symbol_to_entrezgene")


    norm_methods = ["RPDN", "prop_1_from_all", "the other weird one I don't quite get yet (TM)"]



def intersect_top_propagated_with_pvalue(results: list[str], k: int=50, output_file: str=None):
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


def get_p_value(num_proteins: int, sorted_intersection_scores: list[float]):
    root_folder = Path(r'D:\data\propagations\stefan_cross_validation\random_subgroups')
    folders = list((root_folder / f"{num_proteins}").glob("*"))
    better_than_real = 0
    for folder in folders:
        with open(str(folder / "detailed_intersection_data.json"), 'r') as handler:
            folder_data = json.load(handler)
        chunk_scores = sorted([chunk["size"] for chunk in folder_data["chunks"].values()])
        # better_than_real += int(
        #     sum(
        #         [int(chunk_scores[i] >= sorted_intersection_scores[i]) for i in range(len(chunk_scores))]
        #     ) >= 1
        # )
        # better_than_real += len([1 for i in range(len(chunk_scores)) if chunk_scores[i] > sorted_intersection_scores[i]])
        if sum(chunk_scores) > sum(sorted_intersection_scores):
            better_than_real += 1
    # return better_than_real / (len(sorted_intersection_scores) * len(folders))
    return better_than_real / len(folders)



def get_intersection_quality_metrics(root_folder: str):
    protein_map = {}
    for protein_folder in Path(root_folder).glob("*"):
        protein_name = protein_folder.name
        with open(str(protein_folder / "detailed_intersection_data.json"), 'r') as handler:
            folder_data = json.load(handler)
        protein_map[protein_name] = sorted([chunk["size"] for chunk in folder_data["chunks"].values()])

    return protein_map


def get_protein_interactor_num(root_folder: str):
    protein_map = dict()
    for file in [str(f) for f in Path(f"D:\configurations\stefan1\specific_proteins").glob("*.json")]:
        protein_name = Path(file).name.split(".")[0].split("_")[-1]
        with open(file, 'r') as handler:
           d = json.load(handler)
        size = len(d["prior_set"][0]["nodes"]) + len(d["prior_set"][1]["nodes"])
        print(f"{file}: {size}")
        protein_map[protein_name] = size
    return protein_map


def calc_protein_p_value(conf_root_folder: str, prop_root_folder: str):
    prot_to_num_interactors = get_protein_interactor_num(conf_root_folder)
    prot_to_intersection_metrics = get_intersection_quality_metrics(prop_root_folder)
    p_values = dict()
    for protein in prot_to_num_interactors:
        num_interactors = prot_to_num_interactors[protein]
        intersection_metric = prot_to_intersection_metrics[protein]
        p_values[protein] = get_p_value(num_interactors, intersection_metric)

    return p_values



def intersection_bar_plot(chunk_intersection_dict: dict[str, list[float]]):
    protein_names = []
    chunk_values = []

    # for k,v in chunk_intersection_dict.items():
    #     protein_names.append(k)
    #     for i in range(len(v)):
    #         chunk_values[i].append(v[i])
    #
    # num_chunks = 5
    # for i in range(num_chunks):
    #     for k, v in chunk_intersection_dict.items():

    for k, v in chunk_intersection_dict.items():
        protein_names.append(k)
        chunk_values.append(v)
    chunk_values_as_cols = list(np.array(chunk_values).T)
    d = {"interactors of": protein_names}
    for i in range(len(chunk_values_as_cols)):
        d[f"random split {i}"] = chunk_values_as_cols[i]
    # df = pd.DataFrame([protein_names] + chunk_values_as_cols,
    #                   columns=["interactors of"] + [f"chunk {i}" for i in range(len(chunk_values_as_cols))])
    df = pd.DataFrame(d)
    ax = df.plot(x="interactors of", y=[f"random split {i}" for i in range(len(chunk_values_as_cols))], kind="bar")
    ax.set_ylabel("% intersection out of top 50")
    plt.title("intersections of top scoring proteins when propagated from randomly split sets of interactors with covid proteins")
    plt.show()

