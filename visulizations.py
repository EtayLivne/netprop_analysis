from cProfile import label
from random import shuffle
import matplotlib.pyplot as plt
from utils.utils import load_json, dump_json
from pathlib import Path
import numpy as np
from crispr import get_huh_crispr
import pandas as pd


def _extract_intersection_values(intersection_results: dict[str, dict[str, float]]) -> list[float]:
    results = []
    for intersection_data in intersection_results.values():
        results.extend(intersection_data.values())
    if not results:
        results = [0]
    return results


def _merge_split_results(folder: str, file_name_pattern: str=None) -> tuple[list, list]:
    if file_name_pattern is None:
        file_name_pattern = "*"
    
    files = [str(f) for f in Path(folder).glob(file_name_pattern)]
    all_nodes, all_values = [], []
    for f in files:
        intersections = load_json(f)
        for intersection_data in intersections.values():
            all_nodes.extend(intersection_data.keys())
            all_values.extend(intersection_data.values())
    all_nodes = list(set(all_nodes))
    all_values = list(set(all_values))
    return all_nodes, all_values




def visualize_intersection_statistics(root_folder: str, t_values: list[int]):
    root = Path(root_folder)
    folders = list(root.glob("*"))
    th_values = {v: [[] for i in range(10)] for v in t_values}
    bins = [(0.05 * i, 0.05 * (i+1)) for i in range(20)]
    for t in th_values.keys():
        for j in range(len(folders)):
            folder = folders[j]
            t_files_in_folder = list(folder.glob(f"*{str(t)}*"))
            if len(t_files_in_folder) != 1:
                raise Exception("incorrect num files in folder")
            file = t_files_in_folder[0]
            intersections = load_json(str(file))
            intersection_values = sorted(_extract_intersection_values(intersections))
            local_bins = [0 for i in range(20)]
            for i in range(20):
                while intersection_values and intersection_values[0] >= bins[i][0] and intersection_values[0] <= bins[i][1]:
                    local_bins[i] += 1
                    intersection_values.pop(0)
            if len(intersection_values) > 0:
                raise Exception("how could a value not be between 0 and 1?")
                
            th_values[t][j] = local_bins
    
    fig, axes = plt.subplots(len(th_values))
    for i in range(len(t_values)):
        ax = axes[i]
        t = t_values[i]
        list_of_val_lists = th_values[t]
        boxes = [[lst[k] for lst in list_of_val_lists] for k in range(len(list_of_val_lists[0]))]
        bin_vals = boxes
        ax.boxplot(bin_vals)
        ax.set_xticks(np.arange(len(bins)) + 1, labels=[f"{b[1]:.2f}" for b in bins])   
        ax.set_title(f"Top {t} nodes by propagation score")
        if i == len(t_values) - 1:
            ax.set_xlabel("appearance in % of all intersections")
        ax.set_ylabel("num of nodes")
    
    plt.suptitle("SARS COV 1: Frequency of a node appearing in multiple intersections")
    plt.show()



def visualize_intersection_crispr_prediction(root_folder: str, th_values: list[int], crispr_path: str):
    nodes_in_graph = 19100
    crispr_nodes = set(get_huh_crispr(crispr_path))
    root = Path(root_folder)
    folders = [str(folder) for folder in root.glob("*")]
    for t in th_values:
        hits_array = []

        for folder in folders:
            nodes, values = _merge_split_results(folder, file_name_pattern=f"*{str(t)}*")
            num_nodes = len(nodes)
            num_hits = len(crispr_nodes & set(nodes))
            expected_hits = num_nodes/nodes_in_graph
            hits_array.append(num_hits)
    
        print(f"{t}: {hits_array}")



def visualize_correlation(corr_files: list[str], corr_type):
    proteins = dict()
    for f in corr_files:
        values = pd.read_csv(f).drop('Unnamed: 0', axis=1).fillna(0).to_numpy()
        if len(values) < 10:
            continue
        if corr_type == "inter":
            values = np.triu(values)
        proteins[Path(f).stem] = values.flatten()
    
    x_labels = [l.replace("_", ", ") for l in list(proteins.keys())]
    x = np.arange(len(x_labels))
    data = list(proteins.values())
    plt.boxplot(data)
    plt.xticks(x + 1, labels=x_labels, rotation=90, fontsize=6)
    plt.title("cross protein pairwise correlation values in real network")
    plt.ylabel("spearman correlation")
    plt.show()





# TEMPORARY LOCATION MOVE TO DIVIDE&CONQUER2!
def randomly_allocate_interactors(prot_to_interactor: str) -> dict[str, list[str]]:
    p_to_t = load_json(prot_to_interactor)
    all_interactors = list(set().union(*[v for v in p_to_t.values()]))
    shuffle(all_interactors)
    new_p_to_t = dict()
    for k,v in p_to_t.items():
        new_p_to_t[k] = all_interactors[:len(v)]
        all_interactors = all_interactors[len(v):]
    
    return new_p_to_t


if __name__ == "__main__":
    root_path = r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\randomized_sets\intersection_results\inter"
    th_values = [20, 100]
    randomized_root = Path(r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\randomized_sets\intersection_results")
    real_root = Path(r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\real\intersection_results")
    randomized_sets_root = Path(r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\randomized_sets\intersection_results"
                                )
    randomized_root_cross = str(randomized_root / "cross")
    randomized_root_inter = str(randomized_root / "inter")
    real_root_cross = str(real_root / "cross")
    real_root_inter = str(real_root / "inter")
    real_root_by_size = str(real_root / "by_size")
    real_root_all_interactors = str(real_root / "all")
    randomized_sets_root_cross = str(randomized_sets_root / "cross")
    randomized_sets_root_inter = str(randomized_sets_root / "inter")
    randomized_sets_root_by_zise = str(randomized_sets_root / "by_size")


    # SARS cov 1 visualization
    cov1_root = Path(r"D:\data\propagations\sars_cov_1_individual_interactors\intersections_data\intersection_results")
    cov1_inter_root = str(cov1_root / "inter")
    cov1_size_root = str(cov1_root / "by_size")


    crispr_path =r"D:\data\other\huh7_crispr_translated.csv"

    visualize_intersection_statistics(cov1_size_root, th_values)



    # visualize_intersection_crispr_prediction(real_root_inter, [20, 100], crispr_path)

    # randomly_allocated = randomly_allocate_interactors(r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json")
    # dump_json(randomly_allocated, r"D:\data\propagations\krogan_interactors\individual_interactors\randomly_allocated_metadata.json")


    # real_inter_corr_files = [str(f) for f in Path(r"D:\data\propagations\krogan_interactors\individual_interactors\correlations\automated\cross_protein_top_100").glob("*")]
    # visualize_correlation(real_inter_corr_files, "cross")

