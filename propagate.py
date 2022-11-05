import json
from netprop.propagation import propagate_from_config
from utils.new_file_loaders import CovToHumanMeta, NetpropCovToHumanLoader
from netprop.networks.randomization import randomize_network
from netprop.networks.loaders import NetpropNetwork, HSapiensNetworkLoader
from netprop.models import ConfigModel, PropagationResultModel
from copy import deepcopy
from metrics.correlation import CorrelationMetric
import pandas as pd
from scripts.divide_and_conquer2 import create_intersections_data, kinda_correlate_interesection_to_prop, \
    intersect_top_propagated_with_cripsr_kos, crispr_targets_in_top_intersections, randomly_split, PropResNewFormat, random_cross_proteins_splits, split_by_size
from pathlib import Path
from utils.queue_managers import dump_json, load_json
from itertools import product, chain
from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
from functools import partial

def histogram_to_json(path, outpath):
    with open(path, 'r') as handler:
        lines = [line.replace('[', '').replace(']', '').replace(',', '').replace('\'', '').strip() for line in handler.readlines()]
    dictogram = {}
    for i in range(len(lines)):
        knockouts = [e.split('_')[0] for e in lines[i].split()]
        dictogram[i] = knockouts

    with open(outpath, 'w') as handler:
        json.dump(dictogram, handler, indent=4)


def construct_network():
    nw = CovToHumanMeta(r"D:\data\networks\H_sapiens_aug_2020.net",
                        r"D:\data\networks\metastudy\protein_interactome_translated.csv",
                        r"D:\data\networks\metastudy\rna_interactome_translated.csv",
                        rna_cell_lines=['Huh7', 'Huh7.5'], protein_cell_lines="all").load()
    NetpropNetwork.record_network(nw, r"D:\data\networks\metastudy\complete_networks\huh7_only.json")




def liq_per_cov_protein_all_interactors(krogan_separated_conf: str, output_path: str) -> None:
    ks_conf = ConfigModel.parse_file(krogan_separated_conf)
    output_conf = deepcopy(ks_conf)

    output_conf.prior_set = []
    combined_ps = {"id": "all", "confidence": 0.7, "nodes": []}
    nodes_to_liquids = {}
    for ps in ks_conf.prior_set:
        for node in ps.nodes:
            if node not in nodes_to_liquids:
                nodes_to_liquids[node.id] = set()
            nodes_to_liquids[node.id].add(ps.id.split("_")[0])
            nodes_to_liquids[node.id].add(node.id)

    combined_ps["nodes"] = [{"id": node, "source_of": list(data)} for node, data in nodes_to_liquids.items()]
    output_conf.prior_set.append(combined_ps)
    with open(output_path, 'w') as handler:
        json.dump(output_conf.dict(exclude_unset=True), handler, indent=4)


def prop_result_to_big_ass_df(res_path: str, output_path: str):
    res = PropagationResultModel.parse_file(res_path).prop_scores_as_series(by_liquids=None)
    res = res.set_index("nodes")
    res = res.rename(columns={c: c.split(".")[-1] for c in list(res.columns)})
    res.to_csv(output_path)


def metadata_file(conf_path: str, metadata_path:str):
    conf = ConfigModel.parse_file(conf_path)
    metadata = dict()
    for ps in conf.prior_set:
        nodes = [elem.id for elem in ps.nodes]
        protein = ps.id.split("_")[0]
        metadata[protein] = nodes
    with open(metadata_path, 'w') as handler:
        json.dump(metadata, handler, indent=4)


def metadata_file_new(conf_path: str, metadata_path: str):
    conf = ConfigModel.parse_file(conf_path)
    prior_nodes = conf.prior_set[0].nodes
    prot_to_nodes = dict()
    for n in prior_nodes:
        proteins = [x for x in n.source_of if not x.isnumeric()]
        for p in proteins:
            if p not in prot_to_nodes:
                prot_to_nodes[p] = []
            prot_to_nodes[p].append(n.id)

    dump_json(prot_to_nodes, metadata_path)


def redo_csv_from_metadata(metadata_path: str, csv_path: str, new_csv_path: str):
    metadata = load_json(metadata_path)
    df = pd.read_csv(csv_path)
    for prot_name in metadata:
        if prot_name not in df.columns:
            raise Exception(f"metadata protein {prot_name} not in df")
        df[prot_name] = df[metadata[prot_name]].sum(axis=1)
    df.to_csv(new_csv_path)


def inter_protein_correlation(prop_csv_path: str, output_dir_path: str, threshold: int=None):
    # prot_to_nodes = load_json(r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json")
    prot_to_nodes = load_json(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_metadata.json")
    df = pd.read_csv(prop_csv_path)
    inter_protein_correlations = dict()
    if threshold is None:
        threshold = len(df) + 1
    for prot, interactors in prot_to_nodes.items():
        prot_df = df[["nodes", prot] + interactors].sort_values(prot, ascending=False).iloc[:int(threshold)]
        prot_df = prot_df.drop(prot, axis=1)
        inter_protein_correlations[prot] = CorrelationMetric(prot_df).calc_metric()
    #
    #
    #
    #
    # inter_protein_correlations = {
    #     prot: CorrelationMetric(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\all\no_knockouts.csv", interactors).calc_metric()
    #     for prot, interactors in prot_to_nodes.items()
    # }

    root_folder = Path(output_dir_path)
    root_folder.mkdir(exist_ok=True)
    for prot, df in inter_protein_correlations.items():
        df.to_csv(root_folder / f"{prot}.csv")


def _cross_correlation_to_csv(prop_res_path: str, output_folder_path: str, threshold: int, prots):
    root_folder = Path(output_folder_path)
    root_folder.mkdir(exist_ok=True)
    prot1, prot2 = prots[0]
    interactors1, interactors2 = prots[1]
    interactors1, interactors2 = prots[1]
    df = pd.read_csv(prop_res_path)
    df["combined_score"] = df[prot1] + df[prot2]
    if threshold is None:
        threshold = len(df) + 1
    df = df[["nodes", "combined_score"] + interactors1 + interactors2].sort_values("combined_score", ascending=False).iloc[:int(threshold)]
    df = df.drop("combined_score", axis=1)
    corr_table = CorrelationMetric(df,
                                   interactors1 + interactors2).calc_metric()
    only_cross_section = corr_table.loc[interactors1, interactors2]

    only_cross_section.to_csv(root_folder / f"{prot1}_{prot2}.csv")


def cross_protein_correlation(metadata_path: str, prop_res_path: str, output_folder_path: str, threshold: int=None):
    prot_to_nodes = load_json(metadata_path)
    prot_pairs = []
    proteins = list(prot_to_nodes.keys())
    for i in range(len(proteins) - 1):
        prot1 = proteins[i]
        for j in range(i + 1, len(proteins)):
            prot2 = proteins[j]
            prot_pairs.append(((prot1, prot2), (prot_to_nodes[prot1], prot_to_nodes[prot2])))

    map_func = partial(_cross_correlation_to_csv, prop_res_path, output_folder_path, threshold)
    with Pool(7) as p:
        p.map(map_func, prot_pairs)


def very_basic_table_diff(csv1: str, csv2: str, table_type: str) -> float:
    t1 = pd.read_csv(csv1, index_col=0).sort_index().sort_index(axis=1)
    t2 = pd.read_csv(csv2, index_col=0).sort_index().sort_index(axis=1)
    mask = ((t1 - t2) > 0).fillna(0)
    if table_type == "inter":

        distinct_elements_in_table = (len(t1) ** 2 - len(t1)) / 2   # num elements above diagonal
        return mask.to_numpy().sum() / (2*distinct_elements_in_table)   # divide by 2 to not count lower than diagonal
    elif table_type == "cross":
        return mask.to_numpy().sum()/ (len(t1)*len(t1.columns))


def compare_tables(file_batch_1: list[str], file_batch_2: list[str], table_type: str) -> pd.Series:
    s = pd.Series(dtype="float64")
    fb1 = sorted(file_batch_1)
    fb2 = sorted(file_batch_2)
    for i in range(len(file_batch_1)):
        try:
            f1, f2 = fb1[i], fb2[i]
        except:
            x = 7
        s[Path(f1).stem] = very_basic_table_diff(f1, f2, table_type)
    return s


def _scores_from_file_list(file_list: list[str], abs_vals: bool=False) -> np.ndarray:
    scores = []
    for f in file_list:
        df = pd.read_csv(f, index_col=0).fillna(0)
        val_list = chain(*df.values.tolist())
        scores.extend(val_list)

    scores = np.array(sorted(scores))
    scores = scores[scores != 0]
    if abs_vals:
        scores = np.abs(scores)

    return scores


def _calc_integral(scores: np.ndarray) -> np.ndarray:
    step_size = 1e-5
    latest_index = 0
    integral = np.array([0] * int(1e5))
    ticks = np.arange(0, 1, 1e-5)
    for i in range(len(ticks)):
        threshold = step_size * i
        while latest_index < len(scores) and scores[latest_index] < threshold:
            latest_index += 1
        integral[i] = latest_index

    integral = integral / len(scores)
    return integral


def _calc_bins(scores: np.ndarray, num_bins: int):
    thresholds = [(-1 + i*(2/num_bins), -1 + (i+1)*(2/num_bins)) for i in range(num_bins)]
    bins = np.array([np.sum(np.logical_and(t[0] < scores, scores < t[1])) / len(scores) for t in thresholds])
    return bins


def compare_inter_to_cross(inter_files: list[str], cross_files: list[str], method: str) -> None:
    inter_scores = _scores_from_file_list(inter_files)
    cross_scores = _scores_from_file_list(cross_files)
    plt.title("spearman correlation value distribution [top 1k]")
    plt.xlabel("correlation value")
    plt.ylabel("% of all pairwise correlations")
    if method == "integral":
        inter_scores_integral = _calc_integral(inter_scores)
        cross_scores_integral = _calc_integral(cross_scores)
        ticks = np.arange(0, 1, 1e-5)


        plt.plot(inter_scores_integral, ticks, c="b", label="inter protein")
        plt.plot(cross_scores_integral, ticks, c="r", label="cross protein")

    elif method == "bins":
        num_bins = 20
        thresholds = [-1 + i*2/num_bins for i in range(num_bins + 1)]
        X = [f"{thresholds[i]:.1f} - {thresholds[i + 1]:.1f}" for i in range(len(thresholds) - 1)]
        X_axis = np.arange(len(X))
        plt.bar(X_axis + 0.2, _calc_bins(inter_scores, num_bins), 0.4, label="inter protein")
        plt.bar(X_axis - 0.2, _calc_bins(cross_scores, num_bins), 0.4, label="cross protein")
        plt.xticks(X_axis, X, rotation=90)

    plt.legend()
    plt.show()




def visualize_intersections(parent_folders: list[str], categroy_labels: list[str], thresholds: list[str]):
    nrows, ncols = len(thresholds), len(parent_folders)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
    axes = axes.reshape((nrows, ncols))
    if len(thresholds) * len(parent_folders) == 1:
        axes = np.array([axes])
    split_results = dict()
    for j in range(len(parent_folders)):
        parent_folder = parent_folders[j]
        split_results[parent_folder] = dict()
        for t in thresholds:
            res = load_json(str(Path(parent_folder) / f"split_results_as_percentages_top_{t}.json"))
            keys_for_all_proteins = set().union(*[set(v.keys()) for v in res.values()])

            split_results[parent_folder].update({t: {"res": res, "intersected": keys_for_all_proteins}})

    x_labels_defined = False
    for i in range(len(thresholds)):
        t = thresholds[i]
        y_labels = list(set().union(*[split_results[pf][t]["intersected"] for pf in parent_folders]))
        y_label_to_idx = {y_labels[i - 1]: i for i in range(1, len(y_labels) + 1)}


        for j in range(len(parent_folders)):
            pf = parent_folders[j]
            ax = axes[i,j]
            res = split_results[pf][t]["res"]

            # if not x_labels_defined:
            x_labels = list(split_results[pf][t]["res"].keys())
            x_label_to_idx = {x_labels[i - 1]: i for i in range(1, len(x_labels) + 1)}
                # x_labels_defined = True

            x, y, z = [], [], []
            for x_label,data in res.items():
                for y_label, value in data.items():
                    y.append(y_label_to_idx[y_label])
                    x.append(x_label_to_idx[x_label])
                    z.append(value)
            z = np.array(z)
            ax.scatter(x,y,s=z*1000, alpha=0.5)
            ax.set_yticks(np.arange(len(y_labels) + 1), ["-"] + y_labels, fontsize=5)
            ax.set_xticks(np.arange(len(x_labels) + 1), ["-"] + x_labels, fontsize=10)


    plt.figtext(0.05, 0.95, "x axis: protein for which interactors were split\n"
                            "y axis: protein that appeared in intersections", multialignment="left", color="b")
    for i in range(axes.shape[0]):
        pos = axes[i, 0].get_position()
        row_center_y = (pos.y0 + pos.y1) / 2
        plt.figtext(0, row_center_y, f" intersections in top {thresholds[i]}")

    for i in range(axes.shape[1]):
        pos = axes[0, i].get_position()
        col_center_x = (pos.x0 + pos.x1) / 2
        plt.figtext(col_center_x, pos.y1 + 0.02, categroy_labels[i])
    plt.suptitle("bubble size proportional to % of intersections in which protein appeared")
    plt.show()



def create_corr_files(threshold: int):
    root = Path(r"D:\data\propagations\krogan_interactors\individual_interactors")
    name = f"top_{threshold}" if threshold is not None else "max"
    inter_non_random = (str(root/ "all.csv"), str(root / f"correlations\\inter_protein_{name}"))
    cross_non_random = (str(root/ "all.csv"),
                        str(root / f"correlations\\cross_protein_{name}"))
    inter_random = (str(root / r"randomized90\all\no_knockouts.csv"),
                    str(root / f"correlations\\automated\\inter_protein_{name}_randomized90"))

    cross_random = (str(root / r"randomized\all\no_knockouts.csv"),
                    str(root / f"correlations\\cross_protein_{name}_randomized"))

    inter_randomized_interactor_sets = (str(root / "randomized_interactor_sets/all/no_knockouts.csv"),
                                        str(root / f"correlations\\automated\\inter_protein_{name}_randomized_sets"))

    # tuples = (inter_non_random, cross_non_random, inter_random, cross_random)
    # inter_tuples = (inter_non_random, inter_random)
    inter_tuples = (inter_randomized_interactor_sets,)

    cross_tuples = (cross_non_random, cross_random)
    for tup in inter_tuples:
        inter_protein_correlation(tup[0], tup[1], threshold)


    # for tup in cross_tuples:
    #     cross_protein_correlation(r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json", tup[0], tup[1], threshold)



def plot_correlation_dominance(thresholds: list[int]):
    # for t in thresholds:
    #     create_corr_files(t)

    inter_res_per_threshold, cross_res_per_threshold = [], []

    for t in thresholds:
        name = f"top_{int(t)}" if t is not None else "max"

        corr_folders_root = Path(r"D:\data\propagations\krogan_interactors\individual_interactors\correlations\automated")
        non_random_inter = corr_folders_root / f"inter_protein_{name}"
        non_random_cross= corr_folders_root / f"cross_protein_{name}"
        random_inter = corr_folders_root / f"inter_protein_{name}_randomized"
        random_cross = corr_folders_root / f"cross_protein_{name}_randomized"

        get_files = lambda folder: [str(f) for f in folder.glob("*")]
        s_inter = compare_tables(get_files(non_random_inter), get_files(random_inter), "inter")
        s_cross = compare_tables(get_files(non_random_cross), get_files(random_cross), "cross")
        try:
            inter_res_per_threshold.append(len(s_inter[s_inter >= 0.5]) / len(s_inter))
            cross_res_per_threshold.append(len(s_cross[s_cross >= 0.5]) / len(s_cross))
        except:
            print(f"problem: {name}")
            x = 7

    thresholds[-1] = "no thresholds (all)"
    x = np.arange(len(thresholds))
    plt.bar(x + 0.2, inter_res_per_threshold, label="inter",  width=0.4)
    plt.bar(x - 0.2, cross_res_per_threshold, label="cross", width=0.4)
    plt.legend()
    plt.xticks(x, thresholds)
    plt.xlabel("# top propagated results")
    plt.ylabel("% proteins where real correlation was higher than random")
    plt.title("for any threshold, pairwise correlation is weaker in real graph than in randomized one")
    plt.show()




def main():
    # root = Path(r"D:\data\propagations\krogan_interactors\individual_interactors\correlations")
    # corr_files = [
    #         str(f) for f in
    #         (root / "inter_protein_top_100").glob("*")
    #     ]
    # randomized_corr_files = [
    #         str(f) for f in
    #         (root / "inter_protein_top_100_randomized").glob("*")
    #     ]
    # #
    # cross_corr_files = [
    #         str(f) for f in
    #         (root / "cross_protein_top_100").glob("*")
    #     ]
    #
    # randomized_cross_corr_files = [
    #         str(f) for f in
    #         (root / "cross_protein_top_100_randomized").glob("*")
    #     ]
    #
    # create_corr_files(10000)
    # create_corr_files(None)
    # plot_correlation_dominance([10, 100, 500, 1000, int(1e4), None])
    create_corr_files(100)

    # s = compare_tables(cross_corr_files, randomized_cross_corr_files, "cross")
        # s = compare_tables(corr_files, randomized_corr_files)

    # cross_protein_correlation(r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json",
    #                           r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\all\no_knockouts.csv",
    #                           r"D:\data\propagations\krogan_interactors\individual_interactors\correlations\cross_protein_1k_randomized")

    # compare_inter_to_cross(corr_files, cross_corr_files, "bins")

    x = 7





    # nw = HSapiensNetworkLoader( r"D:\data\networks\randomized_h_sapiens\h_sapiens_90").load()
    # NetpropNetwork.record_network(nw,  r"D:\data\networks\randomized_h_sapiens\h_sapiens_90.json")
    # propagate_from_config(r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_randomized_conf.json")

    # prop_result_to_big_ass_df(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized90\all\no_knockouts.json",
    #                           r"D:\data\propagations\krogan_interactors\individual_interactors\randomized90\all\no_knockouts.csv")





if __name__ == "__main__":
    # metadata_file_new(r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_conf.json",
    #                   r"D:\data\propagations\krogan_interactors\individual_interactors\fixed_metadata.json")

    #
    # df = pd.read_csv(r"D:\data\propagations\krogan_interactors\individual_interactors\randomly_allocated_all.csv")
    # df = df.drop("Unnamed: 0", axis=1)
    # df.to_csv(r"D:\data\propagations\krogan_interactors\individual_interactors\randomly_allocated_all.csv", index=False)

    # split_by_size(r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
    #               [10, 15, 20,30,40],
    #               r"D:\data\propagations\krogan_interactors\individual_interactors\interactors_from_all_proteins\metadata.json")
    # redo_csv_from_metadata(r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json",
    #                         r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
    #                         r"D:\data\propagations\krogan_interactors\individual_interactors\randomly_allocated_all.csv")
    # intersect_top_propagated_with_cripsr_kos(r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
    #                                         r"D:\data\other\huh7_crispr_translated.csv", "no")
    # prop_result_to_big_ass_df(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\no_knockouts.json",
    #                           r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\no_knockouts.csv")
    # res = pd.read_csv(r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv")
    #
    # proteins = [c for c in res.columns if not c.isnumeric() and c != "nodes"]
    # res["all"] = res.apply(lambda row: row.sum(), axis=1)
    #
    # propr_res = PropResNewFormat(r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
    #                              r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json")
    # splits = randomly_split(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_metadata.json",
    #                         None)
    # splits = random_cross_proteins_splits(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\no_knockouts.csv",
    #                                       r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_metadata.json")
    # dump_json(splits, r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\rand_sets_cross_splits.json")
    # visualize_intersections([r"D:\data\propagations\krogan_interactors\individual_interactors\split_results",
    #                          r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\split_results",
    #                          r"D:\data\propagations\krogan_interactors\individual_interactors\cross_protein_split_results"], ["real", "randomized", "cross protein"], [20, 100])
    
    main()


    # nw = HSapiensNetworkLoader(r"D:\\data\\networks\\randomized_h_sapiens\\h_sapiens_3").load()
    # NetpropNetwork.record_network(nw, "D:\\data\\networks\\randomized_h_sapiens\\h_sapiens_3.json")


    # propagate_from_config(r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_randomized_conf.json")
    # prop_result_to_big_ass_df(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\all\no_knockouts.json",
    #                           r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\all\no_knockouts.csv")
    # d = {}
    # for intersection_file in Path(r"D:\data\propagations\krogan_interactors\individual_interactors\split_results").glob("*.json"):
    #     top_k = int(intersection_file.stem.split("_")[-1])
    #     d[top_k] = crispr_targets_in_top_intersections(str(intersection_file), r"D:\data\other\huh7_crispr_translated.csv")
    #
    # dump_json(d, r"D:\data\propagations\krogan_interactors\individual_interactors\splits_vs_crispr.json")



    # res_dict = \
    # kinda_correlate_interesection_to_prop(r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
    #                                       "D:\data\propagations\krogan_interactors\individual_interactors\split_results_as_percentages.json",
    #                                       50)

    #
    # res_dict = \
    #     intersect_top_propagated_with_cripsr_kos(r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
    #                                              r"D:\data\other\huh7_crispr_translated.csv",
    #                                              10)
    # for k, v in res_dict.items():
    #     print(f"{k}: {v}")