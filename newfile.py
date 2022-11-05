import pandas as pd
import numpy as np
from crispr import get_huh_crispr
import  matplotlib.pyplot as plt
from scripts.divide_and_conquer2 import create_intersections_data, randomly_split, random_cross_proteins_splits, split_by_size
from multiprocessing import pool, Process, cpu_count
from functools import partial
from random_sets import randomly_replace_interactors, patch_conf_with_prot_to_interactor
from random_sets import randomly_replace_interactors
from utils.queue_managers import dump_json
from pathlib import Path


def intersect_top_k_hits(crispr_path: str, res_df_path: str):
    crispr_hits = [int(node) for node in get_huh_crispr(crispr_path)]

    res_df = pd.read_csv(res_df_path, index_col="nodes")
    res_df = res_df[[c for c in res_df.columns if not c.isnumeric()]]
    top_hits_per_metric = dict()
    for i in range(1, 5):
        if i < 4:
            s = res_df.apply(lambda row: sum(abs(row) ** i) ** (1. / i), axis=1).sort_values(ascending=False)
        else:
            s = res_df.apply(lambda row: max(row), axis=1).sort_values(ascending=False)

        sorted_hits = [hit for hit in s.index if hit in crispr_hits]
        label = f"L{i}" if i < 4 else "L_inf"
        top_hits_per_metric[label] = sorted_hits

    return top_hits_per_metric


def plot_intersecting_hits(top_hits_per_metric: dict):
    fig, axes = plt.subplots(nrows=2, ncols=3)
    axes = axes.flatten()
    labels = list(top_hits_per_metric.keys())
    series_length = len(top_hits_per_metric[labels[0]])
    x = np.arange(series_length) / series_length
    ax_loc = 0
    for i in range(len(labels) - 1):
        label_1 = labels[i]
        s1 = top_hits_per_metric[label_1]
        for j in range(i + 1, len(labels)):
            label_2 = labels[j]
            s2 = top_hits_per_metric[label_2]
            y = []
            intersections = 0
            discovered_hits = {n: 0 for n in s1}
            for i in range(series_length):
                for node in (s1[i], s2[i]):
                    discovered_hits[node] += 1
                    if discovered_hits[node] == 2:
                        intersections += 1
                y.append(intersections)
            y = np.array(y) / series_length

            ax = axes[ax_loc]
            ax.plot(x, y)
            ax.set_title(f"{label_1}, {label_2}")
            ax.set_xlabel("% hits discovered")
            ax.set_ylabel("% intersection")
            ax_loc += 1
    axes = axes.reshape((2, 3))
    plt.show()




## DEPRECATED ##
# def multiprocess_create_intersection_data(res_file: str, splits_file: str, output_prefix: str, thresholds: list[int]):
#     _create_funcs = [partial(create_intersections_data, res_file, splits_file, output_prefix + f"_{t}.json", t)
#                      for t in thresholds]
#     # processes = [Process(target=f) for f in _create_funcs]
#     # for p in processes:
#     #     p.start()
#     # for p in processes:
#     #     p.join()
#     for f in _create_funcs:
#         f()


def intersection_data():

    res = r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv"
    output_prefix = r"D:\data\propagations\krogan_interactors\individual_interactors\split2_results\split_results_as_percentages_top"
    rand_res_0 = r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\all\no_knockouts.csv"
    rand_0_output_prefix = r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\split2_results\split_results_as_percentages_top"
    rand_res_90 = r"D:\data\propagations\krogan_interactors\individual_interactors\randomized90\all\no_knockouts.csv"
    rand_90_output_prefix = r"D:\data\propagations\krogan_interactors\individual_interactors\randomized\split2_results\split_results_as_percentages_top"
    splits_file = r"D:\data\propagations\krogan_interactors\individual_interactors\splits.json"
    thresholds = [10*i for i in range(1, 11)] + [150]

    random_deg_preserving_sets_res = r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\no_knockouts.csv"
    random_deg_preserving_output_prefix = r"D:\data\propagations\krogan_interactors\individual_interactors\all_split_results\random_set_split_results\split_results_as_percentages_top"
    random_deg_preserving_splits_file = r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\rand_sets_splits.json"

    # multiprocess_create_intersection_data(res, splits_file, output_prefix, thresholds)
    # multiprocess_create_intersection_data(rand_res_0, splits_file, rand_0_output_prefix, thresholds)
    # multiprocess_create_intersection_data(rand_res_90, splits_file, rand_90_output_prefix, thresholds)
    create_intersections_data(random_deg_preserving_sets_res,
                              random_deg_preserving_splits_file,
                              random_deg_preserving_output_prefix,
                              [20, 100])


def pipeline(metadata_file: str, res_file: str, output_root: str, split_repetitions: int, intersection_thresholds: list[int]):
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
        # (inter_output_dir, partial(randomly_split, metadata_file, None), "inter"),
        # (cross_output_dir, partial(random_cross_proteins_splits, res_file, metadata_file), "cross"),
        # (by_size_output_dir, partial(split_by_size, res_file, [10, 15, 20, 30, 40]), "by_size"),
        (all_interactors_output_dir, partial(randomly_split, metadata_file, ["all"]), "all")
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



if __name__ == "__main__":
    input_tuples = [
        # (r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_metadata.json",
        #  r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactor_sets\all\no_knockouts.csv",
        #  r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\randomized_sets",
        #  20,
        #  [20, 100]),
        (r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json",
         r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
         r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\real",
         20,
         [20, 100])
        # (
        #     r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json",
        #     r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv",
        #     r"D:\data\propagations\krogan_interactors\individual_interactors\splits_statistics\from_all_proteins",
        #     20,
        #     [20, 100]
        # ),
    ]

    for tup in input_tuples:
        pipeline(*tup)


    # intersection_data()
    # top_hits_metric = intersect_top_k_hits(r"D:\data\other\huh7_crispr_translated.csv",
    #                                        r"D:\data\propagations\krogan_interactors\individual_interactors\all.csv")
    # plot_intersecting_hits(top_hits_metric)
    # randomly_replace_interactors(r"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json",
    #                              r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_metadata.json")

    # patch_conf_with_prot_to_interactor(r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_metadata.json",
    #                                    r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_conf.json",
    #                                    r"D:\data\propagations\krogan_interactors\individual_interactors\randomized_interactors",
    #                                    r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_randomized_sets_conf.json")




