import json
import numpy as np
import pandas as pd
from netprop.models import PropagationResultModel
from pathlib import Path
import matplotlib.pyplot as plt
from multiprocessing import Queue, Pool
from time import sleep
from typing import Callable


def propagation_results_to_df(prop: PropagationResultModel):
    names = np.array(prop.nodes.keys())
    propagation_scores = np.array([prop.nodes[n]["liquids"]["info"] for n in names])
    return pd.DataFrame([names, propagation_scores], columns=["name", "propagation_score"])


def propagation_diff_df(prop1: PropagationResultModel, prop2: PropagationResultModel,
                        positive_filter: Callable = None) -> pd.Series:
    names = list(prop1.nodes.keys())
    if not positive_filter:
        positive_filter = lambda x: True
    try:
        return pd.Series({n: prop1.nodes[n].liquids["info"] - prop2.nodes.get(n, prop1.nodes[n]).liquids["info"]
                          for n in names if positive_filter(n)})
    except:
        print(prop1.network["suppressed_nodes"][0])
        return pd.Series([])


def get_scores_df(reference_propagation: str, tested_propagation: list[str], sort=True): # crispr_screen_csv: str
    reference_prop = PropagationResultModel.parse_file(reference_propagation)
    # crispr_df = get_crispr_enrichments(crispr_screen_csv, (r"D:\complete_translation.json"))

    for propagation in tested_propagation:
        prop = PropagationResultModel.parse_file(propagation)
        knockout_node_id = prop.network["suppressed_nodes"][0]
        prop_id = f"{knockout_node_id}_knockout"

        prop_series = propagation_diff_df(prop, reference_prop)
        if sort:
            prop_series.sort_values(inplace=True, ascending=True)

        yield prop_id, prop_series

def get_crispr_enrichments(ranking_path: str, translation_path: str):
    with open(translation_path, 'r') as handler:
        trans_dict = json.load(handler)

    df = pd.read_csv(ranking_path)[["id", "pos|rank"]]
    df.id = df.id.apply(lambda x: trans_dict.get(x, None))
    df.sort_values(by=["pos|rank"])
    return df["id"].to_numpy()


# not needed for now
def get_crispr_rankings(rankning_path, translator_path):
    with open(translator_path, 'r') as handler:
        trans_dict = json.load(handler)
    series = pd.read_csv(rankning_path, index_col="id", usecols=["id", "pos|score"])["pos|score"]
    series.rename("crispr_score", inplace=True)
    to_drop = [symbol for symbol in series.index if symbol.upper() not in trans_dict]
    series.drop(labels=to_drop, inplace=True)
    series.rename(index={symbol: trans_dict.get(symbol.upper()) for symbol in series.index}, inplace=True)
    for dup in series[series.index.duplicated()].index.to_list():
        dup_mean = series[dup].mean()
        series.drop(dup, inplace=True)
        series[dup] = dup_mean

    return series







# def perform_gsea(reference_propagation: str, tested_propagation: str, root__outdir: str

# def calc_spearman(df: pd.DataFrame):
#     return stats.spearmanr(df["propagation_diff"].to_numpy(), df["crispr_score"].to_numpy())
#
#
# def rank_prop_by_spearman(reference_propagation: str, propagations: list[str], crispr_screen_csv):
#     return {prop_id: calc_spearman(df) for prop_id, df in
#             get_scores_df(reference_propagation, propagations, crispr_screen_csv)}


# def knockout_ranking(knockout_prop: PropagationResultModel, original_prop: PropagationResultModel):
#     diff = diff_between(knockout_prop, original_prop)
#     return np.array(sorted(diff.keys(), key=lambda node: diff[node]))


def worker_main(diff_exp_genes: list[str], task_q: Queue, output_q: Queue):
    histogram = [list() for i in range(len(diff_exp_genes) + 1)]
    while True:
        task = task_q.get(block=True)
        if task == "HALT":
            break

        prop_id, prop_df = task
        negatives = prop_df < 0
        histogram[sum([1 for g in diff_exp_genes if negatives.get(g, False)])].append(prop_id)

    output_q.put(histogram)



if __name__ == "__main__":
    print("ya")
    propagation_files = Path(r"D:\propagations\new_vanilla\merged_covid").glob("*knockout.json")
    reference = r"D:\propagations\new_vanilla\merged_covid\no_knockouts.json"
    crispr_screen_csv = r"D:\mmc1.csv"
    # with open(r"D:\differential_gene_expression\https_doi.org_10_1038_s41586_020_2332_7\upregulated_genes_24h_entrez.txt",'r') as handler:
    #     up_reg_genes = [gene.strip() for gene in handler.readlines()]
    #
    # with open(r"D:\differential_gene_expression\https_doi.org_10_1038_s41586_020_2332_7\downregulated_genes_24h_entrez.txt",'r') as handler:
    #     down_reg_genes = [gene.strip() for gene in handler.readlines()]


    # diff_exp_genes = up_reg_genes + down_reg_genes

    with open(r"D:\differential_gene_expression\https_doi.org10.1038s41586-021-03493-4\diff_exp_unified_translated.json", 'r') as handler:
        diff_exp_genes = list(json.load(handler).keys())


    histogram = [list() for i in range(len(diff_exp_genes) + 1)]
    counter = 0
    task_queue = Queue()
    output_queue = Queue()
    num_workers = 14
    pool = Pool(num_workers, worker_main, [diff_exp_genes, task_queue, output_queue])
    for prop_id, prop_df in get_scores_df(reference, propagation_files):
        task_queue.put((prop_id, prop_df))
        counter += 1
        if counter % 50 == 0:
            sleep(1)

    for i in range(num_workers):
        task_queue.put("HALT")

    histogram = [list() for i in range(len(diff_exp_genes) + 1)]
    for i in range(num_workers):
        worker_histogram = output_queue.get()
        for i in range(len(histogram)):
            histogram[i].extend(worker_histogram[i])

    with open(r"D:\histogram_new_data.txt", 'w') as handler:
        handler.write("\n".join([str(i) for i in histogram]))

    plt.figure()
    plt.bar(list(range(len(histogram))), [len(i) for i in histogram])
    plt.show()