from netprop.models.results_models import PropagationResultModel
from propagation_diff import get_crispr_rankings
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
def aggregate_node_propagations(files: list[str], nodes: list[str]):
    aggregate_prop = {n: list() for n in nodes}
    for file in files:
        res = PropagationResultModel.parse_file(file)
        for n in nodes:
            try:
                aggregate_prop[n].append(res.nodes[n].liquids["info"])
            except:
                print(f"file {file} did not conribute to aggregate of {n}")
                continue
    return aggregate_prop


def calc_p_value(node: str, real_prop_value: float, randomized_props: dict):
    all_values_sorted = sorted(randomized_props[node] + [real_prop_value])
    return (1 + all_values_sorted.index(real_prop_value)) / len(all_values_sorted)


def gene_set_p_value(gene_names: list[str], real_propagation_file: str, randomized_propagations_files: list[str]):
    real_res = PropagationResultModel.parse_file(real_propagation_file)
    real_res_values = {gene: real_res.nodes[gene].liquids["info"] for gene in gene_names}
    randomized_res_values = aggregate_node_propagations(randomized_propagations_files, gene_names)
    p_values = {gene: calc_p_value(gene, real_res.nodes[gene].liquids["info"], randomized_res_values) for gene in gene_names}
    significants = {gene: p_value for gene, p_value in p_values.items() if p_value <= 0.05}
    percentage = [len([1 for g in gene_names[:i] if g in significants])/(i+1) for i in range(len(gene_names))]
    print(f"{len(significants)} significant genese found:")
    for k,v in significants.items():
        print(f"{k}: {v}")
    plt.figure()
    plt.plot(percentage)
    plt.show()

def main():
    crispr_genes = get_crispr_rankings(r"D:\mmc1.csv", r"D:\complete_translation_upper.json")
    gene_names = crispr_genes.index.to_list()[5000:6000]

    real_propagation_file = r"D:\propagations\vanilla\merged_covid\no_knockouts.json"
    randomized_propagation_files = list(Path(r"D:\propagations\no_knockouts\no_knockout_new").glob("*"))
    real_res = PropagationResultModel.parse_file(real_propagation_file)
    gene_names = list(real_res.nodes.keys())
    for gene in gene_names:
        if gene not in real_res.nodes:
            print(gene)

    to_remove = []
    for g in gene_names:
        if g not in real_res.nodes:
            to_remove.append(g)
    for g in to_remove:
        gene_names.remove(g)
    gene_set_p_value(gene_names, real_propagation_file, randomized_propagation_files)


if __name__ == "__main__":
    main()