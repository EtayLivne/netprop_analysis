from prepare_for_gsea import *
from crispr import get_crispr_rankings
import pandas as pd
from scipy import stats
from gsea.gsea_data import GSEAData
from gsea.gsea_analysis import GSEAAnalysis

from utils.utils import load_json, dump_json

def crispr_scores_for_spearman(cripsr_scores_path: str, translation_path: str)-> pd.DataFrame:
    crispr_ranked_nodes = get_crispr_rankings(cripsr_scores_path, translation_path)
    return 1/crispr_ranked_nodes



def gsea_result_spearman(scored_knockouts_path: str, crispr_scores_path: str, translation_path: str):
    # scored_knockouts_path = r"D:\analysis_output.json"
    # crispr_scores_path = r"D:\data\other\crispr_screen.csv"
    # translation_path = r"D:\data\other\symbol_to_entrezgene.json"
    with open(scored_knockouts_path, 'r') as handler:
        gsea_data = json.load(handler)
    gsea_data = {k.split("_")[0]: v[0] for k, v in gsea_data.items()}
    crispr_ranked_nodes = crispr_scores_for_spearman(crispr_scores_path, translation_path)
    gsea_input = []
    crispr_input = []
    exclude_counter = 0
    for gene in gsea_data:
        if gene not in crispr_ranked_nodes.iloc[:18671].index:
            print(f"excluding {gene}")
            exclude_counter += 1
            continue
        gsea_input.append(gsea_data[gene])
        crispr_input.append(crispr_ranked_nodes[gene])
    # df = pd.merge(pd.DataFrame.from_dict(gsea_data, orient="index", columns="genes"), crispr_ranked_nodes.to_frame(), how="inner")
    print(f"{exclude_counter} genes have been excluded")
    return stats.spearmanr(gsea_input, crispr_input)


def gsea_scores(ref_prop_file: str, histogram_file: str, props_root_folder: str, diff_exp_genes: str, outpath: str):
    data = GSEAData()
    data.reference_propagation =ref_prop_file
    data.propagation_files_from_histogram(histogram_file, props_root_folder)
    # data.set_target_pathway_to_diff_exp_genes(diff_exp_genes)
    data.set_target_pathway_to_top_prop_scores(ref_prop_file)
    analysis = GSEAAnalysis()
    analysis.analyze(data, outpath)

if __name__ == "__main__":
    # print(gsea_result_spearman())


    pass









# def single_sample_gsea_analysis(gct_file: str, gmt_file: str, outdir):
#     return gp.ssgsea(gct_file, gmt_file, outdir=outdir, sample_norm_method='custom', no_plot=True, permutation_num=1000)