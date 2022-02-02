from score_from_netprop import propagation_to_df
from gsea_analysis import gsea_result_spearman
from gsea.gsea_analysis import gene_set_enrichment
from netprop.models import PropagationResultModel
from load_diff_exp_genes import load
from crispr import get_crispr_rankings


def single_prop_analysis(prop_file: str, diff_exp_file: str):
    prop = PropagationResultModel.parse_file(prop_file)
    prop_id, prop_series = propagation_to_df(prop, sort=True)
    diff_exp = load(diff_exp_file, 'stukalov')
    gsea_scores, leading_edge = gene_set_enrichment(prop_series, diff_exp)
    print(f"score: {gsea_scores}, leading edge size: {len(leading_edge)}")


def multi_prop_analysis(prop_files: list[str], diff_exp_file: str, id_of_interest: str,
                        merge_liquids: bool=False, by_liquid: str="info"):
    diff_exp = load(diff_exp_file, 'stukalov')
    gsea_scores = {}
    for prop_file in prop_files:
        prop = PropagationResultModel.parse_file(prop_file)
        try:
            prop_id, prop_series = propagation_to_df(prop, sort=True, by_liquid=by_liquid)
        except:
            print(f"problem with {prop_file}")
            continue
        gsea_scores[prop_id] = gene_set_enrichment(prop_series, diff_exp)
    sorted_props = sorted(list(gsea_scores.keys()), key=lambda k: gsea_scores[k][0])
    print(f"sorted {len(sorted_props)} props, and {id_of_interest} is at index {sorted_props.index(id_of_interest)}")
    for k in sorted_props:
        print(f"{k}: {gsea_scores[k][0]}")
    print(1 + sorted_props.index(id_of_interest) / len(sorted_props))