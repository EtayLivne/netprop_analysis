from netprop.models import PropagationResultModel
from crispr import get_crispr_rankings
import pandas as pd
from scipy.stats import spearmanr

def correlate(df: pd.DataFrame, col1, col2):
    df = df.dropna()
    return spearmanr(df[col1], df[col2])

def correlate_prop_to_crispr(crispr_path: str, prop_path: str, by_liquid="info"):
    prop_df = PropagationResultModel.parse_file(prop_path).prop_scroes_as_series()
    crispr_df = get_crispr_rankings()
    crispr_df.reset_index()
    crispr_df.rename({"index": "nodes"})
    merged = pd.merge([prop_df, crispr_df], on="nodes", how="inner")
    relevant_props = [c for c in merged.columns if c.split(".")[1] == by_liquid]

    return {c: correlate(merged, c, "crispr_score") for c in relevant_props}