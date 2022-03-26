from utils.new_file_loaders import NetpropCovToHumanLoader
from  score_from_netprop import *
from crispr import get_crispr_rankings
from scipy.stats import spearmanr
from pathlib import Path
from metrics.roc_auc import HuhCrisprROC


if __name__ == "__main__":
    roc = HuhCrisprROC("prop_file_path", "hits_file_path")
    roc_df = roc.calc_metric()
