from utils.new_file_loaders import NetpropCovToHumanLoader
from  score_from_netprop import *
from crispr import get_crispr_rankings
from scipy.stats import spearmanr

if __name__ == "__main__":
    # network = NetpropCovToHumanLoader(r"D:\data\networks\try.json",
    #                                   r"D:\data\networks\metastudy\protein_interactome_translated.csv",
    #                                   r"D:\data\networks\metastudy\rna_interactome_translated.csv",
    #                                   merged_covid=False).load()
    # NetpropCovToHumanLoader.record_network(network, r"D:\data\networks\try_with_covid.json")
    df = propagations_to_df([r"D:\data\propagations\new_network_loader\combined_human_merged_covid.json"])
    crispr_scores = get_crispr_rankings(r"D:\data\other\crispr_screen.csv", r"D:\data\other\symbol_to_entrezgene.json")
    crispr_df = crispr_scores.to_frame()
    crispr_df.reset_index(inplace=True)
    crispr_df.rename(inplace=True, columns={"id": "nodes"})
    final_df = pd.merge(df, crispr_df, how="inner", on="nodes")
    print(spearmanr(final_df["combined_human_merged_covid.info"], final_df["crispr_score"]))
