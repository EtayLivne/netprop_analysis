from utils.new_file_loaders import NetpropCovToHumanLoader
from  score_from_netprop import *
from crispr import get_crispr_rankings
from scipy.stats import spearmanr
from pathlib import Path

if __name__ == "__main__":
    # network = NetpropCovToHumanLoader(r"D:\data\networks\try.json",
    #                                   r"D:\data\networks\metastudy\protein_interactome_translated.csv",
    #                                   r"D:\data\networks\metastudy\rna_interactome_translated.csv",
    #                                   merged_covid=False).load()
    # NetpropCovToHumanLoader.record_network(network, r"D:\data\networks\try_with_covid.json")
    #df = propagations_to_df([str(e) for e in list(Path(r"D:\data\propagations\randomized_new_network_loader").glob("*"))])
    df = pd.read_csv("temp.csv")
    print("we have a df!")
    # df.to_csv("temp.csv")
    p_values = infer_p_value(df, "ref")

    crispr_scores = get_crispr_rankings(r"D:\data\other\crispr_screen.csv", r"D:\data\other\symbol_to_entrezgene.json")
    crispr_df = crispr_scores.to_frame()
    crispr_df.reset_index(inplace=True)
    crispr_df.rename(inplace=True, columns={"id": "nodes"})
    crispr_df.sort_values(by=["crispr_score"], inplace=True, ascending=True)
    top_crispr = set(crispr_df.head(25)["nodes"])
    #sorted_df = df.sort_values(by=["ref"], ascending=False)
    p_values.sort_values(by=["p_values"], inplace=True, ascending=True)

    #for i in range(50):
    #    n1 = df.iloc[i]["nodes"]
    #    ref1 = df.iloc[i]["ref"]
    #    n2 = sorted_df.iloc[i]["nodes"]
    #    ref2 = sorted_df.iloc[i]["ref"]
    #    print(f"{i}: {n1}: {ref1}        {n2}: {ref2}")


    i = 0
    for index, row in p_values.iterrows():
        node = row["nodes"]
        if node in top_crispr:
            print(f"{node}: {i}")
        i += 1

    # final_df = pd.merge(p_values, crispr_df, how="inner", on="nodes")
    #
    # print(spearmanr(final_df["p_values"], final_df["crispr_score"]))
