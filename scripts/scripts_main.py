from utils.new_file_loaders import NetpropCovToHumanLoader
from  score_from_netprop import *
if __name__ == "__main__":
    # network = NetpropCovToHumanLoader(r"D:\data\networks\try.json",
    #                                   r"D:\data\networks\metastudy\protein_interactome_translated.csv",
    #                                   r"D:\data\networks\metastudy\rna_interactome_translated.csv",
    #                                   merged_covid=False).load()
    # NetpropCovToHumanLoader.record_network(network, r"D:\data\networks\try_with_covid.json")
    df = propagations_to_df([r"D:\data\propagations\new_network_loader\combined_human_merged_covid.json"])
    y = 4