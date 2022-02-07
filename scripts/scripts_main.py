import pandas as pd
from scripts.correlate_prop_scores_to_diff_exp import single_prop_analysis, multi_prop_analysis
from pathlib import Path
from scripts.split_csv import split_csv, to_entrez, extract_uniprot, ncbi_query_strings, back_to_entrez, update_translation_dict
from utils.new_file_loaders import CovToHumanMeta
from scripts.prop_config_from_network import patch_prior_set_by_network_tag, transfer_suppressed_sets

if __name__ == "__main__":
    # multi_prop_analysis(list(Path(r"D:\data\propagations\randomized\merged_covid").glob("*.json")) + [r"D:\data\propagations\new_vanilla\merged_covid\no_knockouts.json"],
    #                     r"D:\data\diff_genes\stukalov\diff_exp_unified_translated.json",
    #                     "no_knockouts")


    # root_dir = Path(r"D:\data\networks\metastudy")
    # new_dfs = split_csv(r"D:\data\networks\metastudy\all_interactions_raw.csv", "Assay type")
    # for df_name in list(new_dfs.index):
    #     new_dfs[df_name].to_csv(root_dir / f"{df_name}.csv", index=False)

    # to_entrez(r"D:\data\networks\metastudy\RNA interactome.csv",
    #           r"D:\data\networks\metastudy\rna_interactome_translated.csv",
    #           r"D:\data\other\symbol_to_entrezgene.json")
    # extract_uniprot(r"D:\data\networks\metastudy\rna interactome_translated.csv",
    #                 r"D:\data\networks\metastudy\rna_uniprots.txt")
    # ncbi_query_strings(r"D:\data\networks\metastudy\proteome_uniprots.txt",
    #                    r"D:\data\networks\metastudy\proteome_uniprots_query_strings.txt")
    # nw = CoVToHumanMeta(R"D:\data\networks\H_sapiens_aug_2020.net",
    #                     protein_interactions_path=R"D:\data\networks\metastudy\protein_interactome_translated.csv",
    #                     rna_interactions_path=r"D:\data\networks\metastudy\rna_interactome_translated.csv").load()
    #
    # patch_prior_set_by_network_tag(nw, {"species_id": "sars-cov-2"}, r"D:\configurations\test_config.json",
    #                                r"D:\configurations\new_test_config.json")
    transfer_suppressed_sets(r"D:\configurations\rna_prop_conf.json", r"D:\configurations\all_knockouts_new.json",
                             r"D:\configurations\rna_prop_knockouts_conf.json")

    # new_translations_proteins = back_to_entrez(r"D:\data\networks\metastudy\proteome_uniprots_query_strings naama.txt", None,
    #                                       r"D:\data\networks\metastudy\proteome_outputs.txt")
    #
    # new_translations_rna = back_to_entrez(r"D:\data\networks\metastudy\rna_uniprots_query_strings_naama.txt"
    #                                       ,None,r"D:\data\networks\metastudy\rna_outputs.txt")
    #
    # new_tranlations = new_translations_proteins | new_translations_rna
    # update_translation_dict(r"D:\data\other\symbol_to_entrezgene.json",new_tranlations,r"D:\data\other\symbol_to_entrezgene_new.json" )
    # x = 7

    # multi_prop_analysis(list(Path(r"D:\data\propagations\randomized\merged_covid").glob("*.json")) + [r"D:\data\propagations\new_vanilla\merged_covid\no_knockouts.json"],
    #                     r"D:\data\diff_genes\stukalov\diff_exp_unified_translated.json",
    #                     "no_knockouts")