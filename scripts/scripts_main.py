from utils.new_file_loaders import NetpropCovToHumanLoader
from score_from_netprop import *
from crispr import get_crispr_rankings
from scipy.stats import spearmanr
from pathlib import Path
from netprop.networks.loaders import NetpropNetwork, CombinedHumanCovidNetworkLoader
from netprop.models import ConfigModel
from netprop.propagation import propagate_from_config
import networkx as nx
from random import shuffle, sample
from math import ceil
import json
from metrics.roc_auc import *
from crispr import get_huh_crispr
from scripts.divide_and_conquer import split_to_randomized_component, merge_results_as_different_liquids, get_nodes_in_intersection, intersect_crispr_ko_with_node_intersections



if __name__ == "__main__":
    # g = nx.cycle_graph(24)
    # np.random.seed(100)
    # g = nx.connected_watts_strogatz_graph(20, 5, 0.5, seed=100)
    # labels = {}
    # scores = np.random.randint(0, 255, 20)
    # node_scores = {f"{i}": scores[i] for i in range(20)}
    # draw_network(g, labels, node_scores)
    # df = pd.DataFrame({"node": ["CHD6", "MBPTS2", "ACE2", "VAC14", "EXOC2"], "score": [0.86, 0.85, 0.63, 0.57, 0.52]})
    # print(df)

    # from netprop.networks.loaders import CombinedHumanCovidNetworkLoader
    # krogan_nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net", r"D:\data\networks\krogan_cov_to_human.cx", r"D:\data\other\symbol_to_entrezgene.json", merge_covid=False).load()
    # #huh_7_nw = NetpropNetwork(r"D:\data\networks\metastudy\complete_networks\huh7_only_fixed.json").load()
    # base_conf = r"C:\studies\thesis\code\cicd\data_cache\metastudy_conf.json"
    # new_conf_separated = r"C:\studies\thesis\code\cicd\data_cache\krogan_separated_all_interactors_conf.json"
    # new_conf_merged = r"C:\studies\thesis\code\cicd\data_cache\krogan_merged_all_interactors_conf.json"
    # prop_output_dir_separated = "/output/props/krogan_all_interactors_seperated"
    # prop_output_dir_merged = "/output/props/krogan_all_interactors_merged"
    # cov_interactors_prior_set(krogan_nw, base_conf, new_conf_separated, prop_output_dir=prop_output_dir_separated)
    # cov_interactors_prior_set(krogan_nw, base_conf, new_conf_merged, prop_output_dir=prop_output_dir_merged, merge_sources=True)
    # merge_sources: bool = False

    # df = pd.read_csv(r"D:\data\networks\metastudy\protein_interactome_translated.csv")
    # df["Bait"] = df["Bait"].apply(lambda name: name.upper())
    # df.to_csv(r"D:\data\networks\metastudy\protein_interactome_translated_upper_case.csv")

    #
    # m = NormsROC(prop_file_path=r"D:\data\propagations\krogan_interactors\specific_proteins\all_together.json",
    #              by_liquids=None,
    #              hit_set_load_method=get_huh_crispr)
    # # m = SinglePropROC(prop_file_path=r"D:\data\propagations\krogan_interactors\merged\all.json",
    # #                   by_liquids=None,
    # #                   hit_set_load_method=get_huh_crispr)
    # m.load_hit_set(r"D:\data\other\huh7_crispr_translated.csv")
    # res = m.calc_metric()
    # # res.to_csv("temp.json")
    # m.show_results(res)

    # merge_results_as_different_liquids(list(Path(r"D:\data\propagations\krogan_interactors\specific_proteins").glob("*")), r"D:\data\propagations\krogan_interactors\specific_proteins\all_together")

    # multi_prior_set_conf = r"C:\studies\thesis\code\cicd\data_cache\krogan_separated_all_interactors_conf.json"
    # output_dir = r"D:\configurations\krogan_split_interactors"
    # max_splits_per_component = 20
    # split_to_randomized_component(multi_prior_set_conf, output_dir, max_splits_per_component, override_conf_out_path=r"D:\data\propagations\metastudy_split_interactors")

    # files = Path(r"D:\configurations\krogan_split_interactors").glob("*.json")
    # for f in files:
    #     protein = str(f.name).split("_")[0]
    #     with open(str(f), 'r') as handler:
    #         d = json.load(handler)
    #     d.pop("norm")
    #     d.pop("norm_kwargs")
    #
    #     with open(str(f), "w") as handler:
    #         json.dump(d, handler, indent=4)

    # glob_list = [Path(f).glob("*.json") for f in Path(r"D:\data\propagations\krogan_split_interactors").glob("*")]
    # intersection_dicts = []
    # for folder_files in glob_list:
    #     print(f"now handling folder {len(intersection_dicts)}")
    #     intersection_dicts.append(get_nodes_in_intersection([str(f) for f in folder_files], k=200))
    #
    # with open("intersection_dicts_top_200", "w") as handler:
    #     json.dump(intersection_dicts, handler, indent=4)

    with open("intersection_dicts_top_50", "r") as handler:
        intersection_dicts = json.load(handler)
    print(intersect_crispr_ko_with_node_intersections(r"D:\data\other\huh7_crispr_translated.csv", intersection_dicts))







































    #roc = HuhCrisprROC("prop_file_path", "hits_file_path")
    #roc_df = roc.calc_metric()
    #partial_cov_interactors_prior_set()
    # cov_interactors_prior_set(method="iterative")
    # propagate_stefan()
    # compare_results(r"D:\data\propagations\stefan_cross_validation\self test\cov_interactors_analytical.json",
    #                 "D:\data\propagations\stefan_cross_validation\self test\cov_interactors_iterative.json")

    # nw = CombinedHumanCovidNetworkLoader(r"D:\data\networks\H_sapiens_aug_2020.net",
    #                                      r"D:\data\networks\krogan_cov_to_human.cx",
    #                                      r"D:\data\other\symbol_to_entrezgene.json",
    #                                      merge_covid=False).load()
    #
    # nw.remove_nodes_from([n for n, data in nw.nodes(data=True) if data["species_id"] == "sars-cov-2"])
    # NetpropNetwork.record_network(nw, r"D:\data\networks\no_covid.json")
    #
    #
    # with open(r"D:\configurations\stefan1\stefan_conf_backup.json", 'r') as handler:
    #     d = json.load(handler)
    # for prior_set in d["prior_set"]:
    #     for node in prior_set["nodes"]:
    #         node["source_of"] = ["info"]
    # with open(r"D:\configurations\stefan1\stefan_conf.json", 'w') as handler:
    #     json.dump(d, handler, indent=4)
    #propagate_stefan()

    # for conf_score in [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
    #     print(f"***  confidence {conf_score}  ***")
    #     partial_cov_interactors_prior_set(confidence=conf_score)
    #     propagate_stefan()
    #     files = [str(p) for p in Path(r"D:\data\propagations\stefan_cross_validation").glob("*.json")]
    #     intersect_top_propagated(files)
    #     print("\n\n\n")

    # for conf_score in [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]:
    #     print(f"***  confidence {conf_score}  ***")
    #     individual_cov_interactors_prior_set(confidence=conf_score)
    #     propagate_stefan()
    #     tops = top_prop_by_source(r"D:\data\propagations\stefan_unique_liquids\covid_interactors.json")
    #     top_intersected = most_intersected(tops)
    #     print(f"total liquids: {len(tops)}")
    #     histogram = {}
    #     for node, intersection_list in top_intersected.items():
    #         num = len(intersection_list)
    #         if num not in histogram:
    #             histogram[num] = 0
    #         histogram[len(intersection_list)] += 1
    #
    #     for num in sorted(list(histogram.keys())):
    #         print(f"{num}: {histogram[num]}")
    #
    #     print("\n\n\n")

    # specific_protein_partial_cov_interactors_prior_set()
    # stefan_propagations([str(file) for file in Path(r"D:\configurations\stefan1\specific_proteins").glob("*.json")])
    #sizes = set()
    #for file in [str(f) for f in Path(f"D:\configurations\stefan1\specific_proteins").glob("*.json")]:
        #with open(file, 'r') as handler:
        #    d = json.load(handler)
        #size = len(d["prior_set"][0]["nodes"]) + len(d["prior_set"][1]["nodes"])
        #print(f"{file}: {size}")
        #sizes.add(size)

    #nsp4: 8, nsp7: 32, nsp8: 24
    # for size in [9, 11, 15, 20, 26, 30, 32, 40, 47]:
    #     random_subgroups_partial_cov_interactors_prior_set(size, 100)

    #
    # files = []
    # processes = []
    # for triplet in [(40, 47)]:
    #     for folder in [Path(f"D:\\configurations\\stefan1\\random_subgroups\{num}") for num in triplet]:
    #         processes.append(Process(target=stefan_propagations, args=([str(f) for f in Path(folder).glob("*.json")],)))
    #     #     for file in [str(f) for f in Path(folder).glob("*.json")]:
    #     #         files.append(file)
    #     # stefan_propagations(files)
    #     for p in processes:
    #         p.start()
    #     for p in processes:
    #         p.join()



    # for folder in Path(r"D:\data\propagations\stefan_cross_validation\individual_proteins").glob("*"):
    #     files = [str(p) for p in folder.glob("*chunk*")]
    #     print(str(folder))
    #     more_detailed_intersection_data(files, output_file=str(folder / "detailed_intersection_data.json"))


    # for folder in Path(r"D:\data\propagations\stefan_cross_validation\random_subgroups").glob("*"):
    #     for subfolder in folder.glob("*"):
    #         files = [str(p) for p in subfolder.glob("*chunk*")]
    #         print(str(subfolder))
    #         more_detailed_intersection_data(files, output_file=str(subfolder / "detailed_intersection_data.json"))

    # p_value_dict = calc_protein_p_value(r"D:\configurations\stefan1\specific_proteins",
    #                                     r"D:\data\propagations\stefan_cross_validation\individual_proteins")
    # for k, v in p_value_dict.items():
    #     print(f"{k}: {v}")


    # intersection_dict = get_intersection_quality_metrics(r"D:\data\propagations\stefan_cross_validation\individual_proteins")
    # intersection_bar_plot(intersection_dict)

    # folders = [r"D:\data\propagations\stefan_cross_validation\random_subgroups\nsp4",
    #            r"D:\data\propagations\stefan_cross_validation\nsp7",
    #            r"D:\data\propagations\stefan_cross_validation\nsp8"]
    # for folder in folders:
    #     protein = folder.split("\\")[-1]
    #     files = [str(p) for p in Path(folder).glob("*.json")]
    #
    #     intersect_top_propagated(files,
    #                              output_file=f"D:\\data\\propagations\\stefan_cross_validation\\analysis_{protein}.json")
    # root_folders = [folder_path for folder_path in Path(r"D:\data\propagations\stefan_cross_validation\random_subgroups").glob("*") if folder_path.is_dir()]
    # folder_lists = [
    #     [str(subfolder) for subfolder in folder.glob("*")] for folder in root_folders
    # ]
    # for i in range(len(folder_lists)):
    #     print(f"{root_folders[i]}: {average_intersection(folder_lists[i])}")
    # average_intersection()
    # for folder in Path(r"D:\data\propagations\stefan_cross_validation\random_subgroups").glob("*"):
    #     files = [str(p) for p in folder.glob("*.json")]
    #     intersect_top_propagated(files, output_file=r"D:\data\propagations\stefan_cross_validation\random_subgroups\analysis8_16_24.json")

    # from utils.new_file_loaders import CovToHumanMeta
    # fixed_huh7_nw = CovToHumanMeta(r"C:\studies\thesis\code\cicd\data_cache\H_sapiens_aug_2020.net",
    #                                r"D:\data\networks\metastudy\protein_interactome_translated.csv",
    #                                r"D:\data\networks\metastudy\rna_interactome_translated.csv",
    #                                rna_cell_lines=["HUH7", "HUH7.5"], protein_cell_lines="all", min_num_studies=2).load()
    # NetpropNetwork.record_network(fixed_huh7_nw,  r"D:\data\networks\metastudy\huh7_only_fixed.json")



