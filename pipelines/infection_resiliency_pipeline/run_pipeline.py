
from pipeline import NonRepeatingPipeline
from pipelines.infection_resiliency_pipeline.pipeline_objs import select_subsets, resiliency_scores, score_analysis




"""
INDIVIDUAL_INTERACTORS_METADATA = "inidividual_interactors_metadata"
SUBSETS_FILE_PATH = "subsets_file_path"
PROP_RESULTS_FILE = "prop_results_file"
GENERATED_SUBSETS_DIR = "generated_subsets_dir"
GENERATED_RESILIENCY_SCORES_DIR = "generated_resiliency_scores_dir"
GENERATED_TOP_RESILIENT_GENES_PATH = "generated_top_resilient_genes_path"

RATIOS = "ratios"
MIN_SUBSET_SIZE = "min_subset_size"
MAX_SUBSET_SIZE = "max_subset_size"
TOP_PROPAGATED_RANK_THRESHOLD = "top_propagated_rank_threshold"
TOP_RESILIENCE_RANK_THRESHOLD = "top_resilience_rank_threshold"


"""

"D:\data\propagations\krogan_interactors\individual_interactors\metadata.json"

def compose_pipeline(individual_interactors_metadata: str,
                     prop_results_file: str,
                     subset_file_path: str,
                     generated_subsets_dir: str,
                     generated_resiliency_scores_dir: str,
                     generated_top_resilient_genes_path: str,
                     ratios: list[float],
                     min_subset_size: int,
                     max_subset_size: int,
                     top_prop_threshold: int,
                     top_resiliency_threshold: int,
                     descriptor: str,
                     reset_state: bool=False) -> NonRepeatingPipeline:
    print("f1")
    select_subsets_piepline = select_subsets.get_pipeline(ratios,
                                                          min_subset_size,
                                                          max_subset_size,
                                                          individual_interactors_metadata,
                                                          prop_results_file,
                                                          subset_file_path,
                                                          generated_subsets_dir,
                                                          descriptor,
                                                          reset_state=reset_state)
    print("f2")
    resiliency_scores_pipeline = resiliency_scores.get_pipeline(top_prop_threshold,
                                                                generated_resiliency_scores_dir,
                                                                descriptor,
                                                                reset_state=reset_state)
    print("f3")
    score_analysis_pipeline = score_analysis.get_pipeline(top_resiliency_threshold,
                                                          generated_top_resilient_genes_path,
                                                          descriptor,
                                                          reset_state=reset_state)

    print("f4")
    global_pipe_descriptor = descriptor + "_global_pipeline"
    global_pipeline = NonRepeatingPipeline(global_pipe_descriptor, reset_state=reset_state)
    global_pipeline.add_steps(
        [
            select_subsets_piepline,
            resiliency_scores_pipeline,
            score_analysis_pipeline
        ]
    )
    # global_pipeline.add_steps(
    #     [
    #         lambda attrs: print("l1"),
    #         select_subsets_piepline,
    #         lambda attrs: print("l2"),
    #         lambda attrs: print("l3")
    #     ]
    # )
    print("f5")
    return global_pipeline


def run_once(individual_interactors_metadata: str, prop_results_file: str, resiliency_dir: str,
             ratios: list[float], min_subset_size: int, max_subset_size: int,
             top_propagated_rank_threshold: int, top_resilient_rank_threshold: int,
             run_name: str, reset_state: bool=False):
    # resiliency_dir = "/data/resiliency"
    # individual_interactors_metadata = "/data/propagations/krogan_interactors/individual_interactors/metadata.json"
    # prop_results_file = "/data/propagations/krogan_interactors/individual_interactors/all.csv"
    # subsets_master_file_path = resiliency_dir + "/subset_files/1.json"
    # generated_subsets_dir = resiliency_dir + "/subset_files/1"
    # generated_resiliency_scores = resiliency_dir + "/scores/1"
    # generated_top_ranking = resiliency_dir + "/top_ranking/1"
    # resiliency_dir = r"/data/resiliency"
    # individual_interactors_metadata = r"/data/propagations/krogan_interactors/individual_interactors/metadata.json"
    # prop_results_file = r"/data/propagations/krogan_interactors/individual_interactors/all.csv"
    subsets_master_file_path = resiliency_dir + f"/subset_files/{run_name}.json"
    generated_subsets_dir = resiliency_dir + f"/subset_files/{run_name}"
    generated_resiliency_scores = resiliency_dir + f"/scores/{run_name}"
    generated_top_ranking = resiliency_dir + f"/top_ranking/{run_name}"
    # ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    # min_subset_size = 5
    # max_subset_size = 20
    # top_propagated_rank_threshold = 1000
    # top_resilient_rank_threshold = 100
    # run_name = "2"
    # reset_state = True

    # from utils.utils import load_json, dump_json
    # md = load_json(individual_interactors_metadata)
    # dumb_md = {"e": md["e"], "m": md["m"]}
    # dump_json(dumb_md, r"D:\data\propagations\krogan_interactors\individual_interactors\dummy_metadata_for_test.json")

    # exit()
    print("composing pipeline")
    pipeline = compose_pipeline(individual_interactors_metadata,
                                prop_results_file,
                                subsets_master_file_path,
                                generated_subsets_dir,
                                generated_resiliency_scores,
                                generated_top_ranking,
                                ratios,
                                min_subset_size,
                                max_subset_size,
                                top_propagated_rank_threshold,
                                top_resilient_rank_threshold,
                                run_name,
                                reset_state)
    #
    print("executing pipeline")
    pipeline.execute()

# if __name__ == "__main__":
#     main()
# def quack():
#     print("quack")
# main()