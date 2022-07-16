from pathlib import Path
from pipelines.viral_analysis_pipeline.pipeline_objs import prop, intersect, analyze, prop_prep
from pipelines.viral_analysis_pipeline.prop_helpers import prop_from_conf_and_save_new_format

from pipelines.viral_analysis_pipeline.splits_data import generate_intersection_data
from pipeline import NonRepeatingPipeline




def create_prop_from_config(attrs: dict):
    print("creating prop from config")
    prop_config_path = attrs["prop_config_path"]
    random_network_prop_config_path = attrs["random_network_prop_config_path"]

    attrs["real_prop_res_path"] = prop_from_conf_and_save_new_format(str(prop_config_path))
    attrs["rand_prop_res_path"] = prop_from_conf_and_save_new_format(str(random_network_prop_config_path))


def intersection_data_from_prop_results(attrs: dict):
    print("analyzing intersection data from prop results")

    metadata_file_path = attrs["metadata_file_path"]
    virus_res_root = attrs["virus_res_root"]
    real_prop_res_path = attrs["real_prop_res_path"]
    rand_prop_res_path = attrs["rand_prop_res_path"]

    generate_intersection_data(str(metadata_file_path), real_prop_res_path, str(Path(virus_res_root) / "intersections_data"), 10, [20, 100])
    # generate_intersection_data(str(metadata_file_path), rand_prop_res_path, str(Path(virus_res_root) / "intersections_data"),
    #                            10, [20, 100])


def run_virus_pipeline(virus_name: str, conf_to_patch: str, prior_set_source: str, reset_states: bool=False):
    pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_backbone_pipeline_state", reset_state=reset_states)
    pipeline.add_steps(
        [
            prop_prep.get_pipeline(virus_name, conf_to_patch, prior_set_source, reset_states=reset_states),
            prop.get_pipeline(virus_name, reset_state=reset_states),
            intersect.get_pipeline(virus_name, reset_state=reset_states),
            analyze.get_pipeline(virus_name, reset_state=reset_states)
        ]
    )

    pipeline.execute()

    # prepare_virus_root_dir = partial(prepare_root_dirs, virus_name)
    # prepare_virus_config_env = partial(prepare_config_env, conf_to_patch, prior_set_source)
    #
    # conf_pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_conf_pipeline_state", reset_state=reset_states)
    # conf_pipeline.add_steps(
    #     [
    #         prepare_virus_config_env,
    #         copy_and_patch_conf
    #     ])
    #
    # analysis_pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_analysis_pipeline_state", reset_state=reset_states)
    # analysis_pipeline.add_steps(
    #     [
    #         intersection_data_from_prop_results
    #     ])
    #
    # virus_pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_main_pipeline_state", reset_state=reset_states)
    # virus_pipeline.add_steps(
    #     [
    #         prepare_virus_root_dir,
    #         conf_pipeline,
    #         create_prop_from_config,
    #         analysis_pipeline
    #     ]
    # )
    #
    # virus_pipeline.execute()

if __name__ == "__main__":
    # run_virus_pipeline("hep_c",
    #                    r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_conf.json",
    #                    r"D:\data\networks\virus_to_human\hep_c_host_protein_interactome_translated_without_nonexsiting.csv", reset_states=True)

    run_virus_pipeline("influenza_pra1av",
                       r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_conf.json",
                       r"D:\data\networks\virus_to_human\pr81av_host_protein_interactome_translated.csv", reset_states=True)