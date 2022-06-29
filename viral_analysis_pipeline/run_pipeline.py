from pathlib import Path
from netprop.models import ConfigModel
import shutil
from functools import partial

from viral_analysis_pipeline.prop import create_virus_base_conf, create_metadata_file,\
                                         prop_from_conf_and_save_new_format, create_randomized_conf

from viral_analysis_pipeline.splits_data import generate_intersection_data
from pipelines import NonRepeatingPipeline

def prepare_root_dirs(virus_name: str, attrs: dict):
    print("preparing root dirs")
    global_res_root = Path(r"D:\data\propagations")
    virus_res_root = global_res_root / (virus_name + "_individual_interactors")

    config_global_path = Path(r"D:\configurations")
    virus_config_root = config_global_path / (virus_name + "_individual_interactors")

    virus_res_root.mkdir(exist_ok=True)
    virus_config_root.mkdir(exist_ok=True)

    attrs["virus_name"] = virus_name
    attrs["global_res_root"] = global_res_root
    attrs["virus_res_root"] = virus_res_root
    attrs["virus_config_root"] = virus_config_root


def prepare_config_env(conf_to_patch: str, prior_set_source: str, attrs: dict):
    print("preparing config_env")
    virus_res_root = attrs["virus_res_root"]

    virus_config_root = attrs["virus_config_root"]
    base_conf = ConfigModel.parse_file(conf_to_patch)
    network_conf_path = Path(base_conf.networks.path)
    random_network_conf_path = list(Path(network_conf_path.parent).glob("*covid_randomized*"))[0]
    shutil.copy(network_conf_path, virus_config_root / network_conf_path.name)
    shutil.copy(random_network_conf_path, virus_config_root / random_network_conf_path.name)

    prop_config_path = virus_config_root / "individual_interactors_conf.json"
    random_network_prop_config_path = virus_config_root / "individual_interactors_randomized_network_conf.json"
    metadata_file_path = virus_res_root / "metadata.json"

    create_metadata_file(prior_set_source, str(metadata_file_path))

    attrs["conf_to_patch"] = conf_to_patch
    attrs["prop_config_path"] = prop_config_path
    attrs["random_network_prop_config_path"] = random_network_prop_config_path
    attrs["random_network_config_path"] = random_network_conf_path
    attrs["metadata_file_path"] = metadata_file_path


def copy_and_patch_conf(attrs: dict):
    print("copying and patching conf")
    conf_to_patch = attrs["conf_to_patch"]
    metadata_file_path = attrs["metadata_file_path"]
    virus_res_root = attrs["virus_res_root"]
    virus_name = attrs["virus_name"]
    prop_config_path = attrs["prop_config_path"]
    # virus_config_root = attrs["virus_config_root"]
    random_network_prop_config_path = attrs["random_network_prop_config_path"]
    random_network_conf_path = attrs["random_network_config_path"]


    create_virus_base_conf(conf_to_patch, str(metadata_file_path),
                           str(virus_res_root), virus_name + "_individual_interactors",
                           str(prop_config_path))
    create_randomized_conf(str(prop_config_path), random_network_conf_path, str(random_network_prop_config_path)) # irus_config_root / random_network_conf_path.name


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
    generate_intersection_data(str(metadata_file_path), rand_prop_res_path, str(Path(virus_res_root) / "intersections_data"),
                               10, [20, 100])


def run_virus_pipeline(virus_name: str, conf_to_patch: str, prior_set_source: str, reset_states: bool=False):
    prepare_virus_root_dir = partial(prepare_root_dirs, virus_name)
    prepare_virus_config_env = partial(prepare_config_env, conf_to_patch, prior_set_source)

    conf_pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_conf_pipeline_state", reset_state=reset_states)
    conf_pipeline.add_steps(
        [
            prepare_virus_config_env,
            copy_and_patch_conf
        ])

    analysis_pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_analysis_pipeline_state", reset_state=reset_states)
    analysis_pipeline.add_steps(
        [
            intersection_data_from_prop_results
        ])

    virus_pipeline = NonRepeatingPipeline(state_suffix=virus_name + "_main_pipeline_state", reset_state=reset_states)
    virus_pipeline.add_steps(
        [
            prepare_virus_root_dir,
            conf_pipeline,
            create_prop_from_config,
            analysis_pipeline
        ]
    )

    virus_pipeline.execute()

if __name__ == "__main__":
    run_virus_pipeline("sars_cov_1",
                       r"D:\configurations\krogan_invidual_interactors\krogan_individual_interactors_conf.json",
                       r"D:\data\networks\cov_to_human\gordon_sars_cov_1_high_confidence_translated.csv", reset_states=True)