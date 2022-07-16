import shutil
from pathlib import Path

import pandas as pd

from netprop.models import ConfigModel
from pipeline import NonRepeatingPipeline
from pipelines.viral_analysis_pipeline.prop_helpers import create_virus_base_conf, create_randomized_conf

from pipelines.viral_analysis_pipeline.pipeline_objs.pipeline_objs_consts import *


def _prepare_root_dirs(attrs: dict):
    virus_name = attrs[VIRUS_NAME_ATTR]

    print("preparing root dirs")
    global_res_root = Path(r"D:\data\propagations")
    virus_res_root = global_res_root / (virus_name + "_individual_interactors")

    config_global_path = Path(r"D:\configurations")
    virus_config_root = config_global_path / (virus_name + "_individual_interactors")

    virus_res_root.mkdir(exist_ok=True)
    virus_config_root.mkdir(exist_ok=True)

    attrs[GLOBAL_RES_ROOT_ATTR] = global_res_root
    attrs[VIRUS_RES_ROOT_ATTR] = virus_res_root
    attrs[VIRUS_CONFIG_ROOT_ATTR] = virus_config_root


def _prepare_config_env(attrs: dict):
    print("preparing config_env")
    virus_res_root, virus_config_root = attrs[VIRUS_RES_ROOT_ATTR], attrs[VIRUS_CONFIG_ROOT_ATTR]
    conf_to_patch = attrs[CONF_TO_PATCH_ATTR]

    base_conf = ConfigModel.parse_file(conf_to_patch)
    network_conf_path = Path(base_conf.networks.path)
    random_network_conf_path = list(Path(network_conf_path.parent).glob("*covid_randomized*"))[0]
    shutil.copy(network_conf_path, virus_config_root / network_conf_path.name)
    shutil.copy(random_network_conf_path, virus_config_root / random_network_conf_path.name)

    prop_config_path = virus_config_root / "individual_interactors_conf.json"
    random_network_prop_config_path = virus_config_root / "individual_interactors_randomized_network_conf.json"
    metadata_file_path = virus_res_root / "metadata.json"

    attrs[CONF_TO_PATCH_ATTR] = conf_to_patch
    attrs[PROP_CONFIG_PATH_ATTR] = prop_config_path
    attrs[RANDOM_NETWORK_PROP_CONFIG_PATH_ATTR] = random_network_prop_config_path
    attrs[RANDOM_NETWORK_CONFIG_PATH] = random_network_conf_path
    attrs[METADATA_FILE_PATH_ATTR] = metadata_file_path


def _create_metadata_file(attrs: dict):
    prior_set_source, index_col, prior_set_col = attrs[PRIOR_SET_SOURCE_ATTR], attrs[INDEX_COL_ATTR], attrs[PRIOR_SET_COL_ATTR]

    interactors_series = pd.read_csv(prior_set_source, index_col=index_col)[prior_set_col].dropna().astype(int)
    dict_of_lists = dict()
    for item in interactors_series.iteritems():
        viral_protein, interactor = item[0], item[1]
        if viral_protein not in dict_of_lists:
            dict_of_lists[viral_protein] = []
        dict_of_lists[viral_protein].append(str(interactor))


def _copy_and_patch_conf(attrs: dict):
    print("copying and patching conf")
    conf_to_patch = attrs[CONF_TO_PATCH_ATTR]
    metadata_file_path = attrs[METADATA_FILE_PATH_ATTR]
    virus_res_root = attrs[VIRUS_RES_ROOT_ATTR]
    virus_name = attrs[VIRUS_NAME_ATTR]
    prop_config_path = attrs[PROP_CONFIG_PATH_ATTR]
    # virus_config_root = attrs["virus_config_root"]
    random_network_prop_config_path = attrs[RANDOM_NETWORK_PROP_CONFIG_PATH_ATTR]
    random_network_conf_path = attrs[RANDOM_NETWORK_CONFIG_PATH]

    create_virus_base_conf(conf_to_patch, str(metadata_file_path),
                           str(virus_res_root), virus_name + "_individual_interactors",
                           str(prop_config_path))
    create_randomized_conf(str(prop_config_path), random_network_conf_path, str(random_network_prop_config_path)) # irus_config_root / random_network_conf_path.name


def get_pipeline(virus_name: str, conf_to_patch: str, prior_set_source: str, reset_states: bool=False) -> NonRepeatingPipeline:
    init_attrs = {
        VIRUS_NAME_ATTR: virus_name,
        CONF_TO_PATCH_ATTR: conf_to_patch,
        PRIOR_SET_SOURCE_ATTR: prior_set_source
    }
    prop_prep_pipeline = NonRepeatingPipeline(init_attrs=init_attrs, state_suffix=virus_name + "_prop_prep_pipeline_state",
                                              reset_state=reset_states)
    prop_prep_pipeline.add_steps(
        [
            _prepare_root_dirs,
            _prepare_config_env,
            _create_metadata_file,
            _copy_and_patch_conf
        ]
    )

    return prop_prep_pipeline
