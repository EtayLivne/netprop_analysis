from pipeline import NonRepeatingPipeline
from pipelines.viral_analysis_pipeline.prop_helpers import prop_from_conf_and_save_new_format
from pipelines.viral_analysis_pipeline.pipeline_objs.pipeline_objs_consts import *

def _create_prop_from_config(attrs: dict):
    print("creating prop from config")
    prop_config_path = attrs[PROP_CONFIG_PATH_ATTR]
    random_network_prop_config_path = attrs[RANDOM_NETWORK_PROP_CONFIG_PATH_ATTR]

    attrs[REAL_PROP_RES_PATH_ATTR] = prop_from_conf_and_save_new_format(str(prop_config_path))
    attrs[RAND_PROP_RES_PATH_ATTR] = prop_from_conf_and_save_new_format(str(random_network_prop_config_path))


def get_pipeline(virus_name: str, reset_state: bool=False) -> NonRepeatingPipeline:
    pipe = NonRepeatingPipeline(state_suffix=virus_name + "_prop_prep_pipeline_state", reset_state=reset_state)
    pipe.add_step(_create_prop_from_config)
    return pipe
