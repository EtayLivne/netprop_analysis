import pandas as pd

from utils.utils import load_json
from pipeline import NonRepeatingPipeline
from pipelines.infection_resiliency_pipeline.attrs_keys import *
from pipelines.infection_resiliency_pipeline.select_subsets_helpers import subsets_file, rank_by_subset_df


def _gen_subsets_file(attrs: dict) -> None:
    individual_interactors_metadata = attrs[INDIVIDUAL_INTERACTORS_METADATA]
    subsets_file_path = attrs[SUBSETS_FILE_PATH]
    ratios = attrs[RATIOS]
    min_subset_size = attrs[MIN_SUBSET_SIZE]
    max_subset_size = attrs[MAX_SUBSET_SIZE]

    sets = load_json(individual_interactors_metadata)
    subsets_file(sets, ratios, min_subset_size, max_subset_size, subsets_file_path)


def _gen_ranked_dfs_by_subset(attrs: dict) -> None:
    subsets_file_path = attrs[SUBSETS_FILE_PATH]
    prop_results_file = attrs[PROP_RESULTS_FILE]
    generated_subsets_dir = attrs[GENERATED_SUBSETS_DIR]

    rank_by_subset_df(subsets_file_path, prop_results_file, generated_subsets_dir)


def get_pipeline(ratios: list[float],
                 min_subset_size: int,
                 max_subset_size: int,
                 individual_interactors_metadata_path: str,
                 prop_results_path: str,
                 path_to_save_subsets_file: str,
                 path_to_save_subsets_dir: str,
                 descriptor: str,
                 reset_state: bool=False) -> NonRepeatingPipeline:

    init_attrs = {
        RATIOS: ratios,
        MIN_SUBSET_SIZE: min_subset_size,
        MAX_SUBSET_SIZE: max_subset_size,
        INDIVIDUAL_INTERACTORS_METADATA: individual_interactors_metadata_path,
        PROP_RESULTS_FILE: prop_results_path,
        SUBSETS_FILE_PATH: path_to_save_subsets_file,
        GENERATED_SUBSETS_DIR: path_to_save_subsets_dir
    }

    suffix = descriptor + "_select_subsets_state"
    gen_subsets_pipeline = NonRepeatingPipeline(suffix, init_attrs=init_attrs, reset_state=reset_state)
    gen_subsets_pipeline.add_steps(
        [
            _gen_subsets_file,
            _gen_ranked_dfs_by_subset
        ]
    )

    return gen_subsets_pipeline


