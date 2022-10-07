
from pipeline import NonRepeatingPipeline
from pipelines.infection_resiliency_pipeline.attrs_keys import *
from pipelines.infection_resiliency_pipeline.resiliency_scores_helpers import test_resiliency_scores


def gen_resiliency_scores(attrs: dict) -> None:
    prop_res_file = attrs[PROP_RESULTS_FILE]
    generated_subsets_dir = attrs[GENERATED_SUBSETS_DIR]
    output_dir = attrs[GENERATED_RESILIENCY_SCORES_DIR]
    top_propagated_threshold = attrs[TOP_PROPAGATED_RANK_THRESHOLD]

    test_resiliency_scores(top_propagated_threshold, prop_res_file, generated_subsets_dir, output_dir)


def get_pipeline(top_propagated_threshold,
                 path_to_generate_resiliency_scores: str,
                 descriptor: str,
                 reset_state: bool=False) -> NonRepeatingPipeline:
    init_attrs = {
        TOP_PROPAGATED_RANK_THRESHOLD: top_propagated_threshold,
        GENERATED_RESILIENCY_SCORES_DIR: path_to_generate_resiliency_scores
    }
    suffix = descriptor + "_resiliency_scores_pipeline"
    resiliency_scores_pipeline = NonRepeatingPipeline(suffix, init_attrs=init_attrs, reset_state=reset_state)
    resiliency_scores_pipeline.add_steps(
        [
            test_resiliency_scores
        ]
    )

    return resiliency_scores_pipeline
