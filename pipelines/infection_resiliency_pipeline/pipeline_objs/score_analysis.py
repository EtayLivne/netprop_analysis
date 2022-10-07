from pipeline import NonRepeatingPipeline
from pipelines.infection_resiliency_pipeline.attrs_keys import *
from pipelines.infection_resiliency_pipeline.score_analysis_helpers import analyze_by_mean_offset


# top_ranked_threshold: int, result_dir_path: str, output_file_path: str
def analyze(attrs: dict) -> None:
    threshold = attrs[TOP_RESILIENCE_RANK_THRESHOLD]
    generated_resiliency_scores_dir = attrs[GENERATED_RESILIENCY_SCORES_DIR]
    generated_resilient_genes_path = attrs[GENERATED_TOP_RESILIENT_GENES_PATH]

    analyze_by_mean_offset(threshold, generated_resiliency_scores_dir, generated_resilient_genes_path)


def get_pipeline(top_reseilience_rank_threshold: int,
                 path_to_save_gene_list: str,
                 descriptor: str, reset_state: bool=False) -> NonRepeatingPipeline:
    init_attrs = {
        TOP_RESILIENCE_RANK_THRESHOLD: top_reseilience_rank_threshold,
        GENERATED_TOP_RESILIENT_GENES_PATH: path_to_save_gene_list
    }

    suffix = descriptor + "score_analysis_pipeline"
    score_analysis_pipeline = NonRepeatingPipeline(suffix, init_attrs=init_attrs, reset_state=reset_state)
    score_analysis_pipeline.add_steps(
        [
            analyze
        ]
    )

    return score_analysis_pipeline
