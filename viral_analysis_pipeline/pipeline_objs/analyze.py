from visulizations import visualize_intersection_statistics
from viral_analysis_pipeline.pipeline_objs.pipeline_objs_consts import *
from pathlib import Path
from pipelines import NonRepeatingPipeline


def _visualize_intersection_statistics(attrs: dict):
    virus_name = attrs[VIRUS_NAME_ATTR]
    intersetction_res_root = Path(attrs[VIRUS_RES_ROOT_ATTR]) / r"intersections_data/intersection_results"

    intersection_res_root_inter = intersetction_res_root / "inter"
    intersection_res_root_cross = intersetction_res_root / "cross"
    intersection_res_root_by_size = intersetction_res_root / "by_size"
    suptitle_str = f"{virus_name}: frequency of node appearing in multiple_intersection"
    visualize_intersection_statistics(intersection_res_root_inter, [20, 100], suptitle_str)
    visualize_intersection_statistics(intersection_res_root_cross, [20, 100], suptitle_str)
    visualize_intersection_statistics(intersection_res_root_by_size, [20, 100], suptitle_str)


def get_pipeline(virus_name: str, reset_state: bool=False) -> NonRepeatingPipeline:
    pipe = NonRepeatingPipeline(virus_name + "_analysis_state", reset_state=reset_state)
    pipe.add_step(_visualize_intersection_statistics)
    return pipe
