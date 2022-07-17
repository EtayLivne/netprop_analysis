from random import sample
from math import floor
from netprop.models import ConfigModel
from utils.utils import dump_json, load_json
from pathlib import Path
import pandas as pd
import numpy as np
from ray import remote
from ray.util.queue import Queue
from multiprocessing import cpu_count
def _num_subsets(set_size: int, subset_ration: float) -> int:
    return 5    #TODO actual math based calculation


def _gen_subsets(the_set: list[str], subset_ratio: float) -> list[str]:
    set_size = len(the_set)
    sample_size = floor(subset_ratio * set_size)
    for _ in range(_num_subsets(set_size, subset_ratio)):
        yield sample(the_set, sample_size)


#

def subsets_file(sets: dict[str, list[str]], subset_ratios: list[float], output_path: str) -> None:
    results = dict()
    for set_name, set_members in sets.items():
        results[set_name] = {ratio: list(_gen_subsets(set_members, ratio)) for ratio in subset_ratios}
    dump_json(results, output_path)

@remote
class subset_ranker:
    def __init__(self, full_res_path: str):
        self.df = pd.read_csv(full_res_path)

    def specific_subset_ranking(self, subset: list[str]) -> pd.Series:
        self.df["subset"] = self.df.apply(lambda row: row[subset].sum(), axis=1)
        self.df = self.df.sort_values("subset", ascending=False)
        self.df["ranks"] = pd.Series(np.arange(len(self.df)))
        return self.df["ranks"]

def rank_by_subset_df(subsets_file: str, prop_res_file: str, output_dir: str, num_workers=2*(cpu_count() - 2)) -> None:
    subsets = load_json(subsets_file)
    prop_df = pd.read_csv(prop_res_file)
    output_dir_path = Path(output_dir)
    actors = [subset_ranker.remote(prop_res_file) for _ in range(num_workers)]
    q = Queue()

    for set_name, subset_ratios in subsets.items():
        set_output_dir = output_dir_path / set_name
        set_output_dir.mkdir()
        for subset_ratio, subsets in subset_ratios.items():
            set_ratio_df = pd.DataFrame(index=prop_df.index)
            subset_list = list(enumerate(subsets))
            rank_series_list = []
            # nope just do with pool
            while subset_list:
                slice_size = min(len(subset_list), len(rank_series_list))
                rank_series_list.extend([actors[i].specific_subset_ranking.remote(subset_list[i]) for i in range(slice_size)])
                subset_list = subset_list[slice_size:]
            for i, subset in enumerate(subsets):
                prop_df["subset"] = prop_df.apply(lambda row: row[subset].sum(), axis=1)
                prop_df = prop_df.sort_values("subset", ascending=False)
                prop_df["ranks"] = pd.Series(np.arange(len(prop_df)))
                set_ratio_df[i] = prop_df["ranks"]
            set_ratio_df.to_csv(str(set_output_dir/str(subset_ratio)) + ".csv")

