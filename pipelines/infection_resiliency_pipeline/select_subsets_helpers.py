import sys
from os import path
sys.path.append("/code")


from random import sample
from math import floor
from utils.utils import dump_json, load_json
from pathlib import Path
import pandas as pd
import numpy as np
from ray import remote
from ray import init as init_ray
from multiprocessing import cpu_count
from utils.utils import RayQueueManager, QueueWorker, MultiQueueManager
from itertools import chain, combinations
from random import seed as rand_seed
from math import comb


_ACCUMULATOR_KWARGS__CSV_PATH_KEY = "csv_path"

def _powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def _num_subsets(set_size: int, subset_ratio: float,  min_subsets: int, max_subsets: int) -> int:

    num_combs = comb(set_size, round(subset_ratio*set_size))
    num_subsets = min(max_subsets, max(min_subsets, num_combs // 10))
    if num_subsets < 1:
        raise ValueError("HHHOOOOWWWWWwww")
    return min(20, max(5, num_combs // 10))

    # return 20    #TODO actual math based calculation


def _gen_subsets(the_set: list[str], subset_ratio: float, min_subsets: int, max_subsets: int) -> list[str]:
    set_size = len(the_set)
    sample_size = floor(subset_ratio * set_size)
    for _ in range(_num_subsets(set_size, subset_ratio, min_subsets, max_subsets)):
        yield sample(the_set, sample_size)


def _gen_all_subsets(source_set: list[str], shuffle: bool=False) -> list[str]:
    powerset = list(_powerset(source_set))
    if shuffle:
        #shuffle(powerset)
        raise NotImplemented("Note the bug in commented out code!")
    for elem in powerset:
        yield elem


def subsets_file(sets: dict[str, list[str]], subset_ratios: list[float],  min_subsets: int, max_subsets: int,
                 output_path: str,) -> None:
    results = dict()
    for set_name, set_members in sets.items():
        results[set_name] = {ratio: list(_gen_subsets(set_members, ratio, min_subsets, max_subsets)) for ratio in subset_ratios}
    dump_json(results, output_path)


def all_subsets_file(sets: dict[str, list[str]], output_path: str) -> None:
    for set_name, set_elements in sets.items():
        dump_json(list(_gen_all_subsets(set_elements)), output_path + f"\{set_name}.json")
    # dump_json(results, output_path)


@remote(max_restarts=-1)
class SubsetRanker(QueueWorker):
    def __init__(self, task_queue, output_queue, csv_res_path: str):
        super().__init__(task_queue, output_queue)
        self.df = pd.read_csv(csv_res_path, index_col="nodes")

    def _perform_task(self, task_input):
        subset = task_input
        self.df["subset"] = self.df.apply(lambda row: row[subset].sum(), axis=1)
        self.df = self.df.sort_values("subset", ascending=False)
        return pd.Series(np.arange(len(self.df)), index=self.df.index)  # index is protein name,


def _iter_subsets(subsets_file_path: str, output_dir: str) -> tuple:
    all_subsets = load_json(subsets_file_path)
    output_dir_path = Path(output_dir)
    print(f"will do a total of {len(all_subsets)} sets of subsets")
    counter = 0
    stop = False
    for set_name, subset_ratios in all_subsets.items():
        if stop:
            pass
            #break
        set_output_dir = output_dir_path / set_name
        set_output_dir.mkdir(parents=True, exist_ok=True)
        for subset_ratio, subsets in subset_ratios.items():
            if any([len(s) == 0 for s in subsets]):
                print("skipping empty subsets...")
                continue
            if len(subsets) == 0:
                raise ValueError("WHATTT")
            csv_path = str(set_output_dir/str(subset_ratio)) + ".csv"
            yield subsets, {_ACCUMULATOR_KWARGS__CSV_PATH_KEY: csv_path}
            counter += 1
            if counter == 10:
                stop = True
                #break


def _results_to_ratio_df(queue_outputs: list[pd.Series], **kwargs):
    print("Are we even here???????????????? ")
    print(f"got total {len(queue_outputs)} results")
    ratio_df = pd.concat(queue_outputs, axis=1)
    # if len(kwargs) > 0:
    #     raise ValueError(f"kwargs was {kwargs}, but the key {_ACCUMULATOR_KWARGS__CSV_PATH_KEY} of type {type(_ACCUMULATOR_KWARGS__CSV_PATH_KEY)}")
    # else:
    #     raise KeyError
    csv_path = kwargs[_ACCUMULATOR_KWARGS__CSV_PATH_KEY]
    ratio_df.to_csv(csv_path)


def rank_by_subset_df(subsets_file_path: str, prop_res_file: str, output_dir: str, num_workers=min(cpu_count(), 60)) -> None:
    # rand_seed(1042)
    print(f"about to fail with cpu count {cpu_count()} & num workers {num_workers}")
    init_ray(address=None, ignore_reinit_error=True, num_cpus=num_workers)  # configuration for docker to avoid crush on dashboard launch


    mqm = MultiQueueManager(
        6,
        SubsetRanker,
        num_workers // 6,
        _results_to_ratio_df,
        worker_init_args=prop_res_file
    )
    mqm.handle_inputs(_iter_subsets(subsets_file_path, output_dir))



    # all_subsets = load_json(subsets_file_path)
    # output_dir_path = Path(output_dir)

    # qm = RayQueueManager(SubsetRanker, num_workers, prop_res_file)
    # print(f"will do a total of {len(all_subsets)} sets of subsets")
    # counter = 0
    # for set_name, subset_ratios in all_subsets.items():
    #     set_output_dir = output_dir_path / set_name
    #     set_output_dir.mkdir(parents=True, exist_ok=True)
    #     for subset_ratio, subsets in subset_ratios.items():
    #         if any([len(s) == 0 for s in subsets]):
    #             print("skipping empty subsets...")
    #             continue
    #         if len(subsets) == 0:
    #             raise ValueError("WHATTT")
    #         ranked_series = qm.handle_inputs(subsets)
    #         print(f"got total {len(ranked_series)} results")
    #         ratio_df = pd.concat(ranked_series, axis=1)
    #         ratio_df.to_csv(str(set_output_dir/str(subset_ratio)) + ".csv")




