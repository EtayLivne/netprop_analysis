
from multiprocessing import Pool, cpu_count
import pandas as pd
from typing import Union
from abc import abstractmethod, ABC
from netprop.models import PropagationResultModel


class MetricRuntimeError(Exception):
    pass


class Metric(ABC):
    @abstractmethod
    def calc_metric(self):
        raise NotImplemented


class SinglePropMetric(Metric):
    _NODES_COLUMN_NAME = "nodes"
    _SCORE_COLUMN_NAME = "scores"

    def __init__(self, prop_file_path=None):
        if prop_file_path is None:
            self._prop_df = pd.DataFrame(columns=[self._NODES_COLUMN_NAME, self._SCORE_COLUMN_NAME])
        else:
            self.prop_data_to_df(prop_file_path)

    def prop_data_to_df(self, prop_file_path: str):
        res = PropagationResultModel.parse_file(prop_file_path)
        self._prop_data_to_df(res)

    @abstractmethod
    def _prop_data_to_df(self, res: Union[PropagationResultModel, str], by_liquid: str="info"):
        self._prop_df = None


class MultiPropMetric(Metric):
    def __init__(self):
        self._props_df = None

    def props_data_to_df(self, prop_file_paths: list[str], max_proc=min(cpu_count(), 40)):
        with Pool(min(len(prop_file_paths), max_proc)) as p:
            res_list = p.map(PropagationResultModel.parse_file, prop_file_paths)
        self._prop_data_to_df(res_list)

    @abstractmethod
    def _prop_data_to_df(self, res_list: Union[list[PropagationResultModel], list[str]], by_liquid: str="info"):
        self._props_df = None
