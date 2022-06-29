
from multiprocessing import Pool, cpu_count
import pandas as pd
from typing import Union, Callable
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

    def __init__(self, prop_file_path=None, by_liquids: list[str]=None):
        if prop_file_path is None:
            self._prop_df = pd.DataFrame(columns=[self._NODES_COLUMN_NAME, self._SCORE_COLUMN_NAME])
        else:
            self.prop_data_to_df(prop_file_path, by_liquids)

    def prop_data_to_df(self, prop_file_path: str, by_liquids: list[str]):
        # res = PropagationResultModelationResultModel.parse_file(prop_file_path)
        self._prop_data_to_df(prop_file_path, by_liquids=by_liquids)

    def _prop_data_to_df(self, res: Union[PropagationResultModel, str], by_liquids: list[str]):
        if by_liquids is not None and not isinstance(by_liquids, list):
            by_liquids = [by_liquids]

        if type(res) is str:
            extension = res.split(".")[-1]
            if extension == "json":
                res = PropagationResultModel.parse_file(res)
                self._prop_df = res.prop_scores_as_series(by_liquids=by_liquids)
            elif extension == "csv":
                self._prop_df = pd.read_csv(res, usecols=by_liquids)
            else:
                raise ValueError(f"cannot load prop dataframe from unknown file type {extension}")

        elif isinstance(res, PropagationResultModel):
            self._prop_df = res.prop_scores_as_series(by_liquids=by_liquids)

        elif isinstance(res, pd.DataFrame):
            self._prop_df = res

        else:
            raise ValueError("can only load dataframe from a file or from a PropagationResultModel.")

    @abstractmethod
    def calc_metric(self):
        raise NotImplemented


class MultiPropMetric(Metric):
    _NODES_COLUMN_NAME = "nodes"

    def __init__(self, prop_files_path=None):
        if prop_files_path is None:
            self._props_df = pd.DataFrame(columns=[self._NODES_COLUMN_NAME])
        else:
            self.props_data_to_df(prop_files_path)

    def props_data_to_df(self, prop_file_paths: list[str], max_proc=min(cpu_count(), 40)):
        with Pool(min(len(prop_file_paths), max_proc)) as p:
            res_list = p.map(PropagationResultModel.parse_file, prop_file_paths)
        self._prop_data_to_df(res_list)

    @abstractmethod
    def _prop_data_to_df(self, res_list: Union[list[PropagationResultModel], list[str]], by_liquid: str="info"):
        self._props_df = None
