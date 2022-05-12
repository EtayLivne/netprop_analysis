from metrics.metric import SinglePropMetric, MultiPropMetric, MetricRuntimeError
import pandas as pd

from netprop.models import PropagationResultModel
from score_from_netprop import propagation_to_df
from typing import Union, Callable
import matplotlib.pyplot as plt
import numpy as np

# get_huh_crispr for huh_crispr hit set
class SinglePropROC(SinglePropMetric):
    def __init__(self, prop_file_path=None, by_liquids: list[str]=None, hit_set_load_method: Callable=None):
        super().__init__(prop_file_path=prop_file_path, by_liquids=by_liquids)
        self._hit_set_load_method = hit_set_load_method
        self._hit_set = None

    def load_hit_set(self, data_file_path: str, hit_set_load_method: Callable=None):
        if hit_set_load_method:
            self._hit_set_load_method = hit_set_load_method
        if self._hit_set_load_method is None:
            raise MetricRuntimeError(f"Cannot load hit set from {data_file_path} since no loading method defined")
        self._hit_set = self._hit_set_load_method(data_file_path)

    def calc_metric(self):
        if not self._hit_set:
            raise MetricRuntimeError("cannot measure roc before hit set is defined! (call load_hit_set)")
        elif self._prop_df.empty:
            raise MetricRuntimeError("cannot measure roc before prop results are defined! (call prop_data_to_df")

        self._hit_set = [h for h in self._hit_set if h in list(self._prop_df[self._NODES_COLUMN_NAME])]
        results = dict()
        for c in self._prop_df.columns:
            if c == "nodes":
                continue
            ascending = "p_value" in c
            self._prop_df.sort_values(c, ascending=ascending, inplace=True)
            hit_df = self._prop_df[self._NODES_COLUMN_NAME].isin(self._hit_set)
            positives = hit_df.sum()
            negatives = len(hit_df) - positives
            tpr = hit_df.cumsum() / positives
            fpr = (1-hit_df).cumsum() / negatives
            res = pd.concat([tpr, fpr], axis=1)
            res.columns = ["tpr", "fpr"]
            auc = np.trapz(res["tpr"], res["fpr"])
            results[c] = (res, auc)
        # self._plot_roc(roc_df := pd.concat([tpr, fpr], axis=1))
        return results

    @classmethod
    def show_results(cls, results: dict[str, tuple[pd.DataFrame], float], show: bool=True, save: str=None) -> None:
        if len(results) < 6:
            fig, axes = plt.subplots(ncols=len(results))
        else:
            lines = int(input(f"there are {len(results)} plots to show. how many lines?"))
            cols = int(input(f"there are {lines} lines to show. how many cols?"))
            fig, axes = plt.subplots(nrows=lines, ncols=cols)
        axes_shape = axes.shape
        axes = axes.flatten()
        fig.suptitle("ROC AUC: Krogan PPI, protein specific interactors interactors")
        counter = 0
        for liq, data in results.items():
            df, auc = data
            ax = axes[counter]
            x,y = df["fpr"], df["tpr"]
            ax.plot(x, y)
            ax.set_title(liq)
            ax.set_xlabel("fpr")
            ax.set_ylabel("tpr")
            ax.annotate(f'auc: {auc:.4f}', xy=(0.25, 0.1), xycoords='axes fraction')
            counter += 1
        axes = axes.reshape(axes_shape)
        if show:
            plt.show()
        if save is not None:
            plt.savefig(save)


class AugmentedFeaturesROC(SinglePropROC):
    def _prop_data_to_df(self, res: Union[PropagationResultModel, str], by_liquids: list[str]):
        super()._prop_data_to_df(res, by_liquids)
        self._augment_features(by_liquids)

    # adds columns to self._prop_df according to class intention
    def _augment_features(self, by_liquids: Union[str, list[str]]):
        raise NotImplemented


class LinearClassifierROC(AugmentedFeaturesROC):

    def _augment_features(self, by_liquids: Union[str, list[str]]):
        pass


class NormsROC(AugmentedFeaturesROC):
    def _augment_features(self, by_liquids: Union[str, list[str]]):
        # ignores by_liquid, nothing to add that _prop_data_to_df didn't already do
        score_columns = [c for c in self._prop_df if c != "nodes" and "p_value" not in c]
        for i in range(1, 4):
            self._prop_df[f"L{i}"] = self._prop_df.apply(lambda row: sum(abs(row[score_columns])**i)**(1./i), axis=1)
        self._prop_df["L_inf"] = self._prop_df.apply(lambda row: max(row[score_columns]), axis=1)
        self._prop_df.drop(score_columns, inplace=True, axis=1)

# class MultiPropRoc(MultiPropMetric):
#     def __init__(self, prop_files_path=None):
#         super().__init__(prop_files_path=prop_files_path)
#         self._hit_set_load_method = None
#         self._hit_set = None
#
#     def load_hit_set(self, data_file_path: str):
#         self._hit_set = self._hit_set_load_method(data_file_path)
#
#     def calc_metric(self):
#         if not self._hit_set:
#             raise MetricRuntimeError("cannot measure roc before hit set is defined!")
#         elif self._props_df.empty:
#             raise MetricRuntimeError("cannot measure roc before prop results are defined!")
#
#         hit_df = self._props_df[self._NODES_COLUMN_NAME].isin(self._hit_set)
#         positives = hit_df.sum()
#         negatives = len(hit_df) - positives
#         tpr = hit_df.cumsum() / positives
#         fpr = ~hit_df.cumsum() / negatives
#         res = pd.concat([tpr, fpr], axis=1)
#         res.columns = ["tpr", "fpr"]
#         # self._plot_roc(roc_df := pd.concat([tpr, fpr], axis=1))
#         return res
#
#
# class HuhCrisprPValROC(MultiPropRoc):
#     def __init__(self, prop_files_path=None, data_file_path=None):
#         super().__init__(prop_files_path=prop_files_path)
#         self._hit_set_load_method = get_huh_crispr
#         if data_file_path is not None:
#             self.load_hit_set(data_file_path)
#
#
#
#     def _prop_data_to_df(self, res_list: Union[list[PropagationResultModel], list[str]], by_liquid: str = "info"):
#         self._prop_df = propagation_to_df(res_list, by_liquid=by_liquid)
#         liquid_scores = [c for c in self._prop_df.columns if by_liquid in c]
#         if len(liquid_scores) != len(res_list):
#             raise MetricRuntimeError(f"an object of type MultiPropROC cannot handle pds with more or less than one"
#                                      f" score of matching liquid per prop result.\n"
#                                      f"Following columns all match the liquid {by_liquid}: {liquid_scores}")
#
#         reference_liquid = [c for c in liquid_scores if "original" in c]
#         if len(reference_liquid) != 1:
#             raise MetricRuntimeError(f"there needs to be exactly one prop named original and it needs to have by_liquid"
#                                      f"exactly once")
#         reference_liquid = reference_liquid[0]
#
#         def calc_p_value(row):
#             sorted_values = sorted(row[c] for c in liquid_scores)
#             return 1 + (sorted_values.index(row[reference_liquid]) / len(liquid_scores))
#
#         self._prop_df["p_value"] = self._prop_df.apply(calc_p_value)
#         self._prop_df.sort_values("p_value", ascending=False, inplace=True)