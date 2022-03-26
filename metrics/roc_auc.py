from metrics.metric import SinglePropMetric, MetricRuntimeError
import pandas as pd
from crispr import get_huh_crispr
from netprop.models import PropagationResultModel
from score_from_netprop import propagation_to_df
from typing import Union


class SinglePropROC(SinglePropMetric):

    def __init__(self, prop_file_path=None):
        super().__init__(prop_file_path=prop_file_path)
        self._hit_set_load_method = None
        self._hit_set = None

    def load_hit_set(self, data_file_path: str):
        self._hit_set = self._hit_set_load_method(data_file_path)

    def calc_metric(self):
        if not self._hit_set:
            raise MetricRuntimeError("cannot measure roc before hit set is defined!")
        elif self._prop_df.empty:
            raise MetricRuntimeError("cannot measure roc before prop results are defined!")

        hit_df = self._prop_df[self._NODES_COLUMN_NAME].isin(self._hit_set)
        positives = hit_df.sum()
        negatives = len(hit_df) - positives
        tpr = hit_df.cumsum() / positives
        fpr = ~hit_df.cumsum() / negatives

        # self._plot_roc(roc_df := pd.concat([tpr, fpr], axis=1))
        return pd.concat([tpr, fpr], axis=1)

    #TODOD
    # def _plot_roc(self, roc_df):
    #     plt.plot()


class HuhCrisprROC(SinglePropROC):
    def __init__(self, prop_file_path=None, data_file_path=None):
        super().__init__(prop_file_path=prop_file_path)
        self._hit_set_load_method = get_huh_crispr
        if data_file_path is not None:
            self.load_hit_set(data_file_path)

    def _prop_data_to_df(self, res: Union[PropagationResultModel, str], by_liquid: str="info"):
        self._prop_df = propagation_to_df(res, by_liquid=by_liquid)
        liquid_scores = [c for c in self._prop_df.columns if by_liquid in c]
        if len(liquid_scores) != 1:
            raise MetricRuntimeError(f"an object of type SinglePropROC cannot handle pds with more or less than one"
                                     f" score of matching liquid.\n"
                                     f" Follwing columns all match the liquid {by_liquid}: {liquid_scores}")
        score_column = liquid_scores[0]
        self._prop_df.sort_values(score_column, ascending=False, inplace=True)