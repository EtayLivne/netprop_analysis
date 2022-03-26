from metrics.metric import SinglePropMetric, MetricRuntimeError
import pandas as pd
import matplotlib.pyplot as plt
from crispr import get_huh_crispr

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
        elif not self._prop_df:
            raise MetricRuntimeError("cannot measure roc before prop results are defined!")

        hit_df = self._prop_df[self._GENE_COLUMN_NAME].isin(self._hit_set)
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