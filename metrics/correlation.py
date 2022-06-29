from .metric import SinglePropMetric
import pandas as pd
from scipy.stats import spearmanr

# Only works when handed a csv for now
class CorrelationMetric(SinglePropMetric):
    def _single_pair_spearman(self, s1: pd.Series, s2: pd.Series):
        return

    def _calc_pairwise_spearman(self):
        liq_columns = [c for c in self._prop_df.columns if c != "nodes"]
        correlation_df = pd.DataFrame(columns=liq_columns, index=liq_columns)
        for i in range(len(liq_columns) - 1):
            liq1 = liq_columns[i]

            for j in range(i + 1, len(liq_columns)):
                liq2 = liq_columns[j]
                print(f"calculating for {liq1}, {liq2}")
                correlation_df.at[liq1, liq2] = spearmanr(self._prop_df[liq1], self._prop_df[liq2]).correlation
                correlation_df.at[liq2, liq1] = correlation_df.at[liq1, liq2]

        return correlation_df

    def calc_metric(self):
        return self._calc_pairwise_spearman()