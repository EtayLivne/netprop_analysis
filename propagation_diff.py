import json
import pandas as pd
from netprop.models import PropagationResultModel
from typing import Callable


def _propagation_diff_df(prop1: PropagationResultModel, prop2: PropagationResultModel,
                        positive_filter: Callable = None) -> pd.Series:
    names = list(prop1.nodes.keys())
    if not positive_filter:
        positive_filter = lambda x: True
    try:
        return pd.Series({n: prop1.nodes[n].liquids["info"] - prop2.nodes.get(n, prop1.nodes[n]).liquids["info"]
                          for n in names if positive_filter(n)})
    except:
        print(prop1.network["suppressed_nodes"][0])
        return pd.Series([])


def propagation_diff_df(reference_propagation: PropagationResultModel, tested_propagation: PropagationResultModel, sort=False):
    knockout_node_id = tested_propagation.network["suppressed_nodes"][0]
    prop_id = f"{knockout_node_id}_knockout"
    prop_series = _propagation_diff_df(tested_propagation, reference_propagation)
    if sort:
        prop_series.sort_values(inplace=True, ascending=True)
    return prop_id, prop_series


def assert_no_positive_diff(prop_file_list: str, reference_file: str):
    ref = PropagationResultModel.parse_file(reference_file)
    bad_ids = []
    for file in prop_file_list:
        res = PropagationResultModel.parse_file(file)
        p_id, series = propagation_diff_df(ref, res)
        if series[series > 0].any():
            bad_ids.append(p_id)

    print(f"{len[bad_ids]} propagations has positive values")


def scores_iter(ref_file: str, tested_files: list[str], sort=False):
    for f in tested_files:
        ref = PropagationResultModel.parse_file(ref_file)
        test = PropagationResultModel.parse_file(f)
        yield propagation_diff_df(ref, test, sort=sort)




if __name__ == "__main__":
    pass