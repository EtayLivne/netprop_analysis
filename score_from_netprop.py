import json
import pandas as pd
from netprop.models import PropagationResultModel
from typing import Callable
from scripts.merge_fluids import merge_all_liquids

def _propagation_diff_to_df(prop1: PropagationResultModel, prop2: PropagationResultModel,
                            positive_filter: Callable = None, liquid_name: str="info") -> pd.Series:
    names = list(prop1.nodes.keys())
    if not positive_filter:
        positive_filter = lambda x: True
    try:
        return pd.Series({n: prop1.nodes[n].liquids[liquid_name] - prop2.nodes.get(n, prop1.nodes[n]).liquids[liquid_name]
                          for n in names if positive_filter(n)})
    except:
        print(prop1.network["suppressed_nodes"][0])
        return pd.Series([])


# When propagation diff is sorted, knockouts that cause greater negative diff are first, so order is ascending by diff.
def propagation_diff_to_df(reference_propagation: PropagationResultModel, tested_propagation: PropagationResultModel,
                           liquid_name: str="liquid", sort=False):
    knockout_node_id = tested_propagation.network["suppressed_nodes"][0]
    prop_id = f"{knockout_node_id}_knockout"
    prop_series = _propagation_diff_to_df(tested_propagation, reference_propagation, liquid_name=liquid_name)
    if sort:
        prop_series.sort_values(inplace=True, ascending=True)
    return prop_id, prop_series


def assert_no_positive_diff(prop_file_list: str, reference_file: str):
    ref = PropagationResultModel.parse_file(reference_file)
    bad_ids = []
    for file in prop_file_list:
        res = PropagationResultModel.parse_file(file)
        p_id, series = propagation_diff_to_df(ref, res)
        if series[series > 0].any():
            bad_ids.append(p_id)

    print(f"{len[bad_ids]} propagations has positive values")


def scores_iter(ref_file: str, tested_files: list[str],
                liquid_name: str="info", merge_liquids: bool=False, sort=False):
    ref = PropagationResultModel.parse_file(ref_file)
    if merge_liquids:
        merge_all_liquids(ref)
    for f in tested_files:
        test = PropagationResultModel.parse_file(f)
        if merge_liquids:
            merge_all_liquids(test)
        yield propagation_diff_to_df(ref, test, liquid_name=liquid_name, sort=sort)


# When propagation is sorted, nodes with higher liquids will be first (order is descending)
def propagation_to_df(propagation: PropagationResultModel, by_liquid: str="info", sort: bool=False):
    prop_id = propagation.id
    prop_series = pd.Series({n: propagation.nodes[n].liquids[by_liquid] for n in propagation.nodes})
    if sort:
        prop_series.sort_values(ascending=False, inplace=True)
    return prop_id, prop_series


if __name__ == "__main__":
    pass