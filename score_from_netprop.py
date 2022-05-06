import json
import pandas as pd
from netprop.models import PropagationResultModel
from typing import Callable, Union
from scripts.merge_fluids import merge_all_liquids
from multiprocessing import Pool, cpu_count
from functools import reduce


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


def propagation_to_df(propagation: Union[PropagationResultModel, str], by_liquid: str="info", drop_na=True):
    return propagations_to_df([propagation], by_liquid=by_liquid, drop_na=drop_na)


# When propagation is sorted, nodes with higher liquids will be first (order is descending)
def propagations_to_df(propagations: Union[list[PropagationResultModel], list[str]], by_liquid: str="info", drop_na=True):

    # if given list of file paths, convert to list of PropagationResultModels
    if isinstance(propagations[0], str):
        with Pool(max(cpu_count() - 2, 1)) as pool:
            props = pool.map(PropagationResultModel.parse_file, propagations)
    else:
        props = propagations

    how = "inner" if drop_na else "outer"
    dfs = [p.prop_scroes_as_series(by_liquids=by_liquid) for p in props]
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on='nodes', how=how), dfs)
    return df


def infer_p_value(prop_df: pd.DataFrame, ref_prop: str, by_liquid: str="info"):
    #columns = filter(lambda c: c.split(".")[0] != ref_prop and c.split(".")[1] == by_liquid, prop_df.columns)
    columns = filter(lambda c: c.split(".")[0] != ref_prop and "info" in c, prop_df.columns)
    columns = list(columns)

    p_values = prop_df.apply(lambda row: (1 + (row[columns] > row[ref_prop]).sum()) / len(row), axis=1)
    p_values.rename("p_values", inplace=True)
    new_df = pd.DataFrame({"nodes": prop_df["nodes"], "p_values": p_values})
    return new_df


# calcs subtract_from - subtract
def infer_prop_diff(prop_df: pd.DataFrame, subtract_from: str, subtract: str, by_liquid="info"):
    col1 = f"{subtract_from}.{by_liquid}"
    col2 = f"{subtract}.{by_liquid}"
    return col1 - col2


if __name__ == "__main__":
    pass