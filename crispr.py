import json
import pandas as pd


def get_crispr_rankings(rankning_path, translator_path):
    with open(translator_path, 'r') as handler:
        trans_dict = json.load(handler)
    series = pd.read_csv(rankning_path, index_col="id", usecols=["id", "pos|score"])["pos|score"]
    series.rename("crispr_score", inplace=True)
    to_drop = [symbol for symbol in series.index if symbol.upper() not in trans_dict]
    series.drop(labels=to_drop, inplace=True)
    series.rename(index={symbol: trans_dict.get(symbol.upper()) for symbol in series.index}, inplace=True)
    for dup in series[series.index.duplicated()].index.to_list():
        dup_mean = series[dup].mean()
        series.drop(dup, inplace=True)
        series[dup] = dup_mean

    return series