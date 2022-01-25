import json
import numpy as np
import pandas as pd
from netprop.models import PropagationResultModel
from scipy import stats
from pathlib import Path


def make_gct(df: pd.DataFrame, path: str):

    with open(path, 'w') as handler:
        handler.write(f"# 1.3\n"
                      f"{len(df)}\t{len(df.columns)}\n"
                      f"NAME\tDescription\tpropagation_diff\n")
        for i in range(len(df)):
            handler.write(f"{df.iloc[i].names}\tNA\t{df.iloc[i].propagation_diff}\n")


def make_cls(sample_name: str, path:str):
    with open(path, 'w') as handler:
        handler.write(f"1 1 1\n"
                      f"# {sample_name}\n"
                      f"{sample_name}")


def make_gmx(gene_set_name: str, gene_set: list[str], path):
    with open(path, 'w') as handler:
        handler.write(f"{gene_set_name}\n"
                      f"na\n")
        handler.write("\n".join(gene_set))


def make_gmt(gene_set_name: str, gene_set: list[str], path):
    with open(path, 'w') as handler:
        header = f"{gene_set_name}\tna\t"
        handler.write(header + "\t".join(gene_set))