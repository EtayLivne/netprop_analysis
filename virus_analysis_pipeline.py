from utils.utils import load_json, dump_json
from netprop.models import ConfigModel, PropagationResultModel
from netprop.propagation.api import propagate_from_config
import pandas as pd
from pathlib import Path
from random import shuffle
from math import log
from typing import Union
from functools import partial


# calculated so that there's more than a 95% chance for every pair of proteins to be in opposiong halves at least once.
def _calc_num_splits(set_size: int) -> int:
    log_base = (set_size - 2)/(set_size - 1)
    try:
        return int(log(0.1, log_base)) + 1
    except:
        raise ValueError(f"log_base: {log_base}")


def _split_a_list(to_split: list) -> tuple[list, list]:
    shuffled = list(to_split)
    shuffle(shuffled)
    return shuffled[:len(shuffled)//2], shuffled[len(shuffled)//2:]


def _randomly_split(prot_to_interactors: dict, proteins: Union[str, list[str]]=None) -> dict[list[list[str]]]:
    # prot_map = load_json(prot_to_interactors)

    if proteins is None:
        proteins =list(prot_to_interactors.keys())
    elif not isinstance(proteins, list):
        proteins = [proteins]

    prot_to_interactors["all"] = list(set().union(*[v for v in prot_to_interactors.values()]))
    return {
        protein: [_split_a_list(prot_to_interactors[protein]) for
                  i in range(_calc_num_splits(len(prot_to_interactors[protein])))]
        for protein in proteins if len(prot_to_interactors[protein]) >=10
    }

def create_splits_files(viral_prot_to_interactor: str, splits_file_path: str):
    prot_to_interactor = load_json(viral_prot_to_interactor)
    splits_dict = _randomly_split(prot_to_interactor)
    dump_json(splits_dict, splits_file_path)


