from utils.utils import load_json, dump_json
from enum import Enum


class DiffExpGenesSources(Enum):
    STUKALOV = "stukalov"
    BOJKOVA = "bojkova"

def load_bojkova(path: str):
    pass


def load_stukalov(diff_exp_genes_file_path: str):
    return set(load_json(diff_exp_genes_file_path).keys())


_DIFF_EXP_GENES_SOURCE_MAP = {
    "stukalov": load_stukalov,
    "bojkova": load_bojkova
}


def load(filepath: str, source_name: str):
    try:
        return _DIFF_EXP_GENES_SOURCE_MAP[source_name](filepath)
    except KeyError:
        raise ValueError(f"Unsupported source name, choose one of the following: {[e.value for e in DiffExpGenesSources]}")
