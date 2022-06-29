from netprop.models import ConfigModel, PropagationResultModel
from utils.utils import load_json, dump_json
from pathlib import Path
from netprop.propagation.api import propagate_from_config
import pandas as pd

def create_metadata_file(prior_set_source: str, metadata_file_path: str) -> None:
    interactors_series = pd.read_csv(prior_set_source, index_col="Bait")["entrez"]
    dict_of_lists = dict()
    for item in interactors_series.iteritems():
        viral_protein, interactor = item[0], item[1]
        if viral_protein not in dict_of_lists:
            dict_of_lists[viral_protein] = []
        dict_of_lists[viral_protein].append(str(interactor))

    dump_json(dict_of_lists, metadata_file_path)


def create_virus_base_conf(base_conf: str, viral_prot_to_interactor: str, new_output_dir_path: str, new_id: str, new_conf_path: str):
    conf = ConfigModel.parse_file(base_conf)
    interactors_dict = load_json(viral_prot_to_interactor)
    nodes = []
    for viral_protein, interactors in interactors_dict.items():
        for interactor in interactors:
            nodes.append({"id": interactor, "source_of": [viral_protein, interactor]})
    conf.prior_set = [
        {
            "id": "all",
            "confidence": 0.7,
            "nodes": nodes
        }
    ]
    conf.id = new_id
    conf.output_dir_path = new_output_dir_path
    dump_json(conf.dict(exclude_unset=True), new_conf_path)

def create_randomized_conf(base_conf: str, randomized_network_conf: str, new_conf_path: str):
    conf = ConfigModel.parse_file(base_conf)
    conf.networks = {"path": str(randomized_network_conf)}
    conf.output_dir_path = str(conf.output_dir_path + "\\randomized")
    dump_json(conf.dict(exclude_unset=True), new_conf_path)


def prop_from_conf_and_save_new_format(conf_path: str):
    conf = ConfigModel.parse_file(conf_path)
    output_dir_path = conf.output_dir_path
    res_json_path = Path(output_dir_path) / "all.json"
    res_csv_path = Path(output_dir_path) / "all.csv"
    propagate_from_config(conf_path, ordering={"prior_set":100})
    _prop_result_to_big_ass_df(res_json_path, res_csv_path)
    return res_csv_path


def _prop_result_to_big_ass_df(res_path: str, output_path: str):
    res = PropagationResultModel.parse_file(res_path).prop_scores_as_series(by_liquids=None)
    res = res.set_index("nodes")
    res = res.rename(columns={c: c.split(".")[-1] for c in list(res.columns)})
    res.to_csv(output_path)