from netprop.networks.loaders import NetpropNetwork
from netprop.models import ConfigModel
from utils.utils import load_json, dump_json
from random import choice


def randomly_replace_interactors(prot_to_interactor: str, output_path: str):
    prot_to_inter = load_json(prot_to_interactor)
    keys = list(prot_to_inter.keys()) #transforming to list to preserve order over multiple operations
    vals = [prot_to_inter[k] for k in keys]
    vals = deg_preserving_random_sets(vals)
    dump_json({keys[i]: vals[i] for i in range(len(keys))}, output_path)


def deg_preserving_random_sets(list_of_sets: list[list[str]]) -> list[list[str]]:
    path = r"D:\data\networks\no_covid_network.json"
    nw = NetpropNetwork(path).load()
    deg_dict = dict()
    for node, deg in nw.degree():
        if deg not in deg_dict:
            deg_dict[deg] = []
        deg_dict[deg].append(node)

    randomized_sets = []
    for s in list_of_sets:
        new_set = {None}
        for n in s:
            n_deg = nw.degree(n)
            new_node = None
            while new_node in new_set:
                new_node = choice(deg_dict[n_deg])
            new_set.add(new_node)
        new_set.remove(None)
        randomized_sets.append(list(new_set))

    return randomized_sets


def patch_conf_with_prot_to_interactor(prot_to_interactor: str, original_conf: str, new_conf_output_dir: str, new_conf: str):
    d = load_json(prot_to_interactor)
    prior_nodes = []
    for k,v in d.items():
        for n in v:
            prior_nodes.append({"id": n, "source_of": [k, n]})
    conf = ConfigModel.parse_file(original_conf)
    new_prior_set = [
        {
            "id": "all",
            "confidence": 0.7,
            "nodes": prior_nodes
        }
    ]
    conf.prior_set = new_prior_set
    conf.output_dir_path = new_conf_output_dir
    dump_json(conf.dict(exclude_unset=True), new_conf)


