from utils.new_file_loaders import CovToHumanMeta
from utils.queue_managers import dump_json, load_json
from netprop.models import ConfigModel, PropagationSourceModel, PriorSetModel
import networkx as nx

def patch_prior_set_by_network_tag(network: nx.Graph, prior_set_tag: dict[str,str], existing_conf_file: str, output_path: str,
                                   source_of: list[str] = None, confidence: float = 0.7, prior_set_id: str=None):

    conf = ConfigModel.parse_file(existing_conf_file)
    prior_set = conf.prior_set if type(conf.prior_set) is list else [conf.prior_set]
    tag_key = list(prior_set_tag.keys())[0]
    tag_val = prior_set_tag[tag_key]
    not_tag = not bool(tag_val)  # not_tag is promised to have a different value from  tag
    if not prior_set_id:
        prior_set_id = tag_val
    if not source_of:
        source_of = ["info"]
    prior_nodes_by_tag = [PropagationSourceModel(id=n, source_of=source_of)
                          for n, data in network.nodes(data=True)
                          if data.get(tag_key, not_tag) == tag_val]

    prior_set.append(PriorSetModel(confidence=confidence, nodes=prior_nodes_by_tag, id=prior_set_id))
    conf.prior_set = prior_set
    dump_json(conf.dict(exclude_none=True), output_path)

def transfer_suppressed_sets(conf_to, conf_from, output_file: str=None):
    conf1 = ConfigModel.parse_file(conf_to)
    conf2 = load_json(conf_from)
    conf1.suppressed_set = conf2["suppressed_set"]
    output_file = output_file or conf_to
    dump_json(conf1.dict(exclude_none=True), output_file)
