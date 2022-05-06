import networkx as nx
from netprop.models import PropagationResultModel
import matplotlib.pyplot as plt

def draw_prop_network(network: nx.Graph, res_file: str, prior_set_label: str="prior"):
    res = PropagationResultModel.parse_file(res_file)
    prior_set = [n["id"] for n in res.network.prior_set.nodes]
    node_to_label = {n: prior_set_label for n in prior_set}
    node_to_liquid = {
        n: data.liquids["info"] for n, data in res.nodes.items()
    }
    max_liquid = max(node_to_liquid)
    node_to_liquid = {n: v/max_liquid for n, v in node_to_liquid.items()}
    draw_network(network, node_to_label, node_to_liquid)
    # nx.draw(network, node_color=node_to_liquid.values(), cmap=plt.cm.Blues)


def draw_network(network: nx.Graph, label_dct: dict[str, str], node_scores: dict[str, float]):
    nx.draw(network, node_color=list(node_scores.values()), cmap=plt.cm.Blues, labels=label_dct)
    plt.show()