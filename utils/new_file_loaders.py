from netprop.networks.loaders.single import HSapiensNetworkLoader
from netprop.generic_utils.constants import SpeciesIDs

class CovToHumanStukalov(HSapiensNetworkLoader):
    MERGED_COVID_NODE_NAME = "covid"
    CONFIDENCE = 0.9

    def __init__(self, human_network_path: str, covid_to_human_path: str, merged_covid=True):
        super().__init__(human_network_path)
        self.covid_to_human_path = covid_to_human_path
        self.merged_covid = merged_covid

    def load(self, *args, **kwargs):
        network = super().load()
        with open(self.covid_to_human_path, 'r') as handler:
            cov_to_human_edges = [edge.strip().split('\t') for edge in handler.readlines()]

        for edge in cov_to_human_edges:
            cov_protein, interacting_human_proteins = edge[0], edge[1]
            source = self.MERGED_COVID_NODE_NAME if self.merged_covid else cov_protein
            if source not in network:
                network.add_node(source, **self._new_node_attrs(SpeciesIDs.CORONAVIRUS.value))
            target = interacting_human_proteins
            if target not in network:
                print(f"cannot add covid interaction with {target} - it isn't in the network")
                continue
            confidence = self.CONFIDENCE
            network.add_edge(source, target, weight=confidence)
        return network


