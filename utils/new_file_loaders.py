from netprop.networks.loaders.single import HSapiensNetworkLoader
from netprop.generic_utils.constants import SpeciesIDs
from utils.utils import load_json
import pandas as pd

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


class CoVToHumanMeta(HSapiensNetworkLoader):
    _ALL_CELL_LINES = ["Huh7", "Huh7.5", "A549 ACE2", "Calu-3", "Caco-2", "Vero", "Vero-E6",
                       "293T-ACE2", "HEKT293T", "HEKT293T-ACE2","A549"]
    _APPROVED_CELL_LINES = ["HEK293T", "HEK293T-ACE2", "293T-ACE2"]

    def __init__(self, human_network_path: str,
                 protein_interactions_path: str=None, rna_interactions_path: str=None,
                 merged_covid=False):

        super().__init__(human_network_path)

        if not (protein_interactions_path or rna_interactions_path):
            raise ValueError("Cov to human meta loader may load proteins, rna or both, but neither were provided")

        self.protein_interactions_path = protein_interactions_path
        self.rna_interactions_path = rna_interactions_path
        self.merged_covid = merged_covid

    def load(self, *args, **kwargs):
        network = super().load()
        interactions = sum([self._load_rna_interactions(), self._load_protein_interactions()], [])

        for cov_source, human_target in interactions:
            human_target = str(human_target)
            if cov_source not in network:
                network.add_node(cov_source, **self._new_node_attrs(SpeciesIDs.CORONAVIRUS.value))
            if human_target not in network:
                print(f"cannot add covid interaction with {human_target} - it isn't in the network")
                continue
            confidence = 0.9
            network.add_edge(cov_source, human_target, weight=confidence)
        return network

    def _load_protein_interactions(self):
        if not self.protein_interactions_path:
            return []


        id_col = "entrez"
        df = pd.read_csv(self.protein_interactions_path,
                         usecols=["Assay cell line", "Bait", "Reference", id_col], index_col="Assay cell line")
        df.dropna(inplace=True)
        df[id_col] = df[id_col].astype("int32")
        interactions = []
        # if self.merged_covid:
        #     df["Bait"] = df["Bait"].apply(lambda name: "covid_proteins")

        validation_counter = dict()
        for cell_line in set(df.index):
            if cell_line not in self._APPROVED_CELL_LINES:
                continue
            for _, row in df.loc[cell_line].iterrows():
                key = (row["Bait"], row[id_col])
                validation_counter[key] = validation_counter.get(key, 0) + 1

            interactions.extend([k for k in validation_counter if validation_counter[k] > 1])

        if self.merged_covid:
            interactions = [("covid_proteins", x[1]) for x in interactions]

        return interactions

    def _load_rna_interactions(self):
        if not self.rna_interactions_path:
            return []
        id_col = "entrez"
        df = pd.read_csv(self.rna_interactions_path,
                         usecols=["Assay cell line", id_col])
        df[id_col].astype("int32")
        df.set_index(id_col)
        interactions = []
        for cell_line in set(df.index):
            if cell_line not in self._APPROVED_CELL_LINES:
                continue
            interactions.extend([("covid_rna", row[id_col]) for _, row in df.loc[cell_line].iterrows()])

        return interactions
