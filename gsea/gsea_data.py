from pathlib import Path
from utils.utils import load_json, dump_json

def get_gsea_data(histogram_path: str, propagations_root_dir: str, diff_exp_genes_path: str, reference_propgation: str):
    data = GSEAData()
    data.reference_propagation = reference_propgation
    data.set_target_pathway_to_diff_exp_genes(diff_exp_genes_path)
    data.propagation_files_from_histogram(histogram_path, propagations_root_dir)
    return data

class GSEAData:
    propagation_files: list[str] = None
    reference_propagation: str = None
    target_pathway: set[str] = None

    @staticmethod
    def _load_histogram(histogram_file_path: str) -> dict[str, str]:
        histogram = load_json(histogram_file_path)
        # Analysis is only performed on nodes that affect at least one differentially expressed gene
        histogram.pop("0")
        return histogram

    @staticmethod
    def _load_diff_exp_genes(diff_exp_genes_file_path) -> set[str]:
        return set(load_json(diff_exp_genes_file_path).keys())

    @staticmethod
    def _propagation_files_from_root_histogram(root_folder, relevent_knockouts: list[str]) -> list[str]:
        return [Path(root_folder) / f"{n}_knockout.json" for n in relevent_knockouts]

    def propagation_files_from_histogram(self, histogram_path: str, propgations_root_dir: str):
        histogram = self._load_histogram(histogram_path)
        self.propagation_files = self._propagation_files_from_root_histogram(propgations_root_dir,
                                                                             sum(list(histogram.values()), []))

    def set_target_pathway_to_diff_exp_genes(self, diff_exp_genes_path: str):
        self.target_pathway = self._load_diff_exp_genes(diff_exp_genes_path)

