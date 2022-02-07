from abc import ABCMeta, abstractmethod
from gsea.gsea_data import GSEAData
from utils.utils import dump_json


class AbstractGSEAAnalysis(metaclass=ABCMeta):
    @abstractmethod
    def _analyze(self, data: GSEAData) -> dict:
        raise NotImplemented

    def _validate_data(self, data: GSEAData):
        return data.propagation_files and data.reference_propagation and data.target_pathway

    def analyze(self, data: GSEAData, output_path: str):
        if not self._validate_data(data):
            print(f"prop files: {data.propagation_files}\n "
                  f"ref prop: {data.reference_propagation}\n"
                  f"pathway: {data.target_pathway}")
            raise ValueError("YOUR DATA IS BAD ETAY")
        try:
            output = self._analyze(data)
        except Exception as e:
            print("error while attempting analysis")
            raise e

        dump_json(output, output_path)