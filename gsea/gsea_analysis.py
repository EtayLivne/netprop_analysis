from gsea.abstract_gsea import AbstractGSEAAnalysis
from gsea.gsea_data import GSEAData
from score_from_netprop import scores_iter
import pandas as pd
import numpy as np
from multiprocessing import Pool, Queue, Process, cpu_count
from time import sleep


def gene_set_enrichment(sorted_genes: pd.Series, pathway: set[str], ro=1):
    len_genes, len_pathway = len(sorted_genes), len(pathway)

    ro_vec = ro ** np.arange(len_genes)
    score_values = {True: np.sqrt((len_genes - len_pathway)/len_pathway),
                     False: -1 * np.sqrt(len_pathway/(len_genes - len_pathway))}
    hits_to_scores = lambda gene: score_values[gene in pathway]
    positional_scores = np.vectorize(hits_to_scores)(sorted_genes.index.to_numpy()) * ro_vec

    running_score, max_score = 0, 0
    leading_edge_index = 0
    for index, score in enumerate(positional_scores):
        running_score += score
        if running_score > max_score:
            max_score = running_score
            leading_edge_index = index

    leading_edge_hits = positional_scores[:leading_edge_index + 1] > 0
    return max_score, sorted_genes[:leading_edge_index + 1][leading_edge_hits].index.to_list()


class GSEAAnalysis(AbstractGSEAAnalysis):
    def _analyze(self, data: GSEAData, liquid: str = "info") -> dict:
        output = dict()
        for prop_id, prop_series in scores_iter(data.reference_propagation, data.propagation_files,
                                                liquid_name=liquid, sort=True):
            output[prop_id] = gene_set_enrichment(prop_series, data.target_pathway)
        return output


class MultiprocessGSEAAnalysis(AbstractGSEAAnalysis):
    def __init__(self, num_processes: int = cpu_count()):
        super().__init__()
        self.num_processes = num_processes

    @staticmethod
    def _score_producing_worker(reference_file, tested_files: list[str], liquid: str, task_q: Queue):
        for prop_id, prop_series in scores_iter(reference_file, tested_files, liquid_name=liquid,  sort=True):
            task_q.put((prop_id, prop_series))
        sleep(60)
        task_q.put("HALT")

    @staticmethod
    def _analysis_worker(diff_ex_genes: set[str], task_q: Queue, out_q: Queue):
        counter = 0
        while True:
            task = task_q.get(block=True)
            if task == "HALT":
                break
            prop_id, prop_series = task
            out_q.put({prop_id: gene_set_enrichment(prop_series, diff_ex_genes)})
            counter += 1
            if counter % 10 == 0:
                print(f"look at me! printed {counter} in total!")
        out_q.put("DONE")

    def _analyze(self, data: GSEAData, liquid: str = "info") -> dict:
        task_q, out_q = Queue(), Queue()
        num_workers = (self.num_processes // 2) - 1
        num_score_producing_workers = num_workers
        chunk_size = len(data.propagation_files) // num_score_producing_workers
        score_producing_workers = [Process(target=self._score_producing_worker,
                                           args=(data.reference_propagation,
                                                 data.propagation_files[i * chunk_size: (i + 1) * chunk_size],
                                                 liquid,
                                                 task_q))
                                   for i in range(num_score_producing_workers)]
        for worker in score_producing_workers:
            worker.start()

        analysis_workers = num_workers
        pool = Pool(self.num_processes, self._analysis_worker, (data.target_pathway, task_q, out_q))
        output = dict()
        done_count = analysis_workers
        while done_count > 0 and len(output) < len(data.propagation_files):
            message = out_q.get()
            if message == "DONE":
                done_count -= 1
                continue
            else:
                output.update(message)

        return output
