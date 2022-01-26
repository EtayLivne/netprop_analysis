from netprop.models import PropagationResultModel
from multiprocessing import Queue
from utils.utils import load_json, dump_json
from pathlib import Path
from multiprocessing import Pool, Queue, cpu_count
from propagation_diff import propagation_diff_df


def produce_df(reference_file: str, task_q: Queue, output_q: Queue):
    reference = PropagationResultModel.parse_file(reference_file)
    while True:
        task = task_q.get(block=True)
        if task == "HALT":
            break
        tested_propagation_path = task
        tested = PropagationResultModel.parse_file(tested_propagation_path)
        output_q.put(propagation_diff_df(reference, tested))


def process_df(diff_exp_genes: list[str], task_q: Queue, output_q: Queue):
    histogram = [list() for i in range(len(diff_exp_genes) + 1)]
    while True:
        task = task_q.get(block=True)
        if task == "HALT":
            break

        prop_id, prop_df = task
        negatives = prop_df < 0
        histogram[sum([1 for g in diff_exp_genes if negatives.get(g, False)])].append(prop_id)

    output_q.put(histogram)


def produce_histogram(propagation_files_root, reference_file, diff_genes_file, output_path, num_workers):
    print("ya")
    propagation_files = Path(propagation_files_root).glob("*knockout.json")

    diff_exp_genes = list(load_json(diff_genes_file).keys())

    producing_workers = min(1, num_workers // 2 - 1)
    processing_workers = min(1, num_workers // 2 - 1)
    producer_queue = Queue()
    processor_queue = Queue()
    output_queue = Queue()

    for file in propagation_files:
        producer_queue.put(file)

    for i in range(producing_workers):
        producer_queue.put("HALT")

    producer_queue.close()

    producer_pool = Pool(producing_workers, produce_df, [reference_file, producer_queue, processor_queue])
    processor_pool =  Pool(processing_workers, process_df, [diff_exp_genes, processor_queue, output_queue])

    producer_pool.join()
    for i in range(processing_workers):
        processor_pool.put("HALT")
    processor_queue.close()

    for i in range(num_workers):
        processor_queue.put("HALT")
    processor_queue.close()
    processor_pool.join()

    histogram = [list() for i in range(len(diff_exp_genes) + 1)]
    for i in range(processing_workers):
        worker_histogram = output_queue.get()
        for i in range(len(histogram)):
            histogram[i].extend(worker_histogram[i])

    dictogram = {i: histogram[i] for i in range(len(histogram))}
    print("***histogram***")
    for k, v in dictogram.items():
        print(f"{f}: {len(v)}")
    dump_json(dictogram, output_path)

    return histogram
    # plt.figure()
    # plt.bar(list(range(len(histogram))), [len(i) for i in histogram])
    # plt.show()