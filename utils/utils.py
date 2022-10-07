import json
from ray.util.queue import Queue as RayQueue
from multiprocessing import Queue as MPQueue
from ray import remote
import pandas as pd
import numpy as np
from typing import Union
from abc import ABC, abstractmethod
from ray.util.queue import Queue, Empty

def load_json(file_path):
    with open(file_path, 'r') as handler:
        d = json.load(handler)
    return d


def dump_json(obj: dict, file_path: str):
    with open(file_path, 'w') as handler:
        json.dump(obj, handler, indent=4)


class QueueWorker:
    _STOP_TASK_QUEUE_MESSAGE = "HALT"

    def __init__(self, task_queue: Union[RayQueue, MPQueue], output_queue: Union[RayQueue, MPQueue]):
        self.task_queue = task_queue
        self.output_queue = output_queue
        self._stop_task = False

    def _get_queue_message(self):
        m = self.task_queue.get()
        if isinstance(m, str) and m == self._STOP_TASK_QUEUE_MESSAGE:
            self._stop_task = True
        return m

    def _perform_task(self, task_input):
        raise NotImplemented

    def start(self):
        while True:
            task_input = self._get_queue_message()
            if self._stop_task:
                self._stop_task = False
                break
            else:
                self.output_queue.put(self._perform_task(task_input))


class QueueManager(ABC):
    def __init__(self, worker_class, num_workers: int, worker_init_args: Union[list, any]=None, use_ray: bool=True):
        if worker_init_args and not isinstance(worker_init_args, list):
            worker_init_args = [[worker_init_args]] * num_workers
        elif worker_init_args is None:
            worker_init_args = [[]] * num_workers

        self.task_queue = self._init_queue()
        self.output_queue = self._init_queue()
        print(f"worker init args: {len(worker_init_args)}")
        self.workers = self._init_workers(worker_class, worker_init_args)

        self._stop_task_queue_message = worker_class._STOP_TASK_QUEUE_MESSAGE

    @abstractmethod
    def _init_queue(self):
        raise NotImplemented

    @abstractmethod
    def _init_workers(self, worker_class, worker_init_args):
        raise NotImplemented

    @abstractmethod
    def _start_workers(self):
        raise NotImplemented

    @abstractmethod
    def _get_single_output(self):
        raise NotImplemented

    def fill_task_queue(self, task_inputs):
        for task_input in task_inputs:
            self.task_queue.put(task_input)
        for _ in range(len(self.workers)):
            self.task_queue.put(self._stop_task_queue_message)

    def handle_inputs(self, task_inputs: list):
        self.fill_task_queue(task_inputs)
        self._start_workers()
        task_outputs = []
        while len(task_outputs) < len(task_inputs):
            task_outputs.append(self._get_single_output())
        return task_outputs


class RayQueueManager(QueueManager):
    def _init_queue(self):
        return RayQueue()

    def _init_workers(self, worker_class: QueueWorker, worker_init_args):
        return [worker_class.remote(self.task_queue, self.output_queue, *worker_args) #, max_restarts=-1
                for worker_args in worker_init_args]

    def _start_workers(self):
        for worker in self.workers:
            worker.start.remote()

    def _get_single_output(self):
        while True:
            try:
                output = self.output_queue.get()
                break
            except Empty:
                pass
        return output


# class MPQueueManager(QueueManager):
#     def _init_queue(self):
#         return MPQueue()
#
#     def _start_workers(self):
#         for worker in self.workers:
#             worker.start()  #???
#
#     def _get_single_output(self):
#         pass #????

# import ray
# @ray.remote
# class Duck(QueueWorker):
#     def _perform_task(self, task_input):
#         return f"{task_input}! QUACK! {task_input}!"