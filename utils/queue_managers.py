from time import sleep
from typing import Union
from typing import Callable
from functools import partial
from abc import ABC, abstractmethod
from threading import Lock, Thread, Event
from multiprocessing import Process, Manager
from multiprocessing import Queue as MPQueue
from queue import Empty as PythonEmptyQueueException

from ray.util.queue import Queue as RayQueue
from ray.util.queue import Empty as RayEmptyQueueException


class QueueWorker(ABC):
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

    @abstractmethod
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


def queue_worker_wrapper(worker_class: QueueWorker, task_queue: MPQueue, output_queue: MPQueue,
                         start_event_proxy: Event, *worker_init_args):
    obj = worker_class(task_queue, output_queue, *worker_init_args)
    # start_event_proxy.wait()
    obj.start()


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

    def _fill_task_queue(self, task_inputs):
        for task_input in task_inputs:
            self.task_queue.put(task_input)
        for _ in range(len(self.workers)):
            self.task_queue.put(self._stop_task_queue_message)

    def handle_inputs(self, task_inputs: list):
        self._fill_task_queue(task_inputs)
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
            except RayEmptyQueueException:
                sleep(0.05) # if no output wait a moment before checking again
        return output


class PythonQueueManager(QueueManager):
    def __init__(self, worker_class, num_workers: int, worker_init_args: Union[list, any] = None, use_ray: bool = True):
        self._sync_manager = Manager()
        self._shared_event = self._sync_manager.Event()
        self._shared_event.clear()
        super().__init__(worker_class, num_workers, worker_init_args=worker_init_args, use_ray=use_ray)

    def _init_queue(self):
        return MPQueue()

    def _init_workers(self, worker_class, worker_init_args):
        return \
            [
                Process(target=queue_worker_wrapper,
                        args=(worker_class, self.task_queue, self.output_queue, *worker_args))
                for worker_args in worker_init_args
            ]

    def _start_workers(self):
        for worker in self.workers:
            worker.start()

    def _get_single_output(self):
        while True:
            try:
                output = self.output_queue.get()
                break
            except PythonEmptyQueueException:
                sleep(0.05) # if no output wait a moment before checking again
        return output


def _activate_queue(lock: Lock, queue: QueueManager, inputs, accumulator_func: Callable):

    with lock:  # Will only happen after main thread releases this lock
        queue_outputs = queue.handle_inputs(inputs)
    return accumulator_func(queue_outputs) if accumulator_func else queue_outputs


# Currently assuming Ray implementation, may make it more generic in the future
class MultiQueueManager(ABC):
    def __init__(self, num_queues: int, worker_class: QueueWorker, workers_per_queue: int, accumulator_func=None,
                 worker_init_args: Union[list, any]=None):
        self._queues = [RayQueueManager(worker_class, workers_per_queue, worker_init_args=worker_init_args)
                        for _ in range(num_queues)]
        self._accumulator_func = accumulator_func
        self._availability_flags = [Lock() for _ in self._queues]

    def _get_available_queue(self):
        flag_num = 0
        while True:
            l = self._availability_flags[flag_num]
            acq = l.acquire(timeout=0.01)
            if acq:
                print(f"flag {flag_num} was available!")
                return flag_num
            flag_num = (flag_num + 1) % len(self._availability_flags)
            sleep(1)

    def _activate_single_queue(self, inputs, accumulator_func, queue_num) -> Thread:
        try:
            lock = self._availability_flags[queue_num]
            q = self._queues[queue_num]
            t = Thread(target=_activate_queue, args=(lock, q, inputs, accumulator_func)).start()

        finally:
            """
            the main program loops sequentially between here and _get_available_queue so we know no other thread is
            waiting on this. The thread will acquire lock again
            """
            lock.release()
            return t

    def handle_inputs(self, input_iterator) -> None:
        threads = []
        for inputs, accumulator_args in input_iterator:
            queue_num = self._get_available_queue()
            accumulator_func = partial(self._accumulator_func, **accumulator_args) #csv_path
            new_thread = self._activate_single_queue(inputs, accumulator_func, queue_num)
            threads.append(new_thread)


        sleep(6)
        for t in threads:
            t.join()
        print("joined all multi queue threads")

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