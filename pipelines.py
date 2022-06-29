from typing import Union, Callable
from functools import partial
from utils.utils import load_json, dump_json
from json.decoder import JSONDecodeError
from pathlib import Path

class AbstractPipeline:
    def __init__(self):
        self.attrs = dict()
        self._steps = []
        self.steps = self._steps

    @property
    def steps(self):
        for step in self._steps:
            yield step

    @steps.setter
    def steps(self, value):
        self._steps = value

    def execute(self):
        raise NotImplemented


class Pipeline(AbstractPipeline):
    def add_step(self, step: Union[AbstractPipeline, Callable]):
        self._steps.append(step)

    def add_steps(self, steps: list[Union[AbstractPipeline, Callable]]):
        for step in steps:
            self.add_step(step)

    def _execute_step(self, step):
        if isinstance(step, AbstractPipeline):
            step.attrs = self.attrs
            step.execute()
        elif callable(step):
            try:
                step(self.attrs)
            except TypeError as err:
                print("if a step is a callable it must accept exactly a single argument of type dict (self.attrs)")
                raise err
        else:
            raise TypeError("each step must either be a pipeline or a callable that modifies self.attrs")

    def execute(self):
        for step in self.steps:
            self._execute_step(step)


class NonRepeatingPipeline(Pipeline):
    _ATTARS__STATE_KEY = "attrs"
    _EXECUTION_STATE__STATE_KEY = "execution_state"

    def __init__(self, state_suffix: str, reset_state: bool=False):
        super().__init__()
        self._execution_state = []
        self._state_path = "_" + state_suffix + ".json"
        self._current_step_index = 0

        if reset_state:
            self.reset_state()

        if Path(self._state_path).exists():
            self.init_state()

    @property
    def steps(self):
        if not self._execution_state:
            self._execution_state = [False] * len(self._steps)

        self._current_step_index = 0
        for step in self._steps:
            yield step
            self._current_step_index += 1

    @steps.setter
    def steps(self, value):
        self._steps = value

    def _execute_step(self, step):
        if not self._execution_state[self._current_step_index]:
            super()._execute_step(step)
            self._execution_state[self._current_step_index] = True
        else:
            print(f"skipping step {self._current_step_index}")

    def execute(self):
        try:
            super().execute()
        finally:
            self.record_state()

    def init_state(self):
        try:
            state_dict = load_json(self._state_path)
        except (FileNotFoundError, JSONDecodeError) as file_read_exception:
            print(f"could  not load state from path {self._state_path}")
            raise file_read_exception

        self.attrs = state_dict.get(self._ATTARS__STATE_KEY, dict())
        self._execution_state = state_dict.get(self._EXECUTION_STATE__STATE_KEY, [])

    def record_state(self):
        serializable_attrs = dict()
        for attr_name, value in self.attrs.items():
            if isinstance(value, Path):
                value = str(value)
            serializable_attrs[attr_name] = value

        dump_json({self._ATTARS__STATE_KEY: serializable_attrs, self._EXECUTION_STATE__STATE_KEY: self._execution_state},
                  self._state_path)

    def reset_state(self):
        Path(self._state_path).unlink(missing_ok=True)
        self._execution_state = []
        self._current_step_index = 0
        self._steps = []
        self.attrs = dict()


p1 = NonRepeatingPipeline("p1_state")
p2 = NonRepeatingPipeline("p2_state", reset_state=True)


def add_key(key, attrs):
    attrs[key] = 1


def print_attrs(attrs):
    print(attrs)


# add_one = partial(add_key, "one")
# add_two = partial(add_key, "two")
# p1.add_step(add_one)
# p1.add_step(print_attrs)
# p2.add_step(add_two)
# p2.add_step(p1)
# p2.execute_pipeline()