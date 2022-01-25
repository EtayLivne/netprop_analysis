import json


def load_json(file_path):
    with open(file_path, 'r') as handler:
        d = json.load(handler)
    return d


def dump_json(obj: dict, file_path: str):
    with open(file_path, 'w') as handler:
        json.dump(obj, handler, indent=4)
