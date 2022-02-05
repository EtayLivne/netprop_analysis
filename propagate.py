import json
from netprop.propagation import propagate_from_config, NETWORK_ORDERING_KEYWORD
from utils.new_file_loaders import CoVToHumanMeta

def histogram_to_json(path, outpath):
    with open(path, 'r') as handler:
        lines = [line.replace('[', '').replace(']', '').replace(',', '').replace('\'', '').strip() for line in handler.readlines()]
    dictogram = {}
    for i in range(len(lines)):
        knockouts = [e.split('_')[0] for e in lines[i].split()]
        dictogram[i] = knockouts

    with open(outpath, 'w') as handler:
        json.dump(dictogram, handler, indent=4)






def main():
    # propagate_from_config(r"D:\configurations\temp.json")
    # histogram_to_json(r"D:\histogram_new_data.txt", r"D:\histogram_new_data.json")
    propagate_from_config(r"D:\configurations\test_config.json", ordering={NETWORK_ORDERING_KEYWORD: 1})

if __name__ == "__main__":
    main()