import json
from netprop.propagation import propagate_from_config, NETWORK_ORDERING_KEYWORD
from utils.new_file_loaders import CovToHumanMeta
from netprop.networks.loaders import NetpropNetworkModel
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
    propagate_from_config(r"D:\configurations\rna_prop_conf.json", ordering={NETWORK_ORDERING_KEYWORD: 1})
    #with open(r"D:\data\networks\try_with_covid.json", 'r') as handler:
     #   c = json.load(handler)

    #c = NetpropNetworkModel.load()

if __name__ == "__main__":
    main()