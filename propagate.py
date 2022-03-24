import json
from netprop.propagation import propagate_from_config, NETWORK_ORDERING_KEYWORD
from utils.new_file_loaders import CovToHumanMeta, NetpropCovToHumanLoader
from netprop.networks.randomization import randomize_network
from netprop.networks.loaders import NetpropNetwork
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
    #randomize_network(100, 15,
    #                  NetpropNetwork, [r"D:\data\networks\cov_to_human\metastudy.json"], {},
    #                  r"D:\data\networks\cov_to_human\rand_metastudy")

    propagate_from_config(r"D:\configurations\rna_randomized_prop_config.json", ordering={"network": 100})
    print("yay done <3")



if __name__ == "__main__":
    main()