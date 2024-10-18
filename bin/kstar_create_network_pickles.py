#!/usr/bin/env python3

import pickle
import os
import sys
import pandas as pd

def create_network_pickles(in_net_dir, out_net_dir, phosphoTypes = ['Y','ST']):
    """
    Given network files declared in globals, create pickles of the kstar object that can then be quickly loaded in analysis
    Assumes that the Network structure has two folders Y and ST under the NETWORK_DIR global variable and that
    all .csv files in those directories should be loaded into a network pickle.
    """



    for phosphoType in phosphoTypes:
        network = {}
        if not os.path.isfile(f"{in_net_dir}/network_{phosphoType}.p"):
            directory = f"{in_net_dir}/{phosphoType}/INDIVIDUAL_NETWORKS/"
            #get all csv files in that directory
            for file in os.listdir(directory):
                if file.endswith(".tsv"):
                    #get the value of the network number
                    file_noext = file.strip(".tsv").split('_')
                    key_name = 'nkin'+str(file_noext[1])
                    #print("Debug: key name is %s"%(key_name))
                    network[key_name] = pd.read_csv(f"{directory}{file}", sep='\t')
            print("Loaded %d number of networks for phosphoType %s"%(len(network), phosphoType))
            pickle.dump(network, open(f"{out_net_dir}/kstar_network_{phosphoType}.p", "wb"))
            print(f"Saved pickle file at {out_net_dir}/kstar_network_{phosphoType}.p")
        else:
            print(f"{phosphoType} network pickle already generated")
            
            
if __name__ == "__main__":
    in_net_dir = sys.argv[1]
    out_net_dir = sys.argv[2]
    create_network_pickles(in_net_dir, out_net_dir)