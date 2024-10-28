#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    metadata_tsv = sys.argv[1]
    kinomics_statistic = sys.argv[2]
    kinomics_statistic_up_value = float(sys.argv[3])
    kinomics_statistic_down_value = float(sys.argv[4])
    experiment = sys.argv[5]

    metadata_df = pd.read_csv(metadata_tsv, sep="\t", index_col=None)
    
    df = metadata_df[["Kinase Uniprot ID",kinomics_statistic]]
    
    df.loc[df[kinomics_statistic]>kinomics_statistic_up_value,"Regulation"] = 1
    
    df.loc[df[kinomics_statistic]<kinomics_statistic_down_value,"Regulation"] = -1
    
    df = df.dropna()
    
    df = df[["Kinase Uniprot ID", "Regulation"]]
    
    df['Experiment'] = experiment
    
    cols = ["Experiment", "Kinase Uniprot ID", "Regulation"]
    
    df = df[cols]
    
    print(df.to_csv(sep="\t", index=False, header=False))

if __name__ == '__main__':
    main()