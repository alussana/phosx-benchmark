#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    """
    dataset_csv = 'hernandez2017/benchmark_data.csv'
    metadata_csv = 'hernandez2017/benchmark_metadata.csv'
    rnk_out_prefix = 'datasets/hernandez2017/rnk/'
    metadata_out_tsv = 'datasets/hernandez2017/metadata.tsv'
    """
    dataset_csv = sys.argv[1]
    metadata_csv = sys.argv[2]
    rnk_out_prefix = sys.argv[3]
    metadata_out_tsv = sys.argv[4]

    metadata_df = pd.read_csv(metadata_csv, sep=',', index_col=None)
    metadata_df[['id', 'target','sign']].to_csv(metadata_out_tsv, sep='\t', index=False, header=False)
    
    dataset_df = pd.read_csv(dataset_csv, sep=',', index_col=0)
    
    experiment_id_list = list(dataset_df.columns)
    
    for id in experiment_id_list:
        rnk = dataset_df[id].dropna().sort_values(ascending=False)
        dataset_df[id].to_csv(f'{rnk_out_prefix}{id}.rnk', sep='\t', index=True, header=False)

if __name__ == '__main__':
    main()