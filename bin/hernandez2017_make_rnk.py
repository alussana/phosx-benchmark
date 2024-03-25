#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    """
    fold_tsv = 'input_with_possible_duplicates.rnk'
    """
    rnk_tsv = sys.argv[1]

    rnk = pd.read_csv(rnk_tsv, sep='\t', index_col=0, header=None)
    rnk.columns = ['score']
    genes = list(set(rnk.index))

    # lists with coordinated indexes for genes and fold changes
    g = []
    f = []

    # for each gene pick the most extreme fold change
    for gene in genes:
        fc = rnk.loc[gene, 'score']
        g.append(gene)
        fc_max = fc.max()
        fc_min = fc.min()
        if abs(fc_max) > abs(fc_min):
            f.append(fc_max)
        else:
            f.append(fc_min)

    out_df = pd.DataFrame({'gene':g, 'score':f})
    print(out_df.to_csv(sep='\t', index=False, header=False))


if __name__ == '__main__':
    main()