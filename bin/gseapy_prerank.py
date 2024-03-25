#!/usr/bin/env python3

import sys
import pandas as pd
import gseapy as gp

def main():
    """
    clusters_file = 'input/clusters.tsv'
    gene_list_file = 'genes.rnk'
    out_dir='.'
    threads=1
    min_size=0
    max_size=512
    permutations=1000
    """
    clusters_file = sys.argv[1]
    gene_list_file = sys.argv[2]
    out_dir = sys.argv[3]
    threads = int(sys.argv[4])
    min_size = int(sys.argv[5])
    max_size = int(sys.argv[6])
    permutations = int(sys.argv[7])
    
    try:
        min_size = int(sys.argv[5])
    except:
        min_size = None
        
    

    # read ranked gene list with attributes
    # .META:
    # 1     gene
    # 2     expr/logfc/...
    # the file is supposed to have a header and be tab-delimited
    rnk = pd.read_csv(gene_list_file, sep='\t', index_col=0, header=None)

    """if clusters_file == 'Reactome_2016':
        gene_sets = 'Reactome_2016'
    #elif clusters_file == 'KEGG_2021_Human':
    #    gene_sets = 'KEGG_2021_Human'
    else:
        # read gene sets
        # the file is supposed to have no header and a variable number of
        # tab-delimited fields
        # .META:
        # 1     set name
        # 2     gene1
        # 3     gene2
        # ...   ...
        gene_sets = {}
        with open(clusters_file, 'r') as clusters:
            for line in clusters:
                fields = line.strip().split('\t')
                set_name = fields.pop(0)
                gene_sets[set_name] = fields"""
    gene_sets = {}
    with open(clusters_file, 'r') as clusters:
        for line in clusters:
            fields = line.strip().split('\t')
            set_name = fields.pop(0)
            gene_sets[set_name] = fields

    # check agreement with constraints
    # if no gene set satisfyies them, skip the GSEA step, return empty file
    target_genes = set(list(rnk.index))
    contraints_satisfied = False
    for set_name in gene_sets.keys():
        module_set = set(gene_sets[set_name])
        intersection_len = len(module_set.intersection(target_genes))
        if intersection_len >= min_size and intersection_len <= max_size:
            contraints_satisfied = True
            break        

    if contraints_satisfied:
        pre_res = gp.prerank(
            rnk=rnk, # or a path on disk
            gene_sets=gene_sets, # supported: 'KEGG_2016', 'Reactome_2016', ...
            threads=threads,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutations, # reduce number to speed up testing
            outdir=out_dir, # if None, then don't write to disk
            seed=73,
            #weighted_score_type=0,
            ascending=False, # sorting order of rankings. Default: False.
            verbose=False, # see what's going on behind the scenes
        )
    else:
        pass

if __name__ == '__main__':
    main()