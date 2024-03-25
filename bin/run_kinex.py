#!/usr/bin/env python3

import sys
import pandas as pd
from kinex import Kinex
from kinex.data import get_pssm
#from kinex.data import get_groups
import matplotlib.pyplot as plt

def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word, ][to_col - 1]
    if len(trs) == 0:
        return(None)
    else:
        return(trs.values[0])

def main():
    seqrnk = sys.argv[1]
    out_pdf = sys.argv[2]
    fc_threshold = float(sys.argv[3])
    dict_file = sys.argv[4]
    scoring_matrix_tsv = sys.argv[5]
    
    word_dict = pd.read_csv(dict_file, sep='\t', header=None)
    # remove PAK1 --> PKN1 translation 
    # remove PDHK1 --> PDK1 translation
    # see get_ser_thr_kinase_pssm() process in modules/pssm.nf
    word_dict.loc[word_dict[0]=='PAK1'] = None
    word_dict.loc[word_dict[0]=='PDHK1'] = None
    
    scoring_matrix = pd.read_csv(scoring_matrix_tsv, sep='\t', index_col=0)
    
    # get and translate pssm
    pssm_table = get_pssm()
    # translate gene synonyms into gene names for kinases in pssm_table.index
    tr_index = []
    for i in range(len(pssm_table.index)):
        word = pssm_table.index[i]
        tr_word = translateWord(word, word_dict, from_col=1, to_col=3)
        if tr_word == None:
            tr_word = word
        tr_index.append(tr_word)
    pssm_table.index = tr_index
    
    pssm_table['kinase'] = pssm_table.index
    
    kinex = Kinex(scoring_matrix=scoring_matrix, pssm=pssm_table)
    
    input_sites = pd.read_csv(seqrnk, sep='\t', header=None)
    input_sites.columns = ['Sequence','Score']
    
    input_sites['Sequence'] = input_sites['Sequence'].apply(lambda x: f'{x[:5]}{x[5].lower()}*{x[6:]}')
    
    enrich = kinex.get_enrichment(input_sites, fc_threshold=fc_threshold, phospho_priming=False, favorability=True, method="max")
    
    out_df = enrich.enrichment_table.copy()
    out_df['dominant_direction'].loc[enrich.enrichment_table['dominant_direction']=='downregulated set'] = -1
    out_df['dominant_direction'].loc[enrich.enrichment_table['dominant_direction']=='upregulated set'] = 1
    out_df['Activity Score'] = out_df['dominant_direction'] * out_df['dominant_adjusted_p_value_log10_abs']
    
    print(out_df.to_csv(sep='\t', header=True, index=True))
    
    #fig = enrich.plot(use_adjusted_pval=True)
    #fig.write_image(out_pdf)
    
    
if __name__ == '__main__':
    main()