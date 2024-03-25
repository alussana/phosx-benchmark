#!/usr/bin/env python3

import sys
import pandas as pd

def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word, ][to_col - 1]
    if len(trs) == 0:
        return(None)
    else:
        return(trs.values[0])

def main():
    dict_file = sys.argv[1]
    
    word_dict = pd.read_csv(dict_file, sep='\t', header=None)
    # remove PAK1 --> PKN1 translation 
    # remove PDHK1 --> PDK1 translation
    # see get_ser_thr_kinase_pssm() process in modules/pssm.nf
    word_dict.loc[word_dict[0]=='PAK1'] = None
    word_dict.loc[word_dict[0]=='PDHK1'] = None
    
    scoring_matrix = pd.read_csv("https://zenodo.org/records/10201142/files/kinex_scoring_matrix_82k_sorted.csv.gzip?download=1", compression="gzip")
    
    # translate kinase synonyms into kinase name
    kinase_list = list(scoring_matrix)
    kinase_translated_list = []
    for kinase in kinase_list:
        translated_kinase = translateWord(kinase, word_dict, from_col=1, to_col=3)
        if translated_kinase == None:
            translated_kinase = kinase
        kinase_translated_list.append(translated_kinase)
    scoring_matrix.columns = kinase_translated_list
    
    print(scoring_matrix.to_csv(sep='\t', header=True, index=True))
    
    
if __name__ == '__main__':
    main()