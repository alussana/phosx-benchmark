#!/usr/bin/env python3

import sys
import pandas as pd
import requests
import warnings
from functools import lru_cache

@lru_cache
def get_fasta(uniprot_ac:str):
    url_prefix = 'https://rest.uniprot.org/uniprotkb/'
    r = requests.get(f'{url_prefix}{uniprot_ac}.fasta')
    fasta = r.content.decode()
    lines = fasta.split('\n')
    lines.pop(0)
    sequence = ''.join(lines)
    return sequence
    
def get_phosphosite(uniprot_ac:str, residue:str, pos:int, before:int, after:int):
    if before < 0 or after < 0:
        raise Exception(
            """
            'before' and 'after' have to be positive integers.
            """
        )
    sequence = get_fasta(uniprot_ac)
    seq_length = len(sequence)
    before_padding = '_' * before
    after_padding = '_' * after
    sequence = before_padding + sequence + after_padding
    try:
        phosphosite = sequence[pos-1:pos+before+after]
        if phosphosite[before] != residue:
            warnings.warn(
                f"""
                Phosphosite residue does not match the UniProt sequence position.
                uniprot_ac: {uniprot_ac}
                pos: {pos}
                expected residue: {residue}
                matched residue: {phosphosite[before]}
                """
            )
        else:
            return phosphosite
    except:
        warnings.warn(
            f"""
            Phosphosite position exceedes sequence length.
            pos: {pos}
            sequence length: {seq_length}
            """
        )
    
def make_uniprot_seqrnk_entry(series:pd.Series) -> pd.Series:
    uniprot_ac = series['UniProtAC']
    residue = series['Residue']
    pos = series['Position']
    respos = series['ResiduePosition']
    sequence = get_phosphosite(
        uniprot_ac=uniprot_ac,
        residue=residue,
        pos=pos,
        before=5,
        after=4
    )
    try:
        seq_list = [sequence[i] for i in range(len(sequence))]
        seq_list[5] = seq_list[5].lower()
        sequence = "".join(seq_list)
    except:
        pass
    score = series['Score']
    seqrnk_entry_series = pd.Series(
        {'UniProt_AC': uniprot_ac, 'ResPos': respos, 'Sequence': sequence, 'Score':score},
        name=series.name
    )
    return seqrnk_entry_series

   
def main():
    """
    input_tsv = 'input/file.tsv'
    output_seqrnk = 'datasets/hernandez2017/seqrnk/file.seqrnk'
    """
    input_tsv = sys.argv[1]

    df = pd.read_csv(input_tsv, sep='\t', header=None)
    df.columns = ['UniProtAC','ResiduePosition','Gene_name','ENSG','ENSP','Score']
    
    df['Residue'] = df.apply(lambda x: x['ResiduePosition'][0], axis=1)
    df['Position'] = df.apply(lambda x: x['ResiduePosition'][1:], axis=1)
    
    ac_res_pos_score_df = df.dropna()
    
    ac_res_pos_score_df['Position'] = ac_res_pos_score_df['Position'].astype(int)
    
    seqrnk = ac_res_pos_score_df.apply(make_uniprot_seqrnk_entry, axis=1)
    
    seqrnk = seqrnk.sort_values('Score', ascending=False)
    
    seqrnk = seqrnk.dropna()
        
    print(seqrnk.to_csv(sep='\t', index=False, header=False))
    

if __name__ == '__main__':
    main()