#!/usr/bin/env python3

import sys
import requests
import pandas as pd

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
    except:
        raise Warning(
            f"""
            Phosphosite position exceedes sequence length.
            pos: {pos}
            sequence length: {seq_length}
            """
        )
    if phosphosite[before] != residue:
        raise Warning(
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

def main():
    phosphosite_list_tsv = sys.argv[1]
    
    pd.read_csv(phosphosite_list_tsv, sep='\t')
            
    

if __name__ == '__main__':
    main()