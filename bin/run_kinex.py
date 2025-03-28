#!/usr/bin/env python3

import sys
import pandas as pd
from kinex import Kinex
from kinex.data import get_pssm
#from kinex.data import get_groups
import matplotlib.pyplot as plt


alias_map_txt_gz = sys.argv[3]
alias_map_df = pd.read_csv(alias_map_txt_gz, sep="\t")

Ensembl_HGNC_symbol_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_HGNC_symbol", "alias"
].unique()

Ensembl_HGNC_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_HGNC", "alias"
].unique()

Ensembl_UniProt_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_UniProt", "alias"
].unique()


def find_translation(
    key,
):
    # The following mappings would collide if not for the ad-hoc flow control that follows:
    # Ser/Thr:
    # HGK     MAP4K4
    # NIK     MAP4K4
    # PDHK1   PDK1
    # PDK1    PDK1
    # Tyr:
    # EPHA3   EPHA3
    # ETK     EPHA3
    if key == "NIK":
        return "MAP3K14"
    elif key == "PDK1":
        return "PDPK1"
    elif key == "ETK":
        return "BMX"
    key_parts = key.split("_")
    key = key_parts[0]
    if (
        key in Ensembl_HGNC_symbol_list
        or key in Ensembl_HGNC_list
        or key in Ensembl_UniProt_list
    ):
        if len(key_parts) > 1:
            new_key = f"{key}_TYR"
        else:
            new_key = key
        return new_key
    key_alias_df = alias_map_df.loc[alias_map_df["alias"] == key]
    string_id_array = key_alias_df["#string_protein_id"].unique()
    if len(string_id_array) == 0:
        return key
    else:
        for string_id_str in string_id_array:
            try:
                map_df = alias_map_df.loc[
                    alias_map_df["#string_protein_id"] == string_id_str
                ]
                new_key = map_df.loc[map_df["source"] == "Ensembl_HGNC_symbol", "alias"]
                if len(new_key) == 0:
                    new_key = map_df.loc[map_df["source"] == "Ensembl_HGNC", "alias"]
                if len(new_key) == 0:
                    new_key = map_df.loc[map_df["source"] == "Ensembl_UniProt", "alias"]
                if len(new_key) == 0:
                    next
                new_key = new_key.values[0]
                if len(key_parts) > 1:
                    new_key = f"{new_key}_TYR"
                break
            except Exception as e:
                print(f"{key}: {e}")
                next
        return new_key


def main():
    seqrnk = sys.argv[1]
    fc_threshold = float(sys.argv[2])
    scoring_matrix_tsv = sys.argv[4]
    out_pdf = sys.argv[5]
    
    scoring_matrix = pd.read_csv(scoring_matrix_tsv, sep='\t', index_col=0)
    
    # get and translate pssm
    pssm_table = get_pssm()
    # translate gene synonyms into gene names for kinases in pssm_table.index
    tr_index = []
    for i in range(len(pssm_table.index)):
        word = pssm_table.index[i]
        tr_word = find_translation(word)
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
    
    # the following would work only with the original custom-named kinases
    #fig = enrich.plot(use_adjusted_pval=True)
    #fig.write_image(out_pdf)
    
    
if __name__ == '__main__':
    main()