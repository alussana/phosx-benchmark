#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download Supplementary Table 2 from

Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature 613, 759–766 (2023).

<https://doi.org/10.1038/s41586-022-05575-3>
*/
process get_ser_thr_kinome_2023_suppl_table_2 {

    publishDir "${out_dir}", pattern: "datasets/41586_2022_5575_MOESM4_ESM.xlsx", mode: 'copy'

    output:
        path 'datasets/41586_2022_5575_MOESM4_ESM.xlsx'

    script:
    """
    mkdir -p datasets
    
    wget -P datasets/ "${params.url_ser_thr_kinome_2023_suppl_table_2}"
    """

}

/*
Parse output from get_ser_thr_kinome_2023_suppl_table_2()

Translate kinase name with gene synonym --> gene name id mapping

NOTE: PAK1 is only a synonym, not a gene name, according to the dictionary
      PKN1 is a gene name according to the dictionary
      But in the PSSM source files, both PAK1 and PKN1 are found and they
      correspond to two distinct PSSMs/kinases. Therefore the translation
      of PAK1 into PKN1 is not performed, otherwise two different PSSM would
      result having the same name, i.e. PKN1
      Same for PDHK1 --> PDK1 translation

Save HDF5 datasets with Ser/Thr kinases PSSMs based on the 
"ser_thr_all_norm_scaled_matrice" sheet of 
"input/kinome_2023_suppl_table_2.xlsx" spreadsheet from [1]

Also save each PSSM in text format (pssm/*.tsv) and plot the
sequence logos (pssm/logos/*.svg)

[1]
Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature 613, 759–766 (2023).
<https://doi.org/10.1038/s41586-022-05575-3>
*/
process get_ser_thr_kinase_pssm {

    publishDir "${out_dir}", pattern: "pssm/*tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "pssm/pssm_dict.h5", mode: 'copy'
    publishDir "${out_dir}", pattern: "pssm/logos/*.svg", mode: 'copy'

    input:
        path 'input/kinome_2023_suppl_table_2.xlsx'
        path 'input/gene_synonym_2_gene_name_dict.tsv'

    output:
        path 'pssm/pssm_dict.h5', emit: pssm_dict_h5
        path 'pssm/*.tsv'
        path 'pssm/logos/*.svg'

    script:
    """
    mkdir -p pssm/logos

    get_pssm_ser_thr_kinome_2023_suppl_table_2.py
    """

}