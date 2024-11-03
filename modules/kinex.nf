#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
download kinex scoring matrix containing a priori PSSM score distributions for
the kineses in the human phosphoproteome

translate gene synonyms used for kinases into gene names
*/
process get_kinex_scoring_matrix {

    publishDir "${out_dir}", pattern: "datasets/Kinex/kinex_scoring_matrix_82k_translated.tsv", mode: 'copy'

    input:
        path 'input/gene_synonym_2_gene_name_dict.tsv'

    output:
        path "datasets/Kinex/kinex_scoring_matrix_82k_translated.tsv"

    script:
    """
    mkdir -p datasets/Kinex

    get_kinex_scoring_matrix.py \
        input/gene_synonym_2_gene_name_dict.tsv \
        > datasets/Kinex/kinex_scoring_matrix_82k_translated.tsv
    """

}


/*
run Kinex
*/
process run_kinex {

    publishDir "${out_dir}", pattern: "Kinex/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.seqrnk')
        path 'input/gene_synonym_2_gene_name_dict.tsv'
        path 'input/scoring_matrix.tsv'

    output:
        path "Kinex/${id}.tsv", emit: tsv

    script:
    """
    mkdir -p Kinex/${id}

    run_kinex.py \
        input/input.seqrnk \
        Kinex/${id}/${id}.pdf \
        ${params.kinex_fc_threshold} \
        input/gene_synonym_2_gene_name_dict.tsv \
        input/scoring_matrix.tsv \
        > Kinex/${id}.tsv
    """

}