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
        path 'input/dict.tsv.gz'

    output:
        path "datasets/Kinex/kinex_scoring_matrix_82k_translated.tsv"

    script:
    """
    mkdir -p datasets/Kinex

    get_kinex_scoring_matrix.py \
        input/dict.tsv.gz \
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
        path 'input/9606.protein.aliases.v12.0.txt.gz'
        path 'input/scoring_matrix.tsv'

    output:
        path "Kinex/${id}.tsv", emit: tsv

    script:
    """
    mkdir -p Kinex

    run_kinex.py \
        input/input.seqrnk \
        ${params.kinex_fc_threshold} \
        input/9606.protein.aliases.v12.0.txt.gz \
        input/scoring_matrix.tsv \
        Kinex/pdf/${id}.pdf \
        > Kinex/${id}.tsv
    """

}