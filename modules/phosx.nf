#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
run PhosX given a seqrnk and a collection of PSSMs

NOTE: only S/T phosphosites are considered for now
*/
process run_phosx {

    cpus "${params.n_cores}"

    publishDir "${out_dir}", pattern: "PhosX/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.seqrnk')
        path 'input/pssm_dict.h5'
        path 'input/pssm_bg_scores.tsv.gz'

    output:
        path "PhosX/${id}.tsv", emit: tsv

    script:
    """
    mkdir -p PhosX/${id}

    phosx \
        input/input.seqrnk \
        -n 10000 \
        -c ${params.n_cores} \
        -k ${params.phosx_n_top_kinases} \
        -m ${params.phosx_min_n_hits} \
        > PhosX/${id}.tsv

    #phosx \
        input/input.seqrnk \
        -p input/pssm_dict.h5 \
        -n 10000 \
        -c ${params.n_cores} \
        -k ${params.phosx_n_top_kinases} \
        -m ${params.phosx_min_n_hits} \
        --plot-figures \
        -d PhosX/${id} \
        > PhosX/${id}.tsv

    """

}