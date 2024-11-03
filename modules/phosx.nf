#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
run PhosX given a seqrnk and a collection of PSSMs

NOTE: only S/T phosphosites are considered for now
*/
process run_phosx {

    cpus "${params.n_cores}"
    time '16h'
    memory '16G'

    publishDir "${out_dir}", pattern: "PhosX/${id}.tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "PhosX/${id}/*pdf", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.seqrnk')

    output:
        tuple val(id),
              file("PhosX/${id}.tsv"), emit: tsv
        path "PhosX/${id}/*pdf"

    script:
    """
    mkdir -p PhosX/${id}

    #phosx \
        input/input.seqrnk \
        -n ${params.phosx_n_perm} \
        -c ${params.n_cores} \
        -stk ${params.phosx_s_t_n_top_kinases} \
        -yk ${params.phosx_y_n_top_kinases} \
        -m ${params.phosx_min_n_hits} \
        --plot-figures \
        -d PhosX/${id} \
        > PhosX/${id}.tsv

    phosx \
        input/input.seqrnk \
        -n ${params.phosx_n_perm} \
        -c ${params.n_cores} \
        -stk ${params.phosx_s_t_n_top_kinases} \
        -yk ${params.phosx_y_n_top_kinases} \
        -mh ${params.phosx_min_n_hits} \
        -d PhosX/${id} \
        > PhosX/${id}.tsv

    touch PhosX/${id}/__phony__.pdf

    """

}


/*
translate all the words in the first field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate_phosx_output {

    publishDir "${out_dir}", pattern: "PhosX_tr/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        path "PhosX_tr/${id}.tsv"

    script:
    """
    mkdir -p PhosX_tr

    translator.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        3 \
        1 \
        1 \
        > PhosX_tr/${id}.tsv
    """

}