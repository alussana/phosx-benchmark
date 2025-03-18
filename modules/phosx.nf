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

    publishDir "${out_dir}", pattern: "PhosX/${id}*.tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "PhosX/${id}/*pdf", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.seqrnk')
        path "input/9606.protein.aliases.v12.0.txt.gz"              

    output:
        path "PhosX/${id}.tsv", emit: tsv
        path "PhosX/${id}_untranslated.tsv"
        path "PhosX/${id}/*pdf"

    script:
    """
    mkdir -p PhosX/${id}

    CACHEBUST=42
    
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
        > PhosX/${id}_untranslated.tsv

    translate_phosx_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        PhosX/${id}_untranslated.tsv \
        1 \
        > PhosX/${id}.tsv
        
    touch PhosX/${id}/__phony__.pdf
    """

}