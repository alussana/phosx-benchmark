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
    CACHEBUST=discipline
    
    mkdir -p PhosX/${id}
    
    phosx \
        input/input.seqrnk \
        -n ${params.phosx_n_perm} \
        -c ${params.n_cores} \
        -mp ${params.phosx_min_quantile} \
        -stk ${params.phosx_s_t_n_top_kinases} \
        -yk ${params.phosx_y_n_top_kinases} \
        -mh ${params.phosx_min_n_hits} \
        -astqth ${params.phosx_st_qth} \
        -ayqth ${params.phosx_y_qth} \
        -urt ${params.phosx_upreg_redundancy_threshold} \
        -drt ${params.phosx_downreg_redundancy_threshold} \
        -df1 ${params.phosx_decay_factor} \
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


/*
run PhosX given a seqrnk and a collection of PSSMs

NOTE: only S/T phosphosites are considered for now
*/
process run_phosx_nouae {

    cpus "${params.n_cores}"
    time '16h'
    memory '16G'

    publishDir "${out_dir}", pattern: "PhosXNouae/${id}*.tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "PhosXNouae/${id}/*pdf", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.seqrnk')
        path "input/9606.protein.aliases.v12.0.txt.gz"              

    output:
        path "PhosXNouae/${id}.tsv", emit: tsv
        path "PhosXNouae/${id}_untranslated.tsv"
        path "PhosXNouae/${id}/*pdf"

    script:
    """
    mkdir -p PhosXNouae/${id}
    
    phosx \
        input/input.seqrnk \
        -no-uae \
        -n ${params.phosx_n_perm} \
        -c ${params.n_cores} \
        -mp ${params.phosx_min_quantile} \
        -stk ${params.phosx_s_t_n_top_kinases} \
        -yk ${params.phosx_y_n_top_kinases} \
        -mh ${params.phosx_min_n_hits} \
        -d PhosXNouae/${id} \
        > PhosXNouae/${id}_untranslated.tsv

    translate_phosx_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        PhosXNouae/${id}_untranslated.tsv \
        1 \
        > PhosXNouae/${id}.tsv
        
    touch PhosXNouae/${id}/__phony__.pdf
    """

}