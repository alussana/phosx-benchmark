#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
run PhosX given a seqrnk input
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
    CACHEBUST=42
    
    mkdir -p PhosX/${id}
    
    phosx \
        input/input.seqrnk \
        --n-permutations ${params.phosx_n_perm} \
        --n-proc ${params.n_cores} \
        --s-t-min-quantile ${params.phosx_s_t_min_quantile} \
        --y-min-quantile ${params.phosx_y_min_quantile} \
        --s-t-n-top-kinases ${params.phosx_s_t_n_top_kinases} \
        --y-n-top-kinases ${params.phosx_y_n_top_kinases} \
        --min-n-hits ${params.phosx_min_n_hits} \
        --a-loop-s-t-quantile-threshold ${params.phosx_st_qth} \
        --a-loop-y-quantile-threshold ${params.phosx_y_qth} \
        --upreg-redundancy-threshold ${params.phosx_upreg_redundancy_threshold} \
        --downreg-redundancy-threshold ${params.phosx_downreg_redundancy_threshold} \
        --decay-factor ${params.phosx_decay_factor} \
        --output-dir PhosX/${id} \
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
run PhosX (Ser/Thr kinases only) given a seqrnk input
*/
process run_phosx_serthr_only {

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
    CACHEBUST=42
    
    mkdir -p PhosX/${id}
    
    phosx \
        input/input.seqrnk \
        --n-permutations ${params.phosx_n_perm} \
        --n-proc ${params.n_cores} \
        --ser-thr-only \
        --s-t-min-quantile ${params.phosx_s_t_min_quantile} \
        --y-min-quantile ${params.phosx_y_min_quantile} \
        --s-t-n-top-kinases ${params.phosx_s_t_n_top_kinases} \
        --y-n-top-kinases ${params.phosx_y_n_top_kinases} \
        --min-n-hits ${params.phosx_min_n_hits} \
        --a-loop-s-t-quantile-threshold ${params.phosx_st_qth} \
        --a-loop-y-quantile-threshold ${params.phosx_y_qth} \
        --upreg-redundancy-threshold ${params.phosx_upreg_redundancy_threshold} \
        --downreg-redundancy-threshold ${params.phosx_downreg_redundancy_threshold} \
        --decay-factor ${params.phosx_decay_factor} \
        --output-dir PhosX/${id} \
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
run PhosX without upstream activation evidence given a seqrnk input
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
        --no-upstream-activation-evidence \
        --n-permutations ${params.phosx_n_perm} \
        --n-proc ${params.n_cores} \
        --s-t-min-quantile ${params.phosx_s_t_min_quantile} \
        --y-min-quantile ${params.phosx_y_min_quantile} \
        --s-t-n-top-kinases ${params.phosx_s_t_n_top_kinases} \
        --y-n-top-kinases ${params.phosx_y_n_top_kinases} \
        --min-n-hits ${params.phosx_min_n_hits} \
        --output-dir PhosXNouae/${id} \
        > PhosXNouae/${id}_untranslated.tsv

    translate_phosx_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        PhosXNouae/${id}_untranslated.tsv \
        1 \
        > PhosXNouae/${id}.tsv
        
    touch PhosXNouae/${id}/__phony__.pdf
    """

}


/*
run PhosX (Ser/Thr kinases only) without upstream activation evidence given
a seqrnk input
*/
process run_phosx_serthr_only_nouae {

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
        --no-upstream-activation-evidence \
        --n-permutations ${params.phosx_n_perm} \
        --ser-thr-only \
        --n-proc ${params.n_cores} \
        --s-t-min-quantile ${params.phosx_s_t_min_quantile} \
        --y-min-quantile ${params.phosx_y_min_quantile} \
        --s-t-n-top-kinases ${params.phosx_s_t_n_top_kinases} \
        --y-n-top-kinases ${params.phosx_y_n_top_kinases} \
        --min-n-hits ${params.phosx_min_n_hits} \
        --output-dir PhosXNouae/${id} \
        > PhosXNouae/${id}_untranslated.tsv

    translate_phosx_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        PhosXNouae/${id}_untranslated.tsv \
        1 \
        > PhosXNouae/${id}.tsv
        
    touch PhosXNouae/${id}/__phony__.pdf
    """

}