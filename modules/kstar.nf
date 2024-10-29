#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
Make networks required by KSTAR

Ref:
Crowl, S., Jordan, B.T., Ahmed, H. et al. 
KSTAR: An algorithm to predict patient-specific kinase activities from phosphoproteomic data. 
Nat Commun 13, 4283 (2022). 
<https://doi.org/10.1038/s41467-022-32017-5>
*/
process make_kstar_networks {

    publishDir "${out_dir}", pattern: "datasets/kstar/*.p", mode: 'copy'

    output:
        path "datasets/kstar/kstar_network_ST.p", emit: st_net
        path "datasets/kstar/kstar_network_Y.p", emit: y_net

    script:
    """
    mkdir -p datasets/kstar

    wget -O KSTAR_graphs.tar.gz https://figshare.com/ndownloader/files/28768155 && \
        tar -xzf KSTAR_graphs.tar.gz # this creates NETWORKS/NetworKIN necessary for run_kstar.py

    kstar_create_network_pickles.py \
        NETWORKS/NetworKIN \
        datasets/kstar
    """

}


/*
Run KSTAR

Ref:
Crowl, S., Jordan, B.T., Ahmed, H. et al. 
KSTAR: An algorithm to predict patient-specific kinase activities from phosphoproteomic data. 
Nat Commun 13, 4283 (2022). 
<https://doi.org/10.1038/s41467-022-32017-5>

Compute activity score for each kinase:

The following values are considered for each kinase:

act_upreg = -log2(kstar mann-whitney FPR) in the upregulation tests
act_downreg = log2(kstar mann-whitney FPR) in the downregulation tests

The final kinase activity for a kinase is argmax_{(act_upreg, act_downreg)}{abs(act)}
*/
process run_kstar {

    cpus "${params.n_cores}"
    memory '32G'

    publishDir "${out_dir}", pattern: "KSTAR/untranslated/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.rnk')
        path 'input/kstar_network_ST.p'
        path 'input/kstar_network_Y.p'

    output:
        tuple val(id),
              file("KSTAR/untranslated/${id}.tsv")

    script:
    """
    mkdir -p KSTAR/output
    mkdir -p KSTAR/untranslated/

    run_kstar.py \
        input/input.rnk \
        KSTAR \
        ${id} \
        ${params.n_cores} \
        input \
        ${params.kstar_fc_threshold} \
        | grep -v network \
        > KSTAR/untranslated/${id}.tsv 2> kstar.err
    """

}


/*
Run KSTAR (exception tolerant version)

Ref:
Crowl, S., Jordan, B.T., Ahmed, H. et al. 
KSTAR: An algorithm to predict patient-specific kinase activities from phosphoproteomic data. 
Nat Commun 13, 4283 (2022). 
<https://doi.org/10.1038/s41467-022-32017-5>

Compute activity score for each kinase:

The following values are considered for each kinase:

act_upreg = -log2(kstar mann-whitney FPR) in the upregulation tests
act_downreg = log2(kstar mann-whitney FPR) in the downregulation tests

The final kinase activity for a kinase is argmax_{(act_upreg, act_downreg)}{abs(act)}
*/
process run_kstar_tolerant {

    cpus "${params.n_cores}"
    memory '32G'

    publishDir "${out_dir}", pattern: "KSTAR/untranslated/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.rnk')
        path 'input/kstar_network_ST.p'
        path 'input/kstar_network_Y.p'

    output:
        tuple val(id),
              file("KSTAR/untranslated/${id}.tsv")

    script:
    """
    mkdir -p KSTAR/output
    mkdir -p KSTAR/untranslated/

    run_kstar_tolerant.py \
        input/input.rnk \
        KSTAR \
        ${id} \
        ${params.n_cores} \
        input \
        ${params.kstar_fc_threshold} \
        | grep -v network \
        > KSTAR/untranslated/${id}.tsv 2> kstar.err
    """

}


/*
Translate kinase names from KSTAR output
*/
process translate_kstar_output {

    cpus "${params.n_cores}"

    publishDir "${out_dir}", pattern: "KSTAR/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/input.rnk')
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path "KSTAR/${id}.tsv", emit: tsv

    script:
    """
    mkdir -p KSTAR

    translate_kstar_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        input/input.rnk \
        1 \
        > KSTAR/${id}.tsv
    """

}