#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
unzip precomputed kinase activity scores for PTM-SEA on the hernandez2017
dataset; parse the text files for further processing

precomputed activities are from
Müller-Dott et al 
Comprehensive evaluation of phosphoproteomic-based kinase activity inference
<https://doi.org/10.1101/2024.06.27.601117>
*/
process ptmsea_hernandez2017 {

    publishDir "${out_dir}", pattern: "PTM-SEA/*.tsv", mode: 'copy'

    output:
        path "PTM-SEA/*.tsv"

    script:
    """
    mkdir -p PTM-SEA/

    unzip ${precomp_activity_scores}

    for file in \$(ls 03_activity_scores/ptmsea/hernandez/*tsv); do \
        idtsv=\$(basename \$file); id=\${idtsv::-4}; \
        cat \${file} | sed '1d' | cut -f2,4 \
        > PTM-SEA/\${id}.tsv; \
    done
    """

}


/*
unzip precomputed kinase activity scores for PTM-SEA on the CPTAC
dataset; parse the text files for further processing

precomputed activities are from
Müller-Dott et al 
Comprehensive evaluation of phosphoproteomic-based kinase activity inference
<https://doi.org/10.1101/2024.06.27.601117>
*/
process ptmsea_cptac {

    publishDir "${out_dir}", pattern: "PTM-SEA/*.tsv", mode: 'copy'

    output:
        path "PTM-SEA/*.tsv"

    script:
    """
    mkdir -p PTM-SEA/

    unzip ${precomp_activity_scores}

    echo "brca ccrcc gbm hnscc lscc luad ucec" > cptac_types.txt

    for type in \$(cat cptac_types.txt); do \
        cat 03_activity_scores/ptmsea/cptac/\${type}.tsv | cut -f3 | sort | uniq > \${type}_cond.txt; \
        for cond in \$(cat \${type}_cond.txt); do \
            if [ "\$type" = "brca" ]; then \
                cat 03_activity_scores/ptmsea/cptac/\${type}.tsv \
                | grep -w \${cond} \
                | cut -f2,4 \
                > PTM-SEA/X\${cond}.tsv; \
            else \
                cat 03_activity_scores/ptmsea/cptac/\${type}.tsv \
                | grep -w \${cond} \
                | cut -f2,4 \
                > PTM-SEA/\${cond}.tsv; \
            fi;
        done; \
    done
    """

}


/*
translate all the words in the first field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate_ptmsea_output {

    publishDir "${out_dir}", pattern: "PTM-SEA_tr/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        path "PTM-SEA_tr/${id}.tsv"

    script:
    """
    mkdir -p PTM-SEA_tr

    cat \
        <(echo -e "\\t${params.kinase_activity_metric}") \
        <(  translator.py \
            input/dict.tsv \
            input/file.tsv \
            1 \
            3 \
            1 \
            1 ) \
        > PTM-SEA_tr/${id}.tsv
    """

}