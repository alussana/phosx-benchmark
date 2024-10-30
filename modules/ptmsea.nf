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

    publishDir "${out_dir}", pattern: "PTM-SEA/untranslated/*.tsv", mode: 'copy'

    output:
        path "PTM-SEA/untranslated/*.tsv"

    script:
    """
    mkdir -p PTM-SEA/untranslated/

    unzip ${precomp_activity_scores}

    for file in \$(ls 03_activity_scores/ptmsea/hernandez/*tsv); do \
        idtsv=\$(basename \$file); id=\${idtsv::-4}; \
        cat \${file} | sed '1d' | cut -f2,4 \
        > PTM-SEA/untranslated/\${id}.tsv; \
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

    publishDir "${out_dir}", pattern: "PTM-SEA/untranslated/*.tsv", mode: 'copy'

    output:
        path "PTM-SEA/untranslated/*.tsv"

    script:
    """
    mkdir -p PTM-SEA/untranslated/

    unzip ${precomp_activity_scores}

    echo "brca ccrcc gbm hnscc lscc luad ucec" > cptac_types.txt

    for type in \$(cat cptac_types.txt); do \
        cat 03_activity_scores/ptmsea/cptac/\${type}.tsv | cut -f3 | sort | uniq > \${type}_cond.txt; \
        for cond in \$(cat \${type}_cond.txt); do \
            if [ "\$type" = "brca" ]; then \
                cat 03_activity_scores/ptmsea/cptac/\${type}.tsv \
                | grep -w \${cond} \
                | cut -f2,4 \
                > PTM-SEA/untranslated/X\${cond}.tsv; \
            else \
                cat 03_activity_scores/ptmsea/cptac/\${type}.tsv \
                | grep -w \${cond} \
                | cut -f2,4 \
                > PTM-SEA/untranslated/\${cond}.tsv; \
            fi;
        done; \
    done
    """

}


/*
Translate kinase names from PTM-SEA output
*/
process translate_ptmsea_output {

    cpus "${params.n_cores}"

    publishDir "${out_dir}", pattern: "PTM-SEA/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file("input/${id}.tsv")
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path "PTM-SEA/${id}.tsv", emit: tsv

    script:
    """
    mkdir -p PTM-SEA

    translate_ptmsea_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        input/${id}.tsv \
        1 \
        > PTM-SEA/${id}.tsv
    """

}