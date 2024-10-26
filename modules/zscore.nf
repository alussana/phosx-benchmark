#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
unzip precomputed kinase activity scores for Z-score on the hernandez2017
dataset; parse the text files for further processing

precomputed activities are from
Müller-Dott et al 
Comprehensive evaluation of phosphoproteomic-based kinase activity inference
<https://doi.org/10.1101/2024.06.27.601117>
*/
process zscore_hernandez2017 {

    publishDir "${out_dir}", pattern: "Z-score/untranslated/*.tsv", mode: 'copy'

    output:
        path "Z-score/untranslated/*.tsv"

    script:
    """
    mkdir -p Z-score/untranslated/

    unzip ${precomp_activity_scores}

    for file in \$(ls 03_activity_scores/zscore/hernandez/*tsv); do \
        idtsv=\$(basename \$file); id=\${idtsv::-4}; \
        cat \${file} | sed '1d' | cut -f2,4 \
        > Z-score/untranslated/\${id}.tsv; \
    done
    """

}


/*
unzip precomputed kinase activity scores for Z-score on the cptac
dataset; parse the text files for further processing

precomputed activities are from
Müller-Dott et al 
Comprehensive evaluation of phosphoproteomic-based kinase activity inference
<https://doi.org/10.1101/2024.06.27.601117>
*/
process zscore_cptac {

    publishDir "${out_dir}", pattern: "Z-score/untranslated/*.tsv", mode: 'copy'

    output:
        path "Z-score/untranslated/*.tsv"

    script:
    """
    mkdir -p Z-score/untranslated/

    unzip ${precomp_activity_scores}

    echo "brca ccrcc gbm hnscc lscc luad ucec" > cptac_types.txt

    for type in \$(cat cptac_types.txt); do \
        cat 03_activity_scores/zscore/cptac/\${type}.tsv | cut -f3 | sort | uniq > \${type}_cond.txt; \
        for cond in \$(cat \${type}_cond.txt); do \
            cat 03_activity_scores/zscore/cptac/\${type}.tsv | grep -w \${cond} > Z-score/untranslated/\${cond}.tsv; \
        done; \
    done
    """

}


/*
Translate kinase names from Z-score output
*/
process translate_zscore_output {

    cpus "${params.n_cores}"

    publishDir "${out_dir}", pattern: "Z-score/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file("input/${id}.tsv")
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path "Z-score/${id}.tsv", emit: tsv

    script:
    """
    mkdir -p Z-score

    translate_zscore_output.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        input/${id}.tsv \
        1 \
        > Z-score/${id}.tsv
    """

}