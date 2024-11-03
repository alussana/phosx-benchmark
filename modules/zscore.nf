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

    publishDir "${out_dir}", pattern: "Z-score/*.tsv", mode: 'copy'

    output:
        path "Z-score/*.tsv"

    script:
    """
    mkdir -p Z-score/

    unzip ${precomp_activity_scores}

    for file in \$(ls 03_activity_scores/zscore/hernandez/*tsv); do \
        idtsv=\$(basename \$file); id=\${idtsv::-4}; \
        cat \${file} | sed '1d' | cut -f2,4 \
        > Z-score/\${id}.tsv; \
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

    publishDir "${out_dir}", pattern: "Z-score/*.tsv", mode: 'copy'

    output:
        path "Z-score/*.tsv"

    script:
    """
    mkdir -p Z-score/

    unzip ${precomp_activity_scores}

    echo "brca ccrcc gbm hnscc lscc luad ucec" > cptac_types.txt

    for type in \$(cat cptac_types.txt); do \
        cat 03_activity_scores/zscore/cptac/\${type}.tsv | cut -f3 | sort | uniq > \${type}_cond.txt; \
        for cond in \$(cat \${type}_cond.txt); do \
            cat 03_activity_scores/zscore/cptac/\${type}.tsv \
            | grep -w \${cond} \
            | cut -f2,4 \
            > Z-score/\${cond}.tsv; \
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
process translate_zscore_output {

    publishDir "${out_dir}", pattern: "Z-score_tr/${id}.tsv", mode: 'copy'

    input:
        tuple val(id),
              file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        path "Z-score_tr/${id}.tsv"
         
    script:
    """
    mkdir -p Z-score_tr

    cat \
        <(echo -e "\\t${params.kinase_activity_metric}") \
        <(  translator.py \
            input/dict.tsv \
            input/file.tsv \
            1 \
            3 \
            1 \
            1 ) \
        > Z-score_tr/${id}.tsv
    """

}