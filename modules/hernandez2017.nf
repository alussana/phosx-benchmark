#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Get phosphosite log fold changes in XXX peturbation experiments

[...]

Ref <https://doi.org/10.1093/bioinformatics/btx082>
*/
process parse_hernandez2017_dataset {

    publishDir "${out_dir}",
            pattern: 'datasets/hernandez2017/*.tsv',
            mode: 'copy'

    output:
        path 'datasets/hernandez2017/rnk/*.rnk', emit: rnk
        path 'datasets/hernandez2017/metadata.tsv', emit: metadata

    script:
    """
    mkdir -p datasets/hernandez2017/rnk

    unzip ${hernandez2017_zip}
 
    hernandez2017_parse_dataset.py \
        hernandez2017/benchmark_data.csv \
        hernandez2017/benchmark_metadata.csv \
        datasets/hernandez2017/rnk/ \
        datasets/hernandez2017/metadata.tsv

    sed -i 's/|/\\t/g' datasets/hernandez2017/rnk/*rnk

    for file in \$(ls datasets/hernandez2017/rnk/*.rnk); do \
        cp \${file} file.tmp; \
        cat file.tmp | awk '\$5!=""' > \${file}; \
        rm file.tmp; \
    done
    """

}


/*
translate all the words in the second field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate_hernandez2017_metadata {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv.gz'

    output:
        path 'translated_file.tsv'

    script:
    """
    translate_hernandez2017_metadata.py \
        input/dict.tsv.gz \
        input/file.tsv \
        > translated_file.tsv
    """

}


/*
translate all the words in the first field of input/file.tsv 
specified in the second tab-separated column of input/dict.tsv
with the corresponding word found in the first column

discard untranslated rows

.META:
1   UniProtAC 
2   residuePos
3   Gene Name
4   ENSG
5   ENSP
6   log2fc
*/
process hernandez2017_rnk_translate1col_w_dict {

    publishDir "${out_dir}",
            pattern: 'datasets/hernandez2017/uniprot_rnk/*.tsv',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        tuple val(id), file("datasets/hernandez2017/uniprot_rnk/${id}.tsv")

    script:
    """
    mkdir -p datasets/hernandez2017/uniprot_rnk

    cat input/file.tsv | sed 's/|/\\t/g' | awk '{print \$4"\\t"\$2"\\t"\$1"\\t"\$3"\\t"\$4"\\t"\$5}' > file.tsv

    translator.py \
        input/dict.tsv \
        file.tsv \
        2 \
        1 \
        1 \
        0 \
        > datasets/hernandez2017/uniprot_rnk/${id}.tsv
    """

}

/*
.META: (no header)
1   sequence    
2   log2fc
*/
process hernandez2017_make_seqrnk {

    publishDir "${out_dir}",
            pattern: "datasets/hernandez2017/seqrnk/${id}.seqrnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/hernandez2017/seqrnk/${id}.seqrnk")

    script:
    """
    mkdir -p datasets/hernandez2017/seqrnk

    hernandez2017_make_seqrnk.py \
        input/file.tsv \
        > datasets/hernandez2017/seqrnk/${id}.seqrnk

    touch datasets/hernandez2017/seqrnk/${id}.seqrnk
    """

}

/*
Make the rnk file

After translating ENSP to UniProt AC, uncommon duplicates may occur
For those phosphosites, the maximum absolute value is selected 

.META: (no header)
1   ID (UniProtAC_residuePos)
2   Score (fold change) 
*/
process hernandez2017_make_rnk {

    publishDir "${out_dir}",
            pattern: "datasets/hernandez2017/rnk/${id}.rnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/hernandez2017/rnk/${id}.rnk")

    script:
    """
    mkdir -p datasets/hernandez2017/rnk

    cat input/file.tsv \
        | awk '{print \$1"_"\$2"\\t"\$6}' \
        | sort -rgk2 \
        > rnk_with_possible_duplicates.tsv

    hernandez2017_make_rnk.py \
        rnk_with_possible_duplicates.tsv \
        | sed 's/-[0-9]*_/_/' \
        > datasets/hernandez2017/rnk/${id}.rnk
    """

}

/*
compute and plot distributions of auroc and precision at recall 50%, with the
same procedure as in Hernandez2017
<https://doi.org/10.1093/bioinformatics/btx082>
i.e. take positive kinase-condition pairs, create negative set by drawing the
same number of random  combinations of kinase-conditions pairs that are not in the positive set.
Compute AUROC and precision. Reapeat 60 times.

*/
process benchmark_phosx_hernandez2017 {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/union/*.pdf", mode: 'copy'
    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/intersection_s_t_y/*.pdf", mode: 'copy'
    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/intersection_s_t/*.pdf", mode: 'copy'
    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/union/*.tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/intersection_s_t_y/*.tsv", mode: 'copy'
    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/intersection_s_t/*.tsv", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/hernandez2017/union/*.pdf"
        path "kinase_activity_benchmark/hernandez2017/union/*.tsv", emit: tsv
        path "kinase_activity_benchmark/hernandez2017/intersection_s_t_y/*.pdf"
        path "kinase_activity_benchmark/hernandez2017/intersection_s_t/*.pdf"

    script:
    """
    CACHEBUST=0

    mkdir -p kinase_activity_benchmark/hernandez2017/union/
    mkdir -p kinase_activity_benchmark/hernandez2017/intersection_s_t_y/
    mkdir -p kinase_activity_benchmark/hernandez2017/intersection_s_t/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex

    for file in \$(ls input/phosx/); do \
        linkpath=\$(readlink -f input/phosx/"\$file"); \
        echo "\$linkpath" >> symlink_paths_phosx.txt; \
    done

    for file in \$(cat symlink_paths_phosx.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/phosx/\$name_dot_tsv; \
        echo data/phosx/\$name_dot_tsv >> paths_phosx.txt; \
    done

    cat paths_phosx.txt | sort -g > input_files_phosx.txt

    for file in \$(ls input/gsea/); do \
        linkpath=\$(readlink -f input/gsea/"\$file"); \
        echo "\$linkpath" >> symlink_paths_gsea.txt; \
    done

    for file in \$(cat symlink_paths_gsea.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/gsea/\$name_dot_tsv; \
        echo data/gsea/\$name_dot_tsv >> paths_gsea.txt; \
    done

    cat paths_gsea.txt | sort -g > input_files_gsea.txt

    for file in \$(ls input/kinex/); do \
        linkpath=\$(readlink -f input/kinex/"\$file"); \
        echo "\$linkpath" >> symlink_paths_kinex.txt; \
    done

    for file in \$(cat symlink_paths_kinex.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/kinex/\$name_dot_tsv; \
        echo data/kinex/\$name_dot_tsv >> paths_kinex.txt; \
    done

    cat paths_kinex.txt | sort -g > input_files_kinex.txt

    hernandez2017_phosx_kinase_activity_benchmark.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/hernandez2017/union/

    hernandez2017_phosx_kinase_activity_benchmark_intersection.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/hernandez2017/intersection_s_t/

    hernandez2017_phosx_kinase_activity_benchmark_tyr_intersection.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/hernandez2017/intersection_s_t_y/

    #hernandez2017_phosx_kinase_activity_benchmark_tyronly_intersection.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/hernandez2017/intersection_y/
    """

}