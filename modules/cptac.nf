#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
[...]
*/
process parse_cptac_dataset {

    publishDir "${out_dir}",
            pattern: 'datasets/cptac/*.tsv',
            mode: 'copy'

    output:
        path 'datasets/cptac/data.tsv', emit: data
        path 'datasets/cptac/metadata.tsv', emit: metadata

    script:
    """
    wget "${params.url_cptac_data}"
    wget "${params.url_cptac_metadata}"

    mkdir -p datasets/cptac

    cptac_parse_dataset.R
    """

}


/*
.META: datasets/cptac/samples/*.tsv
1   ENSG                ENSG00000048028.11
2   ENSP                ENSP00000003302.4
3   ResiduePosition     S1053
4   Sequence            TIRPNSPYDL
5   Score               0.812191717100202
*/
process split_cptac_samples {

    publishDir "${out_dir}",
            pattern: 'datasets/cptac/samples/*.tsv',
            mode: 'copy'

    input:
        path 'input/table.tsv'

    output:
        path 'datasets/cptac/samples/*.tsv'

    script:
    """
    mkdir -p datasets/cptac/samples

    cat input/table.tsv \
        | sed '1d' \
        | sed 's/|/\\t/g' \
        | cut -f2-5,7,8 \
        | awk '{ \
            split(\$4,a,""); \
            printf \$1"\\t"\$2"\\t"\$3"\\t"; \
            for (i=3; i<=12; i++) printf "%s", a[i]; \
            print "\\t"\$5"\\t"\$6 \
            }' \
        > parsed_table.tsv

    cat parsed_table.tsv \
        | cut -f5 \
        | sort \
        | uniq \
        > sample_list.txt
    
    for sample_id in \$(cat sample_list.txt); do \
        cat parsed_table.tsv \
            | grep -w \${sample_id} \
            | cut -f1-4,6 \
            > datasets/cptac/samples/\${sample_id}.tsv; \
    done
    """

}


/*
.META: datasets/cptac/rnk/*.rnk
1   ENSP_ResPos     ENSP00000003302.4_S1053
2   Score           0.812191717100202
*/
process make_cptac_rnk {

    publishDir "${out_dir}",
            pattern: 'datasets/cptac/rnk/*.rnk',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file('datasets/cptac/rnk/*.rnk')

    script:
    """
    mkdir -p datasets/cptac/rnk

    cat input/file.tsv \
        | awk '{print \$2"_"\$3"\\t"\$5}' \
        > datasets/cptac/rnk/${id}.rnk
    """

}


/*
.META: datasets/cptac/rnk/*.rnk
1   Sequence        QLGEESEERD 
2   Score           0.812191717100202
*/
process make_cptac_seqrnk {

    publishDir "${out_dir}",
            pattern: 'datasets/cptac/seqrnk/*.seqrnk',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file('datasets/cptac/seqrnk/*.seqrnk')

    script:
    """
    mkdir -p datasets/cptac/seqrnk

    cat input/file.tsv \
        | awk '{print \$4"\\t"\$5}' \
        > datasets/cptac/seqrnk/${id}.seqrnk
    """

}


/*
translate all the words in the first field of input/file.tsv 
specified in the second tab-separated column of input/dict.tsv
with the corresponding word found in the first column

discard untranslated rows

.META: datasets/cptac/rnk/*.rnk
1   UniprotAC_ResPos
2   Score
*/
process make_cptac_uniprot_rnk {

    publishDir "${out_dir}",
            pattern: 'datasets/cptac/uniprot_rnk/*.rnk',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.rnk')
        path 'input/dict.tsv'

    output:
        tuple val(id), file("datasets/cptac/uniprot_rnk/${id}.rnk")

    script:
    """
    mkdir -p datasets/cptac/uniprot_rnk

    cat input/file.rnk \
        | tr '_' '\\t' \
        | sed -r 's/\\.[0-9]+//' \
        > file.rnk

    translator.py \
        input/dict.tsv \
        file.rnk \
        2 \
        1 \
        1 \
        0 \
        | sed 's/\\t/_/' \
        | sed -r 's/-[0-9]+_/_/' \
        > datasets/cptac/uniprot_rnk/${id}.rnk
    """

}


/*
.META: (no header)
1   UniProt AC
2   ResPos
3   sequence    
2   log2fc
*/
process make_cptac_uniprot_seqrnk {

    publishDir "${out_dir}",
            pattern: "datasets/cptac/uniprot_seqrnk/${id}.seqrnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/cptac/uniprot_seqrnk/${id}.seqrnk")

    script:
    """
    mkdir -p datasets/cptac/uniprot_seqrnk

    cat input/file.tsv | tr '_' '\\t' > file.tsv

    cptac_make_uniprot_seqrnk.py \
        file.tsv \
        > datasets/cptac/uniprot_seqrnk/${id}.seqrnk
    """

}


// ===========
// TODO:

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
process cptac_rnk_translate1col_w_dict {

    publishDir "${out_dir}",
            pattern: 'datasets/cptac/uniprot_rnk/*.tsv',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        tuple val(id), file("datasets/cptac/uniprot_rnk/${id}.tsv")

    script:
    """
    mkdir -p datasets/cptac/uniprot_rnk

    cat input/file.tsv | sed 's/|/\\t/g' | awk '{print \$4"\\t"\$2"\\t"\$1"\\t"\$3"\\t"\$4"\\t"\$5}' > file.tsv

    translator.py \
        input/dict.tsv \
        file.tsv \
        2 \
        1 \
        1 \
        0 \
        > datasets/cptac/uniprot_rnk/${id}.tsv
    """

}

/*
.META: (no header)
1   sequence    
2   log2fc
*/
process cptac_make_seqrnk {

    publishDir "${out_dir}",
            pattern: "datasets/cptac/seqrnk/${id}.seqrnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/cptac/seqrnk/${id}.seqrnk")

    script:
    """
    mkdir -p datasets/cptac/seqrnk

    cptac_make_seqrnk.py \
        input/file.tsv \
        > datasets/cptac/seqrnk/${id}.seqrnk

    touch datasets/cptac/seqrnk/${id}.seqrnk
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
process cptac_make_rnk {

    publishDir "${out_dir}",
            pattern: "datasets/cptac/rnk/${id}.rnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/cptac/rnk/${id}.rnk")

    script:
    """
    mkdir -p datasets/cptac/rnk

    cat input/file.tsv \
        | awk '{print \$1"_"\$2"\\t"\$6}' \
        | sort -rgk2 \
        > rnk_with_possible_duplicates.tsv

    cptac_make_rnk.py \
        rnk_with_possible_duplicates.tsv \
        | sed 's/-[0-9]*_/_/' \
        > datasets/cptac/rnk/${id}.rnk
    """

}


/*
translate all the words in the second field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate_cptac_metadata {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv.gz'

    output:
        path 'translated_file.tsv'

    script:
    """
    translate_cptac_metadata.py \
        input/dict.tsv.gz \
        input/file.tsv \
        > translated_file.tsv
    """

}


/*
compute and plot distributions of auroc and precision at recall 50%, with the
same procedure as in cptac
<https://doi.org/10.1093/bioinformatics/btx082>
i.e. take positive kinase-condition pairs, create negative set by drawing the
same number of random  combinations of kinase-conditions pairs that are not in the positive set.
Compute AUROC and AUPR. Reapeat 100 times.
[...]
*/
process benchmark_phosx_cptac {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/cptac/pairwise/*.pdf", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/kstar/*.tsv'
        path 'input/ptmsea/*.tsv'
        path 'input/zscore/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/cptac/pairwise/*.pdf"

    script:
    """
    mkdir -p kinase_activity_benchmark/cptac/pairwise/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex
    mkdir -p data/kstar
    mkdir -p data/ptmsea
    mkdir -p data/zscore

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


    for file in \$(ls input/kstar/); do \
        linkpath=\$(readlink -f input/kstar/"\$file"); \
        echo "\$linkpath" >> symlink_paths_kstar.txt; \
    done

    for file in \$(cat symlink_paths_kstar.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/kstar/\$name_dot_tsv; \
        echo data/kstar/\$name_dot_tsv >> paths_kstar.txt; \
    done

    cat paths_kstar.txt | sort -g > input_files_kstar.txt


    for file in \$(ls input/ptmsea/); do \
        linkpath=\$(readlink -f input/ptmsea/"\$file"); \
        echo "\$linkpath" >> symlink_paths_ptmsea.txt; \
    done

    for file in \$(cat symlink_paths_ptmsea.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/ptmsea/\$name_dot_tsv; \
        echo data/ptmsea/\$name_dot_tsv >> paths_ptmsea.txt; \
    done

    cat paths_ptmsea.txt | sort -g > input_files_ptmsea.txt


    for file in \$(ls input/zscore/); do \
        linkpath=\$(readlink -f input/zscore/"\$file"); \
        echo "\$linkpath" >> symlink_paths_zscore.txt; \
    done

    for file in \$(cat symlink_paths_zscore.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/zscore/\$name_dot_tsv; \
        echo data/zscore/\$name_dot_tsv >> paths_zscore.txt; \
    done

    cat paths_zscore.txt | sort -g > input_files_zscore.txt


    cptac_phosx_kinase_activity_benchmark_pairwise.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input_files_kstar.txt \
        input_files_ptmsea.txt \
        input_files_zscore.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/cptac/pairwise/
    """

}