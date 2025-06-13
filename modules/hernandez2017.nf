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
Add family and specificity annotations to the metadata, from the PhosX metadata HDF5
*/
process link_kinase_annotations_to_metadata {

    publishDir "${out_dir}", pattern: "datasets/hernandez2017/*tsv", mode: 'copy'

    input:
        path "input/metadata.tsv"
        path "input/gene_synomym2gene_name_dict.tsv"
        path "input/kinase_metadata.h5"

    output:
        path "datasets/hernandez2017/metadata_annotated.tsv"

    script:
    """
    mkdir -p datasets/hernandez2017
    
    link_kinase_annotations_to_metadata.py \
        input/metadata.tsv \
        input/gene_synomym2gene_name_dict.tsv \
        input/kinase_metadata.h5 \
        > datasets/hernandez2017/metadata_annotated.tsv
    """

}


/*
[...]
*/
process visualize_kinase_annotations {

    publishDir "${out_dir}", pattern: "datasets/hernandez2017/*pdf", mode: 'copy'

    input:
        path "input/metadata.tsv"

    output:
        path "datasets/hernandez2017/*.pdf"

    script:
    """
    mkdir -p datasets/hernandez2017

    visualize_kinase_annotations.py \
        input/metadata.tsv \
        datasets/hernandez2017/kinase_family_regulation_counts.pdf \
        datasets/hernandez2017/kinase_specificity_regulation_counts.pdf
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
.META: (no header)
1   UniProt AC
2   ResPos
3   sequence    
2   log2fc
*/
process hernandez2017_make_uniprot_seqrnk {

    publishDir "${out_dir}",
            pattern: "datasets/hernandez2017/uniprot_seqrnk/${id}.seqrnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/hernandez2017/uniprot_seqrnk/${id}.seqrnk")

    script:
    """
    mkdir -p datasets/hernandez2017/uniprot_seqrnk

    hernandez2017_make_uniprot_seqrnk.py \
        input/file.tsv \
        > datasets/hernandez2017/uniprot_seqrnk/${id}.seqrnk
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
[...]
*/
process benchmark_phosx_hernandez2017 {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/pairwise/*.pdf", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/kstar/*.tsv'
        path 'input/ptmsea/*.tsv'
        path 'input/zscore/*.tsv'
        path 'input/phosxnouae/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/hernandez2017/pairwise/*.pdf"

    script:
    """
    CACHEBUST=x

    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    mkdir -p kinase_activity_benchmark/hernandez2017/pairwise/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex
    mkdir -p data/kstar
    mkdir -p data/ptmsea
    mkdir -p data/zscore
    mkdir -p data/phosxnouae

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


    for file in \$(ls input/phosxnouae/); do \
        linkpath=\$(readlink -f input/phosxnouae/"\$file"); \
        echo "\$linkpath" >> symlink_paths_phosxnouae.txt; \
    done

    for file in \$(cat symlink_paths_phosxnouae.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/phosxnouae/\$name_dot_tsv; \
        echo data/phosxnouae/\$name_dot_tsv >> paths_phosxnouae.txt; \
    done

    cat paths_phosxnouae.txt | sort -g > input_files_phosxnouae.txt


    hernandez2017_phosx_kinase_activity_benchmark_pairwise.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input_files_kstar.txt \
        input_files_ptmsea.txt \
        input_files_zscore.txt \
        input_files_phosxnouae.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/hernandez2017/pairwise/
    """

}


/*
[...]
*/
process benchmark_phosx_per_kinase_hernandez2017 {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/hernandez2017/pairwise/*.pdf", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/kstar/*.tsv'
        path 'input/ptmsea/*.tsv'
        path 'input/zscore/*.tsv'
        path 'input/phosxnouae/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/hernandez2017/pairwise/*.pdf"

    script:
    """
    CACHEBUST=x
    
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    mkdir -p kinase_activity_benchmark/hernandez2017/pairwise/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex
    mkdir -p data/kstar
    mkdir -p data/ptmsea
    mkdir -p data/zscore
    mkdir -p data/phosxnouae

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


    for file in \$(ls input/phosxnouae/); do \
        linkpath=\$(readlink -f input/phosxnouae/"\$file"); \
        echo "\$linkpath" >> symlink_paths_phosxnouae.txt; \
    done

    for file in \$(cat symlink_paths_phosxnouae.txt); do \
        name_dot_tsv=\$(basename "\$file"); \
        cp \$file data/phosxnouae/\$name_dot_tsv; \
        echo data/phosxnouae/\$name_dot_tsv >> paths_phosxnouae.txt; \
    done

    cat paths_phosxnouae.txt | sort -g > input_files_phosxnouae.txt

    
    hernandez2017_phosx_kinase_activity_benchmark_pairwise_per_kinase_insights.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input_files_kstar.txt \
        input_files_ptmsea.txt \
        input_files_zscore.txt \
        input_files_phosxnouae.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/hernandez2017/pairwise/
    """

}