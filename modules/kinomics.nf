#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
[...]
*/
process parse_kinomics_tremetinib_dataset {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    output:
        path 'datasets/kinomics/tremetinib.seqrnk', emit: seqrnk
        path 'datasets/kinomics/tremetinib.uniprot_seqrnk', emit: uniprot_seqrnk
        path 'datasets/kinomics/tremetinib.rnk', emit: rnk

    script:
    """
    mkdir -p datasets/kinomics

    cat ${kinomics_tmt_mod_sites} \
        | sed 's/"//g' \
        | sed '1d' \
        | cut -f13,15,17 \
        | sed 's/^UniProtKB://' \
        | awk '{print \$1"\\t"\$2"\\t"toupper(\$3)}' \
        > phosphosites.tsv

    for phosphosite in \$(cat ${kinomics_tremetinib_seqrnk} | sed 's/\\t/@/'); do \
        seqrnk_seq=\$(echo \${phosphosite} | sed 's/@.*//'); \
        if [ \${#seqrnk_seq} -ne 10 ]; then \
            continue; \
        fi; \
        seq=\$(echo \${seqrnk_seq} | tr -d "_"); \
        score=\$(echo \${phosphosite} | sed 's/.*@//'); \
        psite=\$(cat phosphosites.tsv | grep "\$seq" || true); \
        if [[ ! -z "\${psite}" ]]; then \
            uniprot_ac=\$(echo \${psite} | cut -d " " -f1); \
            pos=\$(echo \${psite} | cut -d " " -f2); \
            seq_long=\$(echo \${psite} | cut -d " " -f3); \
            if [ \${#seq_long} -eq 14 ]; then \
                seqrnk_seq=\${seq_long:1:-3}; \
            fi; \
            before=\${seqrnk_seq:0:5}; \
            sixth_char=\${seqrnk_seq:5:1}; \
            after=\${seqrnk_seq:6}; \
            sixth_char_lower=\$(echo "\$sixth_char" | tr '[:upper:]' '[:lower:]'); \
            uniprot_seqrnk_seq="\${before}\${sixth_char_lower}\${after}"; \
            if [ \${#seqrnk_seq} -ne 10 ]; then \
                continue; \
            fi; \
            echo -e "\${seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/tremetinib.seqrnk; \
            echo -e "\${uniprot_ac}_\${sixth_char}\${pos}\\t\${score}" \
                >> datasets/kinomics/tremetinib.rnk; \
            echo -e "\${uniprot_ac}\\t\${sixth_char}\${pos}\\t\${uniprot_seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/tremetinib.uniprot_seqrnk; \
        fi; 
    done
    """
}


/*
[...]
*/
process parse_kinomics_tremetinib_metadata {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    output:
        path 'datasets/kinomics/tremetinib_metadata_untranslated.tsv', emit: metadata

    script:
    """
    mkdir -p datasets/kinomics

    parse_kinomics_metadata.py \
        "${kinomics_tremetinib_metadata}" \
        "${params.kinomics_statistic}" \
        "${params.kinomics_statistic_up_value}" \
        "${params.kinomics_statistic_down_value}" \
        tremetinib \
        | awk 'NF' \
        > datasets/kinomics/tremetinib_metadata_untranslated.tsv
    """
}


/*
[...]
*/
process translate_kinomics_tremetinib_metadata {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    input:
        path 'input/file.tsv'
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path 'datasets/kinomics/tremetinib_metadata.tsv', emit: metadata

    script:
    """
    mkdir -p datasets/kinomics

    translate_kinomics_metadata.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        input/file.tsv \
        | awk 'NF' \
        > datasets/kinomics/tremetinib_metadata.tsv
    """
}


/*
[...]
*/
process parse_kinomics_vemurafenib_dataset {

    publishDir "${out_dir}",
            pattern: 'datasets/kinomics/*',
            mode: 'copy'

    output:
        path 'datasets/kinomics/vemurafenib.seqrnk', emit: seqrnk
        path 'datasets/kinomics/vemurafenib.uniprot_seqrnk', emit: uniprot_seqrnk
        path 'datasets/kinomics/vemurafenib.rnk', emit: rnk

    script:
    """
    mkdir -p datasets/kinomics

    cat ${kinomics_tmt_mod_sites} \
        | sed 's/"//g' \
        | sed '1d' \
        | cut -f13,15,17 \
        | sed 's/^UniProtKB://' \
        | awk '{print \$1"\\t"\$2"\t"toupper(\$3)}' \
        > phosphosites.tsv

    for phosphosite in \$(cat ${kinomics_vemurafenib_seqrnk} | sed 's/\\t/@/'); do \
        seqrnk_seq=\$(echo \${phosphosite} | sed 's/@.*//'); \
        seq=\$(echo \${seqrnk_seq} | tr -d "_"); \
        score=\$(echo \${phosphosite} | sed 's/.*@//'); \
        psite=\$(cat phosphosites.tsv | grep \${seq} || true); \
        if [[ ! -z "\${psite}" ]]; then \
            uniprot_ac=\$(echo \${psite} | cut -d " " -f1); \
            pos=\$(echo \${psite} | cut -d " " -f2); \
            seq_long=\$(echo \${psite} | cut -d " " -f3); \
            if [ \${#seq_long} -eq 14 ]; then \
                seqrnk_seq=\${seq_long:1:-3}; \
            fi; \
            before=\${seqrnk_seq:0:5}; \
            sixth_char=\${seqrnk_seq:5:1}; \
            after=\${seqrnk_seq:6}; \
            sixth_char_lower=\$(echo "\$sixth_char" | tr '[:upper:]' '[:lower:]'); \
            uniprot_seqrnk_seq="\${before}\${sixth_char_lower}\${after}"; \
            if [ \${#seqrnk_seq} -ne 10 ]; then \
                continue; \
            fi; \
            echo -e "\${seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/vemurafenib.seqrnk; \
            echo -e "\${uniprot_ac}_\${sixth_char}\${pos}\\t\${score}" \
                >> datasets/kinomics/vemurafenib.rnk; \
            echo -e "\${uniprot_ac}\\t\${sixth_char}\${pos}\\t\${uniprot_seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/vemurafenib.uniprot_seqrnk; \
        fi; 
    done
    """
}


/*
[...]
*/
process parse_kinomics_vemurafenib_metadata {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    output:
        path 'datasets/kinomics/vemurafenib_metadata_untranslated.tsv', emit: metadata

    script:
    """
    mkdir -p datasets/kinomics

    parse_kinomics_metadata.py \
        "${kinomics_vemurafenib_metadata}" \
        "${params.kinomics_statistic}" \
        "${params.kinomics_statistic_up_value}" \
        "${params.kinomics_statistic_down_value}" \
        vemurafenib \
        | awk 'NF' \
        > datasets/kinomics/vemurafenib_metadata_untranslated.tsv
    """
}


/*
[...]
*/
process translate_kinomics_vemurafenib_metadata {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    input:
        path 'input/file.tsv'
        path 'input/9606.protein.aliases.v12.0.txt.gz'

    output:
        path 'datasets/kinomics/vemurafenib_metadata.tsv', emit: metadata

    script:
    """
    mkdir -p datasets/kinomics

    translate_kinomics_metadata.py \
        input/9606.protein.aliases.v12.0.txt.gz \
        input/file.tsv \
        > datasets/kinomics/vemurafenib_metadata.tsv
    """
}


/*
[...]
*/
process join_kinomics_metadata {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    input:
        path 'input/tremetinib_metadata.tsv'
        path 'input/vemurafenib_metadata.tsv'

    output:
        path 'datasets/kinomics/metadata.tsv'

    script:
    """
    mkdir -p datasets/kinomics
    
    cat input/tremetinib_metadata.tsv input/vemurafenib_metadata.tsv \
        | awk 'NF' \
        > datasets/kinomics/metadata.tsv
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
process kinomics_rnk_translate1col_w_dict {

    publishDir "${out_dir}",
            pattern: 'datasets/kinomics/uniprot_rnk/*.tsv',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        tuple val(id), file("datasets/kinomics/uniprot_rnk/${id}.tsv")

    script:
    """
    mkdir -p datasets/kinomics/uniprot_rnk

    cat input/file.tsv | sed 's/|/\\t/g' | awk '{print \$4"\\t"\$2"\\t"\$1"\\t"\$3"\\t"\$4"\\t"\$5}' > file.tsv

    translator.py \
        input/dict.tsv \
        file.tsv \
        2 \
        1 \
        1 \
        0 \
        > datasets/kinomics/uniprot_rnk/${id}.tsv
    """

}


/*
[...]
*/
process benchmark_phosx_kinomics {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/kinomics/pairwise/*.pdf", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/kstar/*.tsv'
        path 'input/phosxnouae/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/kinomics/pairwise/*.pdf"

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    mkdir -p kinase_activity_benchmark/kinomics/pairwise/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex
    mkdir -p data/kstar
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


    kinomics_phosx_kinase_activity_benchmark_pairwise.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input_files_kstar.txt \
        input_files_phosxnouae.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/kinomics/pairwise/
    """

}


/*
[...]
*/
process benchmark_phosx_per_kinase_kinomics {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/kinomics/pairwise/*.pdf", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/kstar/*.tsv'
        path 'input/phosxnouae/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/kinomics/pairwise/*.pdf"

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    mkdir -p kinase_activity_benchmark/kinomics/pairwise/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex
    mkdir -p data/kstar
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

    
    kinomics_phosx_kinase_activity_benchmark_pairwise_per_kinase_insights.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input_files_kstar.txt \
        input_files_phosxnouae.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/kinomics/pairwise/
    """

}