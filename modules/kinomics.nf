#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
[...]
*/
process parse_kinomics_tremetinib_dataset {

    publishDir "${out_dir}", pattern: 'datasets/kinomics/*', mode: 'copy'

    output:
        path 'datasets/kinomics/tremetinib.seqrnk', emit: seqrnk
        path 'datasets/kinomics/tremetinib_uniprot.seqrnk', emit: uniprot_seqrnk
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
        seq=\$(echo \${seqrnk_seq} | tr -d "_"); \
        score=\$(echo \${phosphosite} | sed 's/.*@//'); \
        psite=\$(cat phosphosites.tsv | grep "\$seq" || true); \
        if [[ ! -z "\${psite}" ]]; then \
            uniprot_ac=\$(echo \${psite} | cut -d " " -f1); \
            pos=\$(echo \${psite} | cut -d " " -f2); \
            seq_long=\$(echo \${psite} | cut -d " " -f3); \
            if [ \$(echo "\${seq_long}" | wc -c) -eq 14 ]; then \
                seqrnk_seq=\${seq_long:1:-3}; \
            fi; \
            before=\${seqrnk_seq:0:5}; \
            sixth_char=\${seqrnk_seq:5:1}; \
            after=\${seqrnk_seq:6}; \
            sixth_char_lower=\$(echo "\$sixth_char" | tr '[:upper:]' '[:lower:]'); \
            uniprot_seqrnk_seq="\${before}\${sixth_char_lower}\${after}"; \
            echo -e "\${seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/tremetinib.seqrnk; \
            echo -e "\${uniprot_ac}_\${sixth_char}\${pos}\\t\${score}" \
                >> datasets/kinomics/tremetinib.rnk; \
            echo -e "\${uniprot_ac}\\t\${sixth_char}\${pos}\\t\${uniprot_seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/tremetinib_uniprot.seqrnk; \
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
        path 'datasets/kinomics/vemurafenib_uniprot.seqrnk', emit: uniprot_seqrnk
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
            if [ \$(echo "\${seq_long}" | wc -c) -eq 14 ]; then \
                seqrnk_seq=\${seq_long:1:-3}; \
            fi; \
            before=\${seqrnk_seq:0:5}; \
            sixth_char=\${seqrnk_seq:5:1}; \
            after=\${seqrnk_seq:6}; \
            sixth_char_lower=\$(echo "\$sixth_char" | tr '[:upper:]' '[:lower:]'); \
            uniprot_seqrnk_seq="\${before}\${sixth_char_lower}\${after}"; \
            echo -e "\${seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/vemurafenib.seqrnk; \
            echo -e "\${uniprot_ac}_\${sixth_char}\${pos}\\t\${score}" \
                >> datasets/kinomics/vemurafenib.rnk; \
            echo -e "\${uniprot_ac}\\t\${sixth_char}\${pos}\\t\${uniprot_seqrnk_seq}\\t\${score}" \
                >> datasets/kinomics/vemurafenib_uniprot.seqrnk; \
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



/*
.META: datasets/kinomics/samples/*.tsv
1   ENSG                ENSG00000048028.11
2   ENSP                ENSP00000003302.4
3   ResiduePosition     S1053
4   Sequence            TIRPNSPYDL
5   Score               0.812191717100202
*/
process split_kinomics_samples {

    publishDir "${out_dir}",
            pattern: 'datasets/kinomics/samples/*.tsv',
            mode: 'copy'

    input:
        path 'input/table.tsv'

    output:
        path 'datasets/kinomics/samples/*.tsv'

    script:
    """
    mkdir -p datasets/kinomics/samples

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
            > datasets/kinomics/samples/\${sample_id}.tsv; \
    done
    """

}


/*
.META: datasets/kinomics/rnk/*.rnk
1   ENSP_ResPos     ENSP00000003302.4_S1053
2   Score           0.812191717100202
*/
process make_kinomics_rnk {

    publishDir "${out_dir}",
            pattern: 'datasets/kinomics/rnk/*.rnk',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file('datasets/kinomics/rnk/*.rnk')

    script:
    """
    mkdir -p datasets/kinomics/rnk

    cat input/file.tsv \
        | awk '{print \$2"_"\$3"\\t"\$5}' \
        > datasets/kinomics/rnk/${id}.rnk
    """

}


/*
.META: datasets/kinomics/rnk/*.rnk
1   Sequence        QLGEESEERD 
2   Score           0.812191717100202
*/
process make_kinomics_seqrnk {

    publishDir "${out_dir}",
            pattern: 'datasets/kinomics/seqrnk/*.seqrnk',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file('datasets/kinomics/seqrnk/*.seqrnk')

    script:
    """
    mkdir -p datasets/kinomics/seqrnk

    cat input/file.tsv \
        | awk '{print \$4"\\t"\$5}' \
        > datasets/kinomics/seqrnk/${id}.seqrnk
    """

}


/*
translate all the words in the first field of input/file.tsv 
specified in the second tab-separated column of input/dict.tsv
with the corresponding word found in the first column

discard untranslated rows

.META: datasets/kinomics/rnk/*.rnk
1   UniprotAC_ResPos
2   Score
*/
process make_kinomics_uniprot_rnk {

    publishDir "${out_dir}",
            pattern: 'datasets/kinomics/uniprot_rnk/*.rnk',
            mode: 'copy'

    input:
        tuple val(id), file('input/file.rnk')
        path 'input/dict.tsv'

    output:
        tuple val(id), file("datasets/kinomics/uniprot_rnk/${id}.rnk")

    script:
    """
    mkdir -p datasets/kinomics/uniprot_rnk

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
        > datasets/kinomics/uniprot_rnk/${id}.rnk
    """

}


/*
.META: (no header)
1   UniProt AC
2   ResPos
3   sequence    
2   log2fc
*/
process make_kinomics_uniprot_seqrnk {

    publishDir "${out_dir}",
            pattern: "datasets/kinomics/uniprot_seqrnk/${id}.seqrnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/kinomics/uniprot_seqrnk/${id}.seqrnk")

    script:
    """
    mkdir -p datasets/kinomics/uniprot_seqrnk

    cat input/file.tsv | tr '_' '\\t' > file.tsv

    kinomics_make_uniprot_seqrnk.py \
        file.tsv \
        > datasets/kinomics/uniprot_seqrnk/${id}.seqrnk
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
.META: (no header)
1   sequence    
2   log2fc
*/
process kinomics_make_seqrnk {

    publishDir "${out_dir}",
            pattern: "datasets/kinomics/seqrnk/${id}.seqrnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/kinomics/seqrnk/${id}.seqrnk")

    script:
    """
    mkdir -p datasets/kinomics/seqrnk

    kinomics_make_seqrnk.py \
        input/file.tsv \
        > datasets/kinomics/seqrnk/${id}.seqrnk

    touch datasets/kinomics/seqrnk/${id}.seqrnk
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
process kinomics_make_rnk {

    publishDir "${out_dir}",
            pattern: "datasets/kinomics/rnk/${id}.rnk",
            mode: 'copy'

    input:
        tuple val(id), file('input/file.tsv')

    output:
        tuple val(id), file("datasets/kinomics/rnk/${id}.rnk")

    script:
    """
    mkdir -p datasets/kinomics/rnk

    cat input/file.tsv \
        | awk '{print \$1"_"\$2"\\t"\$6}' \
        | sort -rgk2 \
        > rnk_with_possible_duplicates.tsv

    kinomics_make_rnk.py \
        rnk_with_possible_duplicates.tsv \
        | sed 's/-[0-9]*_/_/' \
        > datasets/kinomics/rnk/${id}.rnk
    """

}


/*
translate all the words in the second field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate_kinomics_metadata {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv.gz'

    output:
        path 'translated_file.tsv'

    script:
    """
    translate_kinomics_metadata.py \
        input/dict.tsv.gz \
        input/file.tsv \
        > translated_file.tsv
    """

}


/*
compute and plot distributions of auroc and precision at recall 50%, with the
same procedure as in kinomics
<https://doi.org/10.1093/bioinformatics/btx082>
i.e. take positive kinase-condition pairs, create negative set by drawing the
same number of random  combinations of kinase-conditions pairs that are not in the positive set.
Compute AUROC and AUPR. Reapeat 100 times.
[...]
*/
process benchmark_phosx_kinomics {

    publishDir "${out_dir}", pattern: "kinase_activity_benchmark/kinomics/pairwise/*.pdf", mode: 'copy'

    input:
        path 'input/phosx/*.tsv'
        path 'input/gsea/*.csv'
        path 'input/kinex/*.tsv'
        path 'input/kstar/*.tsv'
        path 'input/metadata.tsv'

    output:
        path "kinase_activity_benchmark/kinomics/pairwise/*.pdf"

    script:
    """
    mkdir -p kinase_activity_benchmark/kinomics/pairwise/
    mkdir -p data/phosx
    mkdir -p data/gsea
    mkdir -p data/kinex
    mkdir -p data/kstar

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


    kinomics_phosx_kinase_activity_benchmark_pairwise.py \
        input_files_phosx.txt \
        input_files_gsea.txt \
        input_files_kinex.txt \
        input_files_kstar.txt \
        input/metadata.tsv \
        "${params.kinase_activity_metric}" \
        kinase_activity_benchmark/kinomics/pairwise/
    """

}