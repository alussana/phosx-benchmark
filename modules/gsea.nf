#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
perform a gsea prerank analysis given a ranked gene list and gene sets
using the gseapy package https://github.com/zqfang/GSEApy

filter the gene sets by allowing only sets (i.e. pathways, modules) that
include ${params.gseapy_min_size} <= genes <= ${params.gseapy_max_size}

run gsea prerank analysis only for modules/pathways that have an overlap with
the ranked gene list of 
${params.gseapy_min_overlap} <= n <= ${params.gseapy_max_overlap} genes

.META:
1 term: gene set name,
2 es: enrichment score,
3 nes: normalized enrichment score,
4 pval:  Nominal p-value (from the null distribution of the gene set,
5 fdr: FDR qvalue (adjusted False Discory Rate),
6 fwerp: Family wise error rate p-values,
7 tag %: Percent of gene set before running enrichment peak (ES),
8 gene %: Percent of gene list before running enrichment peak (ES),
9 lead_genes: leading edge genes (gene hits before running enrichment peak),
10 matched genes: genes matched to the data,
*/
process run_gsea {

    cpus "${params.gseapy_threads}"
    memory '10G'

    publishDir "${out_dir}",
            pattern: "GSEA/*",
            mode: 'copy'

    input:
        path 'input/clusters.tsv'
        tuple val(genes),
              file('input/genes.tsv')

    output:
        path "GSEA/${genes}/${genes}.csv", emit: csv
        path "GSEA/${genes}/*"

    script:
    """
    mkdir -p GSEA/

    cat input/genes.tsv | grep -v "gene" | grep -v "Gene" | awk 'NF' \
        > genes.rnk

    cat input/clusters.tsv \
        | awk -F "\\t" 'NF>${params.gseapy_min_size}' \
        | awk -F "\\t" 'NF<${params.gseapy_max_size}' \
        > clusters.tsv

    gseapy_prerank.py \
        clusters.tsv \
        genes.rnk \
        GSEA/${genes} \
        ${params.gseapy_threads} \
        ${params.gseapy_min_overlap} \
        ${params.gseapy_max_overlap} \
        ${params.gseapy_n_perm}

    cat GSEA/${genes}/gseapy.gene_set.prerank.report.csv \
        | cut -d ',' -f2- > GSEA/${genes}/${genes}.csv
    """

}