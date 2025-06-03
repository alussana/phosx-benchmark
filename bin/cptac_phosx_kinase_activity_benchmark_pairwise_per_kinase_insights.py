#!/usr/bin/env python3

import pandas as pd
import sys
import gc

from pairwise_comparisons.cptac_per_kinase import make_kinase_activity_df, pairwise_comparison


def main():
    input_list_phosx_txt = sys.argv[1]
    input_list_gsea_txt = sys.argv[2]
    input_list_kinex_txt = sys.argv[3]
    input_list_kstar_txt = sys.argv[4]
    input_list_ptmsea_txt = sys.argv[5]
    input_list_zscore_txt = sys.argv[6]
    input_list_phosxnouae_txt = sys.argv[7]
    metadata_tsv = sys.argv[8]
    kinase_activity_metric_str = sys.argv[9]
    out_prefix = sys.argv[10]
    """
    input_list_phosx_txt = 'input_files_phosx.txt'
    input_list_gsea_txt = 'input_files_gsea.txt'
    input_list_kinex_txt = 'input_files_kinex.txt'
    input_list_kstar_txt = 'input_files_kstar.txt'
    input_list_ptmsea_txt = 'input_files_ptmsea.txt'
    input_list_zscore_txt = 'input_files_zscore.txt'
    metadata_tsv = 'input/metadata.tsv'
    kinase_activity_metric_str = 'Activity Score'
    out_prefix = 'kinase_activity_benchmark/cptac/pairwise/'
    """
    
    
    # metadata - ground truth kinase regulation
    metadata_df = pd.read_csv(metadata_tsv, sep="\t", index_col=None, header=None)
    metadata_df.columns = ["Experiment", "Kinase", "Regulation"]
    metadata_df.index = metadata_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
    )
    

    gsea_kinase_activity_metric_str = "NES"


    # PhosX kinase activity
    kinase_activity_phosx_df = make_kinase_activity_df(
        "PhosX",
        input_list_phosx_txt,
        kinase_activity_metric_str,
        out_prefix,
    )


    # GSEApy kinase activity
    kinase_activity_gsea_df = make_kinase_activity_df(
        "GSEA",  
        input_list_gsea_txt,
        gsea_kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "GSEApy",
        kinase_activity_gsea_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_gsea_df
    gc.collect()


    # Kinex kinase activity
    kinase_activity_kinex_df = make_kinase_activity_df(
        "Kinex",
        input_list_kinex_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "Kinex",
        kinase_activity_kinex_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_kinex_df
    gc.collect()


    # KSTAR kinase activity
    kinase_activity_kstar_df = make_kinase_activity_df(
        "KSTAR",
        input_list_kstar_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "KSTAR",
        kinase_activity_kstar_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_kstar_df
    gc.collect()


    # PTM-SEA kinase activity
    kinase_activity_ptmsea_df = make_kinase_activity_df(
        "PTM-SEA",
        input_list_ptmsea_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "PTM-SEA",
        kinase_activity_ptmsea_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_ptmsea_df
    gc.collect()


    # Z-score kinase activity
    kinase_activity_zscore_df = make_kinase_activity_df(
        "Z-score",
        input_list_zscore_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "Z-score",
        kinase_activity_zscore_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_zscore_df
    gc.collect()
    

    # PhosX (no upstream activation evidence) kinase activity
    kinase_activity_phosxnouae_df = make_kinase_activity_df(
        "PhosX-NoUAE",
        input_list_phosxnouae_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_phosx_df
    gc.collect()


    # Kinex kinase activity
    kinase_activity_kinex_df = make_kinase_activity_df(
        "Kinex",
        input_list_kinex_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "Kinex",
        kinase_activity_kinex_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_kinex_df
    gc.collect()


    # GSEApy kinase activity
    kinase_activity_gsea_df = make_kinase_activity_df(
        "GSEA",  
        input_list_gsea_txt,
        gsea_kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "GSEApy",
        kinase_activity_gsea_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_gsea_df
    gc.collect()


    # KSTAR kinase activity
    kinase_activity_kstar_df = make_kinase_activity_df(
        "KSTAR",
        input_list_kstar_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "KSTAR",
        kinase_activity_kstar_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_kstar_df
    gc.collect()


    # PTM-SEA kinase activity
    kinase_activity_ptmsea_df = make_kinase_activity_df(
        "PTM-SEA",
        input_list_ptmsea_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "PTM-SEA",
        kinase_activity_ptmsea_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_ptmsea_df
    gc.collect()


    # Z-score kinase activity
    kinase_activity_zscore_df = make_kinase_activity_df(
        "Z-score",
        input_list_zscore_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "Z-score",
        kinase_activity_zscore_df,
        metadata_df,
        out_prefix,
    )
    del kinase_activity_zscore_df
    gc.collect()
    

if __name__ == "__main__":
    main()
