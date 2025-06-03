#!/usr/bin/env python3

import pandas as pd
import sys

from pairwise_comparisons.kinomics_per_kinase import make_kinase_activity_df, pairwise_comparison


def main():       
    input_list_phosx_txt = sys.argv[1]
    input_list_gsea_txt = sys.argv[2]
    input_list_kinex_txt = sys.argv[3]
    input_list_kstar_txt = sys.argv[4]
    input_list_phosxnouae_txt = sys.argv[5]
    metadata_tsv = sys.argv[6]
    kinase_activity_metric_str = sys.argv[7]
    out_prefix = sys.argv[8]

   
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


    # KSTAR kinase activity
    kinase_activity_kstar_df = make_kinase_activity_df(
        "KSTAR",
        input_list_kstar_txt,
        kinase_activity_metric_str,
        out_prefix,
    )
    

    # PhosX (no upstream activation evidence) kinase activity
    kinase_activity_phosxnouae_df = make_kinase_activity_df(
        "PhosX-NoUAE",
        input_list_phosxnouae_txt,
        kinase_activity_metric_str,
        out_prefix,
    )

    
    # Pairwise comparisons: PhosX vs <method>
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "Kinex",
        kinase_activity_kinex_df,
        metadata_df,
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
    pairwise_comparison(
        "PhosX",
        kinase_activity_phosx_df,
        "KSTAR",
        kinase_activity_kstar_df,
        metadata_df,
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


    # Pairwise comparisons: PhosX-NoUAE vs <method>
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "Kinex",
        kinase_activity_kinex_df,
        metadata_df,
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
    pairwise_comparison(
        "PhosX-NoUAE",
        kinase_activity_phosxnouae_df,
        "KSTAR",
        kinase_activity_kstar_df,
        metadata_df,
        out_prefix,
    )


if __name__ == "__main__":
    main()
