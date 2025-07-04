#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import gc 
from functools import partial
import dask.dataframe as dd


plt.rcParams['axes.titlesize'] = 8        
plt.rcParams['axes.labelsize'] = 5        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6 


red_shade = "#d62728"
grey_shade = "#7f7f7f"
black = "#000000"


plt.rcParams['axes.prop_cycle'] = plt.cycler(color=[red_shade, black])


def quantile_normalize(df):
    """
    Memory-efficient quantile normalization
    """
    # Process in chunks if dataframe is large
    if df.memory_usage().sum() > 1e9:  # If dataframe > 1GB
        # Save to temporary file
        temp_file = "temp_df_for_qnorm.parquet"
        df.to_parquet(temp_file)
        
        # Process in chunks
        chunk_size = max(1000, len(df) // 10) 
        reader = pd.read_parquet(temp_file, chunksize=chunk_size)
        
        # Initialize result dataframe
        result = pd.DataFrame()
        
        # Collect ranks
        all_ranks = []
        for chunk in reader:
            ranks = chunk.rank(method="first")
            all_ranks.append(ranks)
        
        ranks_df = pd.concat(all_ranks)
        del all_ranks
        gc.collect()
        
        # Determine means
        values = df.values.flatten()
        ranks = ranks_df.values.flatten()
        sorted_idx = np.argsort(ranks)
        sorted_values = values[sorted_idx]
        
        # Compute means for ties
        result = pd.DataFrame(index=df.index, columns=df.columns)
        for i, col in enumerate(df.columns):
            result[col] = np.interp(ranks_df[col], np.sort(ranks), sorted_values)
        
        # Clean up
        try:
            os.remove(temp_file)
        except:
            pass
        
        return result
    else:
        # for smaller dataframes
        ranked = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
        quantiles = df.rank(method="min").stack().astype(int).map(ranked).unstack()
        return quantiles.sort_index()

def scale_01(df):
    """
    Memory-efficient 0-1 scaling
    """
    # Get min and max without creating additional full dataframes
    mins = df.min()
    maxs = df.max()
    
    # Process column by column to reduce memory
    for col in df.columns:
        min_val = mins[col]
        max_val = maxs[col] - min_val
        # In-place operations to avoid creating additional dataframes
        df[col] = (df[col] - min_val) / max_val if max_val > 0 else 0
    
    return df

def process_top_bottom_percentile(kinase_activity_df, regulation_df, metadata_experiments_set, percentile=0.1, top=True):
    """
    Generalized function for top or bottom percentile processing with improved memory usage
    """
    results = []
    
    # Process one experiment at a time to save memory
    for experiment in metadata_experiments_set:
        # Filter the dataframes to reduce memory usage
        reg_kinases = regulation_df.loc[regulation_df["Experiment"] == experiment, "Kinase"].tolist()
        
        if not reg_kinases:
            continue
            
        # Get only the activities for this experiment
        exp_activities = kinase_activity_df.loc[kinase_activity_df["Experiment"] == experiment].copy()
        
        if exp_activities.empty:
            continue
            
        # Calculate unique values and ranks only for this experiment
        uniq_vals = exp_activities["Kinase activity change"].sort_values(ascending=True).unique()
        ranks = np.arange(len(uniq_vals))
        vals2ranks = dict(zip(uniq_vals, ranks))
        
        # Apply ranking in-place
        exp_activities["Rank"] = exp_activities["Kinase activity change"].map(vals2ranks)
        n_ranks = exp_activities["Rank"].max()
        
        # Process kinases for this experiment
        for kinase in reg_kinases:
            kinase_rank = exp_activities.loc[exp_activities["Kinase"] == kinase, "Rank"]
            
            if kinase_rank.empty:
                continue
                
            # Evaluate if this kinase is in the top/bottom percentile
            if top:
                is_match = kinase_rank.values[0] / n_ranks > (1 - percentile)
            else:
                is_match = kinase_rank.values[0] / n_ranks < percentile
                
            if is_match:
                results.append([kinase_rank.index[0], 1])
            else:
                results.append([kinase_rank.index[0], 0])
    
    # Create dataframe from results only once
    if results:
        df = pd.DataFrame(results, columns=["Index", "Rank"])
        df.set_index("Index", inplace=True)
        return df
    else:
        return pd.DataFrame(columns=["Rank"])

top10percentKinases = partial(process_top_bottom_percentile, percentile=0.1, top=True)
bottom10percentKinases = partial(process_top_bottom_percentile, percentile=0.1, top=False)


def make_kinase_activity_df(
    method_name:str,
    input_list_txt:str,
    kinase_activity_metric:str,
    out_prefix:str,
):
    # build dictionary (experiment id --> kinase activity table path)
    data_path_dict = {}
    with open(input_list_txt, "r") as input_list_fh:
        for line in input_list_fh:
            data_path_str = line.strip()
            data_id = os.path.basename(data_path_str)[:-4]
            data_path_dict[data_id] = data_path_str
    # build dataframe of kinase activities
    kinase_activity_df = pd.DataFrame()
    for data_id in data_path_dict.keys():
        df = pd.read_csv(data_path_dict[data_id], sep="\t", index_col=0)
        if len(df.columns) == 0:
            df = pd.read_csv(data_path_dict[data_id], sep=",", index_col=0)
        series = df[kinase_activity_metric]
        series.fillna(0, inplace=True)
        series.name = data_id
        kinase_activity_df = kinase_activity_df.join(series, how="outer")
    # Get columns that can be successfully converted to float
    convertible_columns = []
    for col in kinase_activity_df.columns:
        try:
            pd.to_numeric(kinase_activity_df[col], errors='raise')
            convertible_columns.append(col)
        except (ValueError, TypeError):
            pass
    kinase_activity_df = kinase_activity_df[convertible_columns]
    # normalise and scale
    kinase_activity_df = quantile_normalize(kinase_activity_df)
    kinase_activity_df = scale_01(kinase_activity_df)
    # melt
    kinase_activity_df = pd.melt(
        kinase_activity_df.assign(index=kinase_activity_df.index),
        id_vars=["index"],
    )
    kinase_activity_df.columns = [
        "Kinase",
        "Experiment",
        "Kinase activity change",
    ]
    # remove NAs due to missing phosphosite sequence
    kinase_activity_df = kinase_activity_df.dropna()
    # make index
    kinase_activity_df.index = kinase_activity_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
    )
    kinase_activity_df.to_csv(
        f"{out_prefix}cptac_kinase_activity_{method_name}.tsv",
        header=True,
        index=True,
        sep="\t",
    )
    return kinase_activity_df


def pairwise_comparison(
    method_1_name,
    kinase_activity_method_1_df,
    method_2_name,
    kinase_activity_method_2_df,
    metadata_df,
    out_prefix,
):
    """
    Memory-optimized pairwise comparison
    """
    # Calculate intersection without creating new large dataframes
    intersection_index = set(kinase_activity_method_1_df.index).intersection(
        set(kinase_activity_method_2_df.index)
    )
    
    # Filter dataframes using the intersection
    kinase_activity_method_1_df = kinase_activity_method_1_df.loc[list(intersection_index)]
    kinase_activity_method_2_df = kinase_activity_method_2_df.loc[list(intersection_index)]
    
    # Extract regulation dataframes
    upregulation_df = metadata_df.loc[metadata_df["Regulation"] == 1]
    downregulation_df = metadata_df.loc[metadata_df["Regulation"] == -1]
    
    # Save regulation dataframes
    upregulation_df.to_csv(
        f"{out_prefix}cptac_{method_1_name}_{method_2_name}_upregulated.tsv",
        header=True,
        index=True,
        sep="\t"
    )
    downregulation_df.to_csv(
        f"{out_prefix}cptac_{method_1_name}_{method_2_name}_downregulated.tsv",
        header=True,
        index=True,
        sep="\t",
    )
    
    metadata_kinases_set = set(metadata_df["Kinase"].unique())
    metadata_experiments_set = set(metadata_df["Experiment"].unique())
    
    method_1_top_percentKinases_df = top10percentKinases(
        kinase_activity_method_1_df, upregulation_df, metadata_experiments_set
    )
    method_2_top_percentKinases_df = top10percentKinases(
        kinase_activity_method_2_df, upregulation_df, metadata_experiments_set
    )
    method_1_bottom_percentKinases_df = bottom10percentKinases(
        kinase_activity_method_1_df, downregulation_df, metadata_experiments_set
    )
    method_2_bottom_percentKinases_df = bottom10percentKinases(
        kinase_activity_method_2_df, downregulation_df, metadata_experiments_set
    )
    
    method_1_joined_percentKinases_df = pd.concat([
        method_1_top_percentKinases_df, 
        method_1_bottom_percentKinases_df
    ])
    method_2_joined_percentKinases_df = pd.concat([
        method_2_top_percentKinases_df, 
        method_2_bottom_percentKinases_df
    ])
    
    # Extract kinase names and aggregate results
    method_1_joined_percentKinases_df["Kinase"] = method_1_joined_percentKinases_df.index.str.split("__").str[1]
    method_2_joined_percentKinases_df["Kinase"] = method_2_joined_percentKinases_df.index.str.split("__").str[1]
    
    # Group by kinase and calculate sums
    method_1_hit_count_df = method_1_joined_percentKinases_df.groupby("Kinase")["Rank"].sum().to_frame()
    method_2_hit_count_df = method_2_joined_percentKinases_df.groupby("Kinase")["Rank"].sum().to_frame()
    
    # Rename columns and join results
    method_1_hit_count_df.columns = [method_1_name]
    method_2_hit_count_df.columns = [method_2_name]
    
    hit_count_df = pd.concat([method_1_hit_count_df, method_2_hit_count_df], axis=1)
    hit_count_df = hit_count_df.fillna(0)
    
    # Save results
    #hit_count_df.to_csv(
    #    f"{out_prefix}cptac_{method_1_name}_{method_2_name}_hit_counts.tsv",
    #    header=True,
    #    index=True,
    #    sep="\t"
    #)

    # radial plots
    plt.close()
    labels = list(hit_count_df.index)
    values_1 = list(hit_count_df[method_1_name])
    values_2 = list(hit_count_df[method_2_name])
    num_vars = len(values_1)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    values_1 += values_1[:1]
    values_2 += values_2[:1]
    angles += angles[:1]
    max_value = max(max(values_1), max(values_2))
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
    ax.plot(angles, values_1, marker='o', label=method_1_name)
    ax.plot(angles, values_2, marker='x', linestyle='--', label=method_2_name)
    ax.set_xticklabels([])
    # Desired distance from circle edge to closest point of text
    desired_margin = max_value * 0.3
    for angle, label in zip(angles, labels):
        rotation = np.degrees(angle)
        # Estimate text dimensions (adjust multipliers based on your font)
        char_width = 0.02 * max_value  # approximate character width in data units
        char_height = 0.02 * max_value  # approximate character height in data units
        text_width = len(label) * char_width
        text_height = char_height
        # Calculate the radial extent of the rotated text box
        # When text is rotated by angle θ, the radial extent is:
        # radial_extent = |width * cos(θ)| + |height * sin(θ)|
        angle_rad = np.radians(rotation)
        radial_extent = abs(text_width * np.cos(angle_rad)) + abs(text_height * np.sin(angle_rad))
        # Position the text center so that the closest edge is at desired distance
        radial_position = max_value + desired_margin + radial_extent / 2
        
        if angle > np.pi/2 and angle < 3*np.pi/2:
            rotation = rotation + 180
            ha = 'center'
        else:
            ha = 'center'
        ax.text(angle, radial_position, label, 
                ha=ha, va='center', rotation=rotation)
    fig.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
    ax.legend(loc='upper right', frameon=False, bbox_to_anchor=(1.3, 1.3))
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_topNpercentKinases_radial.pdf")
    
    return hit_count_df
