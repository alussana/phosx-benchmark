#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from itertools import product
from random import sample


plt.rcParams['axes.titlesize'] = 8        
plt.rcParams['axes.labelsize'] = 6        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6 


def quantile_normalize(df: pd.DataFrame):
    """
    Perform quantile normalization on a pandas DataFrame.

    Args:
    df (pandas.DataFrame): DataFrame to be quantile normalized.

    Returns:
    pandas.DataFrame: Quantile normalized DataFrame.
    """
    # Rank the values in each column
    ranked = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    # Sort the ranks and use them as indices to get the quantiles
    quantiles = df.rank(method="min").stack().astype(int).map(ranked).unstack()
    # Sort the indices of the original DataFrame
    quantiles = quantiles.sort_index()
    return quantiles


def scale_01(df):
    # Subtract the minimum in each column
    scaled_df = df.apply(lambda x: x - x.min())
    # Divide by the maximum in each column
    scaled_df = scaled_df.apply(lambda x: x / x.max())
    return scaled_df


def compute_pr(labels, scores, n_points=200):
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    interp_recall = np.linspace(0, 1, n_points)
    interp_precision = np.interp(x=interp_recall, xp=precision, fp=recall)
    roc_auc = auc(recall, precision)
    return (interp_recall, interp_precision, roc_auc)


def compute_roc(labels, scores, n_points=200):
    fpr, tpr, thresholds = roc_curve(labels, scores)
    interp_fpr = np.linspace(0, 1, n_points)
    interp_tpr = np.interp(x=interp_fpr, xp=fpr, fp=tpr)
    roc_auc = auc(fpr, tpr)
    return (interp_fpr, interp_tpr, roc_auc)


def true_kinase_quantile(kinase_activity_df, regulation_df):
    def find_true_kinase_quantile(ex_series):
        kinase_str = ex_series["Kinase"]
        experiment_str = ex_series["Experiment"]
        regulation_int = ex_series["Regulation"]
        exp_activity_df = kinase_activity_df.loc[
            kinase_activity_df["Experiment"] == experiment_str
        ]
        if regulation_int == 1:
            exp_activity_df = exp_activity_df.sort_values(
                by="Kinase activity change", ascending=False
            )
        elif regulation_int == -1:
            exp_activity_df = exp_activity_df.sort_values(
                by="Kinase activity change", ascending=True
            )
        exp_activity_df["Rank"] = range(1, len(exp_activity_df) + 1)
        exp_activity_df["Rank"] = exp_activity_df["Rank"] / len(exp_activity_df)
        if kinase_str in exp_activity_df["Kinase"].values:
            quantile = float(
                exp_activity_df.loc[
                    exp_activity_df["Kinase"] == kinase_str, "Rank"
                ].values[0]
            )
        else:
            quantile = np.nan
        return quantile
    quantile_series = regulation_df.apply(find_true_kinase_quantile, axis=1)
    return quantile_series


def top10percentScoreFraction(
    kinase_activity_df, regulation_df, metadata_experiments_set
):
    top10percentScore = 0
    top10percentScore_total = 0
    for experiment in metadata_experiments_set:
        upreg_kinases = list(
            regulation_df["Kinase"].loc[regulation_df["Experiment"] == experiment,]
        )
        if len(upreg_kinases) == 0:
            next
        else:
            for kinase in upreg_kinases:
                activities = kinase_activity_df.loc[
                    kinase_activity_df["Experiment"] == experiment
                ]
                activity_k = activities["Kinase activity change"].loc[
                    activities["Kinase"] == kinase
                ]
                if len(activity_k) == 0:
                    next
                else:
                    top10percentScore_total = top10percentScore_total + 1
                    if activity_k.values[0] > 0.9:
                        top10percentScore = top10percentScore + 1
    return (top10percentScore, top10percentScore_total)


def bottom10percentScoreFraction(
    kinase_activity_df, regulation_df, metadata_experiments_set
):
    bottom10percentScore = 0
    bottom10percentScore_total = 0
    for experiment in metadata_experiments_set:
        upreg_kinases = list(
            regulation_df["Kinase"].loc[regulation_df["Experiment"] == experiment,]
        )
        if len(upreg_kinases) == 0:
            next
        else:
            for kinase in upreg_kinases:
                activities = kinase_activity_df.loc[
                    kinase_activity_df["Experiment"] == experiment
                ]
                activity_k = activities["Kinase activity change"].loc[
                    activities["Kinase"] == kinase
                ]
                if len(activity_k) == 0:
                    next
                else:
                    bottom10percentScore_total = bottom10percentScore_total + 1
                    if activity_k.values[0] < 0.1:
                        bottom10percentScore = bottom10percentScore + 1
    return (bottom10percentScore, bottom10percentScore_total)


def top10percentKinasesFraction(
    kinase_activity_df, regulation_df, metadata_experiments_set
):
    top10percentScore = 0
    top10percentScore_total = 0
    for experiment in metadata_experiments_set:
        upreg_kinases = list(
            regulation_df["Kinase"].loc[regulation_df["Experiment"] == experiment,]
        )
        if len(upreg_kinases) == 0:
            next
        else:
            activities = kinase_activity_df.loc[
                kinase_activity_df["Experiment"] == experiment
            ]
            # uniq_vals = activities['Kinase activity change'].sort_values(ascending=False).unique()
            uniq_vals = (
                activities["Kinase activity change"]
                .sort_values(ascending=True)
                .unique()
            )
            ranks = np.array([i for i in range(len(uniq_vals))])
            vals2ranks = dict(zip(uniq_vals, ranks))
            activities["Rank"] = [
                vals2ranks[i] for i in activities["Kinase activity change"]
            ]
            n_ranks = activities["Rank"].max()
            for kinase in upreg_kinases:
                rank_k = activities["Rank"].loc[activities["Kinase"] == kinase]
                if len(rank_k) == 0:
                    next
                else:
                    top10percentScore_total = top10percentScore_total + 1
                    if rank_k.values[0] / n_ranks > 0.90:
                        top10percentScore = top10percentScore + 1
    return (top10percentScore, top10percentScore_total)


def bottom10percentKinasesFraction(
    kinase_activity_df, regulation_df, metadata_experiments_set
):
    top10percentScore = 0
    top10percentScore_total = 0
    for experiment in metadata_experiments_set:
        downreg_kinases = list(
            regulation_df["Kinase"].loc[regulation_df["Experiment"] == experiment,]
        )
        if len(downreg_kinases) == 0:
            next
        else:
            activities = kinase_activity_df.loc[
                kinase_activity_df["Experiment"] == experiment
            ]
            # uniq_vals = activities['Kinase activity change'].sort_values(ascending=False).unique()
            uniq_vals = (
                activities["Kinase activity change"]
                .sort_values(ascending=True)
                .unique()
            )
            ranks = np.array([i for i in range(len(uniq_vals))])
            vals2ranks = dict(zip(uniq_vals, ranks))
            activities["Rank"] = [
                vals2ranks[i] for i in activities["Kinase activity change"]
            ]
            n_ranks = activities["Rank"].max()
            for kinase in downreg_kinases:
                rank_k = activities["Rank"].loc[activities["Kinase"] == kinase]
                if len(rank_k) == 0:
                    next
                else:
                    top10percentScore_total = top10percentScore_total + 1
                    if rank_k.values[0] / n_ranks < 0.1:
                        top10percentScore = top10percentScore + 1
    return (top10percentScore, top10percentScore_total)


"""def make_kinase_activity_df(
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
    return kinase_activity_df"""


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
    method_1_name:str,
    kinase_activity_method_1_df:pd.DataFrame,
    method_2_name:str,
    kinase_activity_method_2_df:pd.DataFrame,
    metadata_df:pd.DataFrame,
    out_prefix:str,
):
    # Take only instances for which a kinase activity could be computed by both methods
    intersection_index = list(
        set(kinase_activity_method_1_df.index)
        .intersection(set(kinase_activity_method_2_df.index))
    )
    kinase_activity_method_1_df = kinase_activity_method_1_df.loc[intersection_index,]
    kinase_activity_method_2_df = kinase_activity_method_2_df.loc[intersection_index,]
    # make sets of indexes for positive and negative examples in upregulation and downregulation, separately
    upregulation_df = metadata_df.loc[metadata_df["Regulation"] == 1]
    downregulation_df = metadata_df.loc[metadata_df["Regulation"] == -1]
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
    possible_indexes_list = list(
        product(metadata_kinases_set, metadata_experiments_set)
    )
    possible_indexes_set = set([f"{i[1]}__{i[0]}" for i in possible_indexes_list])
    upregulation_negative_indexes_list = list(
        possible_indexes_set.difference(set(upregulation_df.index))
    )
    downregulation_negative_indexes_list = list(
        possible_indexes_set.difference(set(downregulation_df.index))
    )
    # true kinase quantile
    method_1_upreg_true_kinase_quantile_series = true_kinase_quantile(
        kinase_activity_method_1_df, upregulation_df
    )
    method_1_downreg_true_kinase_quantile_series = true_kinase_quantile(
        kinase_activity_method_1_df, downregulation_df
    )
    method_2_upreg_true_kinase_quantile_series = true_kinase_quantile(
        kinase_activity_method_2_df, upregulation_df
    )
    method_2_downreg_true_kinase_quantile_series = true_kinase_quantile(
        kinase_activity_method_2_df, downregulation_df
    )
    method_1_joined_true_kinase_quantile_series = pd.concat(
        [
            method_1_upreg_true_kinase_quantile_series,
            method_1_downreg_true_kinase_quantile_series,
        ]
    )
    method_2_joined_true_kinase_quantile_series = pd.concat(
        [
            method_2_upreg_true_kinase_quantile_series,
            method_2_downreg_true_kinase_quantile_series,
        ]
    )
    joined_true_kinase_quantile_df = (
        pd.DataFrame(
            {
                method_1_name: method_1_joined_true_kinase_quantile_series,
                method_2_name: method_2_joined_true_kinase_quantile_series,
            }
        )
        .dropna()
        .melt()
    )
    joined_true_kinase_quantile_df.columns = ["Method", "Normalized rank"]
    upreg_true_kinase_quantile_df = (
        pd.DataFrame(
            {
                method_1_name: method_1_upreg_true_kinase_quantile_series,
                method_2_name: method_2_upreg_true_kinase_quantile_series,
            }
        )
        .dropna()
        .melt()
    )
    upreg_true_kinase_quantile_df.columns = ["Method", "Normalized rank"]
    downreg_true_kinase_quantile_df = (
        pd.DataFrame(
            {
                method_1_name: method_1_downreg_true_kinase_quantile_series,
                method_2_name: method_2_downreg_true_kinase_quantile_series,
            }
        )
        .dropna()
        .melt()
    )
    downreg_true_kinase_quantile_df.columns = ["Method", "Normalized rank"]
    n_upreg = len(upreg_true_kinase_quantile_df) // 2
    n_downreg = len(downreg_true_kinase_quantile_df) // 2
    n_joined = len(joined_true_kinase_quantile_df) // 2
    plt.clf()
    plt.figure(figsize=[1.75, 2.5])
    ax = sns.violinplot(
        data=joined_true_kinase_quantile_df,
        x="Method",
        y="Normalized rank",
        cut=0,
        #hue="Method",
    )
    ax.set_title(
        f"Shared positive examples\n(n={n_joined})"
    )
    #ax.set_yticks(np.arange(0, 1.1, 0.1))
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_joined_true_kinase_quantile.pdf")
    plt.clf()
    plt.figure(figsize=[1.75, 2.5])
    ax = sns.violinplot(
        data=upreg_true_kinase_quantile_df,
        x="Method",
        y="Normalized rank",
        cut=0,
        #hue="Method",
    )
    ax.set_title(
        f"Shared positive examples\n(activation, n={n_upreg})"
    )
    #ax.set_yticks(np.arange(0, 1.1, 0.1))
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_upreg_true_kinase_quantile.pdf")
    plt.clf()
    plt.figure(figsize=[1.75, 2.5])
    ax = sns.violinplot(
        data=downreg_true_kinase_quantile_df,
        x="Method",
        y="Normalized rank",
        cut=0,
        #hue="Method",
    )
    ax.set_title(
        f"Shared positive examples\n(inhibition, n={n_downreg})"
    )
    #ax.set_yticks(np.arange(0, 1.1, 0.1))
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_downreg_true_kinase_quantile.pdf")
    # top 10% score true positives performance
    method_1_top10percentScore, method_1_top10percentScore_total = top10percentScoreFraction(
        kinase_activity_method_1_df, upregulation_df, metadata_experiments_set
    )
    method_2_top10percentScore, method_2_top10percentScore_total = top10percentScoreFraction(
        kinase_activity_method_2_df, upregulation_df, metadata_experiments_set
    )
    method_1_bottom10percentScore, method_1_bottom10percentScore_total = (
        bottom10percentScoreFraction(
            kinase_activity_method_1_df, downregulation_df, metadata_experiments_set
        )
    )
    method_2_bottom10percentScore, method_2_bottom10percentScore_total = (
        bottom10percentScoreFraction(
            kinase_activity_method_2_df, downregulation_df, metadata_experiments_set
        )
    )
    method_1_score_total = method_1_top10percentScore_total + method_1_bottom10percentScore_total
    method_2_score_total = method_2_top10percentScore_total + method_2_bottom10percentScore_total
    tp_percentage_df = pd.DataFrame(
        {
            "Upreg.": [
                method_1_top10percentScore / method_1_score_total,
                method_2_top10percentScore / method_2_score_total,
            ],
            "Downreg.": [
                method_1_bottom10percentScore / method_1_score_total,
                method_2_bottom10percentScore / method_2_score_total,
            ],
            "Method": [method_1_name, method_2_name],
        }
    )
    plt.clf()
    plt.figure(figsize=[2, 3])
    ax = tp_percentage_df.set_index("Method").plot(
        kind="bar", stacked=True, color=["red", "blue"], alpha=0.5, figsize=(2.2, 2.6)
    )
    plt.xticks(rotation=0)
    plt.legend(frameon=False)
    plt.ylim([0.0, 1])
    ax.set_ylabel("Fraction of regulated kinases")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_topNpercentScore_stacked_barplot.pdf")
    # top 10% kinases true positives performance
    method_1_top10percentKinases, method_1_top10percentKinases_total = (
        top10percentKinasesFraction(
            kinase_activity_method_1_df, upregulation_df, metadata_experiments_set
        )
    )
    method_2_top10percentKinases, method_2_top10percentKinases_total = (
        top10percentKinasesFraction(
            kinase_activity_method_2_df, upregulation_df, metadata_experiments_set
        )
    )
    method_1_bottom10percentKinases, method_1_bottom10percentKinases_total = (
        bottom10percentKinasesFraction(
            kinase_activity_method_1_df, downregulation_df, metadata_experiments_set
        )
    )
    method_2_bottom10percentKinases, method_2_bottom10percentKinases_total = (
        bottom10percentKinasesFraction(
            kinase_activity_method_2_df, downregulation_df, metadata_experiments_set
        )
    )
    method_1_kinases_total = (
        method_1_top10percentKinases_total + method_1_bottom10percentKinases_total
    )
    method_2_kinases_total = (
        method_2_top10percentKinases_total + method_2_bottom10percentKinases_total
    )
    tp_percentage_df = pd.DataFrame(
        {
            "Upreg.": [
                method_1_top10percentKinases / method_1_kinases_total,
                method_2_top10percentKinases / method_2_kinases_total,
            ],
            "Downreg.": [
                method_1_bottom10percentKinases / method_1_kinases_total,
                method_2_bottom10percentKinases / method_2_kinases_total,
            ],
            "Method": [method_1_name, method_2_name],
        }
    )
    plt.clf()
    plt.figure(figsize=[2, 3])
    ax = tp_percentage_df.set_index("Method").plot(
        kind="bar", stacked=True, color=["red", "blue"], alpha=0.5, figsize=(2.2, 2.6)
    )
    plt.xticks(rotation=0)
    plt.legend(frameon=False)
    plt.ylim([0.0, 1])
    ax.set_ylabel("Fraction of regulated kinases")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_topNpercentKinases_stacked_barplot.pdf")
    # AUROC and AUPR method 1 evaluation 
    violinplots_upreg_df = pd.DataFrame()
    violinplots_downreg_df = pd.DataFrame()
    violinplots_joined_df = pd.DataFrame()
    upreg_auc_list = []
    upreg_apr_list = []
    downreg_auc_list = []
    downreg_apr_list = []
    upregulation_method_1_positive_indexes_list = list(
        set(kinase_activity_method_1_df.index).intersection(set(upregulation_df.index))
    )
    upregulation_method_1_negative_indexes_list = [
        x
        for x in upregulation_negative_indexes_list
        if x in kinase_activity_method_1_df.index
    ]
    downregulation_method_1_positive_indexes_list = list(
        set(kinase_activity_method_1_df.index).intersection(set(downregulation_df.index))
    )
    downregulation_method_1_negative_indexes_list = [
        x
        for x in downregulation_negative_indexes_list
        if x in kinase_activity_method_1_df.index
    ]
    for i in range(256):
        # randomly sample neg ex in equal number as the pos ex, for up- and down- regulation separately
        upreg_neg_idx_list = sample(
            upregulation_method_1_negative_indexes_list,
            len(upregulation_method_1_positive_indexes_list),
        )
        downreg_neg_idx_list = sample(
            downregulation_method_1_negative_indexes_list,
            len(downregulation_method_1_positive_indexes_list),
        )
        # build the evaluation datasets for up- and down- regulation separately
        upreg_df = kinase_activity_method_1_df.loc[
            upregulation_method_1_positive_indexes_list + upreg_neg_idx_list
        ]
        downreg_df = kinase_activity_method_1_df.loc[
            downregulation_method_1_positive_indexes_list + downreg_neg_idx_list
        ]
        upreg_df = upreg_df.merge(upregulation_df, how="left")
        upreg_df["Regulation"].loc[upreg_df["Regulation"].isna()] = 0
        upreg_df.index = upreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
        )
        downreg_df = downreg_df.merge(downregulation_df, how="left")
        downreg_df["Regulation"].loc[downreg_df["Regulation"].isna()] = 0
        downreg_df["Regulation"].loc[downreg_df["Regulation"] == -1] = 1
        downreg_df.index = downreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
        )
        # upregulation
        fpr, tpr, thresholds = roc_curve(
            upreg_df["Regulation"], upreg_df["Kinase activity change"]
        )
        upreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(
            upreg_df["Regulation"], upreg_df["Kinase activity change"]
        )
        upreg_apr_list.append(auc(recall, precision))
        # downregulation
        fpr, tpr, thresholds = roc_curve(
            downreg_df["Regulation"], -downreg_df["Kinase activity change"]
        )
        downreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(
            downreg_df["Regulation"], -downreg_df["Kinase activity change"]
        )
        downreg_apr_list.append(auc(recall, precision))
    n_upreg_method_1 = len(upreg_df.loc[upreg_df["Regulation"] == 1])
    n_downreg_method_1 = len(downreg_df.loc[downreg_df["Regulation"] == 1])
    data = pd.DataFrame({"AUROC": upreg_auc_list, "AUPR": upreg_apr_list})
    data = data.melt()
    data.columns = ["Metric", "Value"]
    data["Method"] = [method_1_name for i in range(len(data))]
    violinplots_upreg_df = pd.concat([violinplots_upreg_df, data])
    data = pd.DataFrame({"AUROC": downreg_auc_list, "AUPR": downreg_apr_list})
    data = data.melt()
    data.columns = ["Metric", "Value"]
    data["Method"] = [method_1_name for i in range(len(data))]
    violinplots_downreg_df = pd.concat([violinplots_downreg_df, data])
    data = pd.DataFrame(
        {
            "AUROC": upreg_auc_list + downreg_auc_list,
            "AUPR": upreg_apr_list + downreg_apr_list,
        }
    )
    data = data.melt()
    data.columns = ["Metric", "Value"]
    data["Method"] = [method_1_name for i in range(len(data))]
    violinplots_joined_df = pd.concat([violinplots_joined_df, data])
    # AUROC and AUPR method 2 evaluation 
    upreg_auc_list = []
    upreg_apr_list = []
    downreg_auc_list = []
    downreg_apr_list = []
    upregulation_method_2_positive_indexes_list = list(
        set(kinase_activity_method_2_df.index).intersection(set(upregulation_df.index))
    )
    upregulation_method_2_negative_indexes_list = [
        x
        for x in upregulation_negative_indexes_list
        if x in kinase_activity_method_2_df.index
    ]
    downregulation_method_2_positive_indexes_list = list(
        set(kinase_activity_method_2_df.index).intersection(set(downregulation_df.index))
    )
    downregulation_method_2_negative_indexes_list = [
        x
        for x in downregulation_negative_indexes_list
        if x in kinase_activity_method_2_df.index
    ]
    for i in range(256):
        # randomly sample neg ex in equal number as the pos ex, for up- and down- regulation separately
        upreg_neg_idx_list = sample(
            upregulation_method_2_negative_indexes_list,
            len(upregulation_method_2_positive_indexes_list),
        )
        downreg_neg_idx_list = sample(
            downregulation_method_2_negative_indexes_list,
            len(downregulation_method_2_positive_indexes_list),
        )
        # build the evaluation datasets for up- and down- regulation separately
        upreg_df = kinase_activity_method_2_df.loc[
            upregulation_method_2_positive_indexes_list + upreg_neg_idx_list
        ]
        downreg_df = kinase_activity_method_2_df.loc[
            downregulation_method_2_positive_indexes_list + downreg_neg_idx_list
        ]
        upreg_df = upreg_df.merge(upregulation_df, how="left")
        upreg_df["Regulation"].loc[upreg_df["Regulation"].isna()] = 0
        upreg_df.index = upreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
        )
        downreg_df = downreg_df.merge(downregulation_df, how="left")
        downreg_df["Regulation"].loc[downreg_df["Regulation"].isna()] = 0
        downreg_df["Regulation"].loc[downreg_df["Regulation"] == -1] = 1
        downreg_df.index = downreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}', axis=1
        )
        # upregulation
        fpr, tpr, thresholds = roc_curve(
            upreg_df["Regulation"], upreg_df["Kinase activity change"]
        )
        upreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(
            upreg_df["Regulation"], upreg_df["Kinase activity change"]
        )
        upreg_apr_list.append(auc(recall, precision))
        # downregulation
        fpr, tpr, thresholds = roc_curve(
            downreg_df["Regulation"], -downreg_df["Kinase activity change"]
        )
        downreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(
            downreg_df["Regulation"], -downreg_df["Kinase activity change"]
        )
        downreg_apr_list.append(auc(recall, precision))
    n_upreg_method_2 = len(upreg_df.loc[upreg_df["Regulation"] == 1])
    n_downreg_method_2 = len(downreg_df.loc[downreg_df["Regulation"] == 1])
    data = pd.DataFrame({"AUROC": upreg_auc_list, "AUPR": upreg_apr_list})
    data = data.melt()
    data.columns = ["Metric", "Value"]
    data["Method"] = [method_2_name for i in range(len(data))]
    violinplots_upreg_df = pd.concat([violinplots_upreg_df, data])
    data = pd.DataFrame({"AUROC": downreg_auc_list, "AUPR": downreg_apr_list})
    data = data.melt()
    data.columns = ["Metric", "Value"]
    data["Method"] = [method_2_name for i in range(len(data))]
    violinplots_downreg_df = pd.concat([violinplots_downreg_df, data])
    data = pd.DataFrame(
        {
            "AUROC": upreg_auc_list + downreg_auc_list,
            "AUPR": upreg_apr_list + downreg_apr_list,
        }
    )
    data = data.melt()
    data.columns = ["Metric", "Value"]
    data["Method"] = [method_2_name for i in range(len(data))]
    violinplots_joined_df = pd.concat([violinplots_joined_df, data])
    plt.clf()
    plt.figure(figsize=[1.75, 2.5])
    ax = sns.violinplot(
        data=violinplots_upreg_df, x="Metric", y="Value", hue="Method", cut=0
    )
    plt.ylim([0.0, 1.0])
    plt.axhline(y=0.5, color="grey", linestyle="--", linewidth=1)
    ax.set_title(
        f"Shared positive examples\n(activation, n={n_upreg_method_1})"
    )
    plt.legend(loc="lower right", frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_upreg_auc_prc_violinplots.pdf")

    plt.clf()
    plt.figure(figsize=[1.75, 2.5])
    ax = sns.violinplot(
        data=violinplots_downreg_df, x="Metric", y="Value", hue="Method", cut=0
    )
    plt.ylim([0.0, 1.0])
    plt.axhline(y=0.5, color="grey", linestyle="--", linewidth=1)
    ax.set_title(
        f"Shared positive examples\n(inhibition, n={n_downreg_method_1})"
    )
    plt.legend(loc="lower right", frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_downreg_auc_prc_violinplots.pdf")
    plt.clf()
    plt.figure(figsize=[1.75, 2.5])
    ax = sns.violinplot(
        data=violinplots_joined_df, x="Metric", y="Value", hue="Method", cut=0
    )
    plt.ylim([0.0, 1.0])
    plt.axhline(y=0.5, color="grey", linestyle="--", linewidth=1)
    ax.set_title(
        f"Shared positive examples\n(n={n_upreg_method_1 + n_downreg_method_1})"
    )
    plt.legend(loc="lower right", frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_joined_auc_prc_violinplots.pdf")
    # Upregulated and Downregulated subset comparisons
    kinase_activity_method_1_df = kinase_activity_method_1_df.rename(
        columns={"Kinase activity change": f"{method_1_name} Activity Score"}
    )
    kinase_activity_method_2_df = kinase_activity_method_2_df.rename(
        columns={"Kinase activity change": f"{method_2_name} Activity Score"}
    )
    regulated_df = (
        metadata_df.merge(
            kinase_activity_method_1_df, on=["Experiment", "Kinase"], how="left"
        )
        .merge(kinase_activity_method_2_df, on=["Experiment", "Kinase"], how="left")
    )
    plt.clf()
    plt.figure(figsize=(3, 3))
    #plt.ylim([0, 30])
    sns.histplot(
        regulated_df[f"{method_1_name} Activity Score"].loc[regulated_df["Regulation"] == 1],
        color="red",
        alpha=0.5,
        binrange=[0, 1],
        binwidth=0.05,
        legend="Upregulated",
    )
    sns.histplot(
        regulated_df[f"{method_1_name} Activity Score"].loc[regulated_df["Regulation"] == -1],
        color="blue",
        alpha=0.5,
        binrange=[0, 1],
        binwidth=0.05,
        legend="Downregulated",
    )
    plt.legend(loc="upper left", frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_in_{method_1_name}_{method_2_name}_score_regulated_kinases.pdf")
    plt.clf()
    plt.figure(figsize=(3, 3))
    #plt.ylim([0, 30])
    sns.histplot(
        regulated_df[f"{method_2_name} Activity Score"].loc[regulated_df["Regulation"] == 1],
        color="red",
        alpha=0.5,
        binrange=[0, 1],
        binwidth=0.05,
        legend="Upregulated",
    )
    sns.histplot(
        regulated_df[f"{method_2_name} Activity Score"].loc[regulated_df["Regulation"] == -1],
        color="blue",
        alpha=0.5,
        binrange=[0, 1],
        binwidth=0.05,
        legend="Downregulated",
    )
    plt.legend(loc="upper left", frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_2_name}_in_{method_1_name}_{method_2_name}_score_regulated_kinases.pdf")
    regulated_df = regulated_df.dropna()
    n_upregulated = len(regulated_df.loc[regulated_df["Regulation"] == 1])
    n_downregulated = len(regulated_df.loc[regulated_df["Regulation"] == -1])
    method_1_fpr, method_1_tpr, method_1_roc_auc = compute_roc(
        regulated_df["Regulation"], regulated_df[f"{method_1_name} Activity Score"]
    )
    method_2_fpr, method_2_tpr, method_2_roc_auc = compute_roc(
        regulated_df["Regulation"], regulated_df[f"{method_2_name} Activity Score"]
    )
    class_imbalance = n_upregulated / (n_upregulated + n_downregulated)
    if class_imbalance < 0.5:
        class_imbalance = 1 / class_imbalance
    plt.clf()
    plt.figure(figsize=(3, 3))
    plt.plot(method_1_fpr, method_1_tpr, lw=2, label=f"{method_1_name} (AUC = {method_1_roc_auc:.2f})")
    plt.plot(method_2_fpr, method_2_tpr, lw=2, label=f"{method_2_name} (AUC = {method_2_roc_auc:.2f})")
    plt.plot([0, 1], [0, 1], color="black", lw=0.5, linestyle="-")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.legend(loc="lower right", frameon=False)
    plt.title(
        f"Shared positive examples\n(activation, n={n_upregulated}; inhibition, n={n_downregulated})"
    )
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_regulated_kinases_ROC.pdf")
    method_1_recall, method_1_precision, method_1_pr_auc = compute_pr(
        regulated_df["Regulation"], regulated_df[f"{method_1_name} Activity Score"]
    )
    method_2_recall, method_2_precision, method_2_pr_auc = compute_pr(
        regulated_df["Regulation"], regulated_df[f"{method_2_name} Activity Score"]
    )
    plt.clf()
    plt.figure(figsize=(3, 3))
    plt.plot(
        method_1_recall, method_1_precision, lw=2, label=f"{method_1_name} (AUC = {method_1_pr_auc:.2f})"
    )
    plt.plot(
        method_2_recall, method_2_precision, lw=2, label=f"{method_2_name} (AUC = {method_2_pr_auc:.2f})"
    )
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.axhline(
        y=class_imbalance,
        lw=0.5,
        color="black",
        linestyle="-",
    )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend(loc="lower left", frameon=False)
    plt.title(
        f"Shared positive examples\n(activation, n={n_upregulated}; inhibition, n={n_downregulated})"
    )
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}cptac_{method_1_name}_{method_2_name}_regulated_kinases_PR.pdf")