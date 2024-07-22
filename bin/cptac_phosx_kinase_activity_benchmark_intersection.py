#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from itertools import product
from random import sample

def quantile_normalize(df:pd.DataFrame):
    """
    Perform quantile normalization on a pandas DataFrame.

    Args:
    df (pandas.DataFrame): DataFrame to be quantile normalized.

    Returns:
    pandas.DataFrame: Quantile normalized DataFrame.
    """
    # Rank the values in each column
    ranked = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    
    # Sort the ranks and use them as indices to get the quantiles
    quantiles = df.rank(method='min').stack().astype(int).map(ranked).unstack()
    
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
    
    return(interp_recall, interp_precision, roc_auc)

def compute_roc(labels, scores, n_points=200):
    fpr, tpr, thresholds = roc_curve(labels, scores)
    
    interp_fpr = np.linspace(0, 1, n_points)
    interp_tpr = np.interp(x=interp_fpr, xp=fpr, fp=tpr)
    
    roc_auc = auc(fpr, tpr)
    
    return(interp_fpr, interp_tpr, roc_auc)


def top5percentScoreFraction(kinase_activity_df, regulation_df, metadata_experiments_set):
    top5percentScore = 0
    top5percentScore_total = 0
    for experiment in metadata_experiments_set:
        upreg_kinases = list(regulation_df['Kinase'].loc[regulation_df['Experiment']==experiment,])
        if len(upreg_kinases) == 0:
            next
        else:
            for kinase in upreg_kinases:
                activities = kinase_activity_df.loc[kinase_activity_df['Experiment']==experiment]
                activity_k = activities['Kinase activity change'].loc[activities['Kinase']==kinase]
                if len(activity_k) == 0:
                    next
                else:
                    top5percentScore_total = top5percentScore_total + 1
                    if activity_k.values[0] >= 0.95:
                        top5percentScore = top5percentScore + 1
    return(top5percentScore, top5percentScore_total)


def bottom5percentScoreFraction(kinase_activity_df, regulation_df, metadata_experiments_set):
    bottom5percentScore = 0
    bottom5percentScore_total = 0
    for experiment in metadata_experiments_set:
        upreg_kinases = list(regulation_df['Kinase'].loc[regulation_df['Experiment']==experiment,])
        if len(upreg_kinases) == 0:
            next
        else:
            for kinase in upreg_kinases:
                activities = kinase_activity_df.loc[kinase_activity_df['Experiment']==experiment]
                activity_k = activities['Kinase activity change'].loc[activities['Kinase']==kinase]
                if len(activity_k) == 0:
                    next
                else:
                    bottom5percentScore_total = bottom5percentScore_total + 1
                    if activity_k.values[0] <= 0.05:
                        bottom5percentScore = bottom5percentScore + 1
    return(bottom5percentScore, bottom5percentScore_total)


def top5percentKinasesFraction(kinase_activity_df, regulation_df, metadata_experiments_set):
    top5percentScore = 0
    top5percentScore_total = 0
    for experiment in metadata_experiments_set:
        upreg_kinases = list(regulation_df['Kinase'].loc[regulation_df['Experiment']==experiment,])
        if len(upreg_kinases) == 0:
            next
        else:
            activities = kinase_activity_df.loc[kinase_activity_df['Experiment']==experiment]
            #uniq_vals = activities['Kinase activity change'].sort_values(ascending=False).unique()
            uniq_vals = activities['Kinase activity change'].sort_values(ascending=True).unique()
            ranks = np.array([i for i in range(len(uniq_vals))])
            vals2ranks = dict(zip(uniq_vals, ranks))
            activities['Rank'] = [vals2ranks[i] for i in activities['Kinase activity change']]
            n_ranks = activities['Rank'].max()
            for kinase in upreg_kinases:
                rank_k = activities['Rank'].loc[activities['Kinase']==kinase]
                if len(rank_k) == 0:
                    next
                else:
                    top5percentScore_total = top5percentScore_total + 1
                    if rank_k.values[0]/n_ranks >= 0.95:
                        top5percentScore = top5percentScore + 1
    return(top5percentScore, top5percentScore_total)


def bottom5percentKinasesFraction(kinase_activity_df, regulation_df, metadata_experiments_set):
    top5percentScore = 0
    top5percentScore_total = 0
    for experiment in metadata_experiments_set:
        downreg_kinases = list(regulation_df['Kinase'].loc[regulation_df['Experiment']==experiment,])
        if len(downreg_kinases) == 0:
            next
        else:
            activities = kinase_activity_df.loc[kinase_activity_df['Experiment']==experiment]
            #uniq_vals = activities['Kinase activity change'].sort_values(ascending=False).unique()
            uniq_vals = activities['Kinase activity change'].sort_values(ascending=True).unique()
            ranks = np.array([i for i in range(len(uniq_vals))])
            vals2ranks = dict(zip(uniq_vals, ranks))
            activities['Rank'] = [vals2ranks[i] for i in activities['Kinase activity change']]
            n_ranks = activities['Rank'].max()
            for kinase in downreg_kinases:
                rank_k = activities['Rank'].loc[activities['Kinase']==kinase]
                if len(rank_k) == 0:
                    next
                else:
                    top5percentScore_total = top5percentScore_total + 1
                    if rank_k.values[0]/n_ranks <= 0.05:
                        top5percentScore = top5percentScore + 1
    return(top5percentScore, top5percentScore_total)


def main():
    input_list_phosx_txt = sys.argv[1]
    input_list_gsea_txt = sys.argv[2]
    input_list_kinex_txt = sys.argv[3]
    metadata_tsv = sys.argv[4]
    kinase_activity_metric_str = sys.argv[5]
    out_prefix = sys.argv[6]
    """
    input_list_phosx_txt = 'input_files_phosx.txt'
    input_list_gsea_txt = 'input_files_gsea.txt'
    input_list_kinex_txt = 'input_files_kinex.txt'
    metadata_tsv = 'input/metadata.tsv'
    kinase_activity_metric_str = 'Activity Score'
    out_prefix = 'phosx/kinase_activity/'
    """
    
    gsea_kinase_activity_metric_str = 'NES'
   
    # PhosX kinase activity 
    # build dictionary (experiment id --> kinase activity table path)
    data_path_dict = {}
    with open(input_list_phosx_txt, 'r') as input_list_fh:
        for line in input_list_fh:
            data_path_str = line.strip()
            data_id = os.path.basename(data_path_str)[:-4]
            data_path_dict[data_id] = data_path_str
    # build dataframe of kinase activities
    kinase_activity_phosx_df = pd.DataFrame()
    for data_id in data_path_dict.keys():
        df = pd.read_csv(
            data_path_dict[data_id],
            sep='\t',
            index_col=0
        )
        series = df[kinase_activity_metric_str]
        series.name = data_id
        kinase_activity_phosx_df = kinase_activity_phosx_df.join(series, how='outer')
    # normalise and scale
    kinase_activity_phosx_df = quantile_normalize(kinase_activity_phosx_df)
    kinase_activity_phosx_df = scale_01(kinase_activity_phosx_df)
    # melt
    kinase_activity_phosx_df = pd.melt(
        kinase_activity_phosx_df.assign(index=kinase_activity_phosx_df.index),
        id_vars=['index']
    )
    kinase_activity_phosx_df.columns = ['Kinase','Experiment','Kinase activity change']
    # remove NAs due to missing phosphosite sequence
    kinase_activity_phosx_df = kinase_activity_phosx_df.dropna()
    # make index
    kinase_activity_phosx_df.index = kinase_activity_phosx_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
        axis=1
    )
    kinase_activity_phosx_df.to_csv(f'{out_prefix}cptac_kinase_activity_phosx.tsv', header=True, index=True, sep='\t')
    ##########
    
    # GSEApy kinase activity
    # build dictionary (experiment id --> kinase activity table path)
    data_path_dict = {}
    with open(input_list_gsea_txt, 'r') as input_list_fh:
        for line in input_list_fh:
            data_path_str = line.strip()
            data_id = os.path.basename(data_path_str)[:-4]
            data_path_dict[data_id] = data_path_str
    # build dataframe of kinase activities
    kinase_activity_gsea_df = pd.DataFrame()
    for data_id in data_path_dict.keys():
        df = pd.read_csv(
            data_path_dict[data_id],
            sep=',',
            index_col=0
        )
        series = df[gsea_kinase_activity_metric_str]
        series.name = data_id
        kinase_activity_gsea_df = kinase_activity_gsea_df.join(series, how='outer')
    # normalise and scale
    kinase_activity_gsea_df = quantile_normalize(kinase_activity_gsea_df)
    kinase_activity_gsea_df = scale_01(kinase_activity_gsea_df)
    # melt
    kinase_activity_gsea_df = pd.melt(
        kinase_activity_gsea_df.assign(index=kinase_activity_gsea_df.index),
        id_vars=['index']
    )
    kinase_activity_gsea_df.columns = ['Kinase','Experiment','Kinase activity change']
    # remove NAs due to no annotated substrate found for given kinases in the data
    kinase_activity_gsea_df = kinase_activity_gsea_df.dropna()
    # make index
    kinase_activity_gsea_df.index = kinase_activity_gsea_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
        axis=1
    )
    kinase_activity_gsea_df.to_csv(f'{out_prefix}cptac_kinase_activity_gsea.tsv', header=True, index=True, sep='\t')
    ##########
    
    # Kinex kinase activity 
    # build dictionary (experiment id --> kinase activity table path)
    data_path_dict = {}
    with open(input_list_kinex_txt, 'r') as input_list_fh:
        for line in input_list_fh:
            data_path_str = line.strip()
            data_id = os.path.basename(data_path_str)[:-4]
            data_path_dict[data_id] = data_path_str
    # build dataframe of kinase activities
    kinase_activity_kinex_df = pd.DataFrame()
    for data_id in data_path_dict.keys():
        df = pd.read_csv(
            data_path_dict[data_id],
            sep='\t',
            index_col=0
        )
        series = df[kinase_activity_metric_str]
        series.name = data_id
        kinase_activity_kinex_df = kinase_activity_kinex_df.join(series, how='outer')
    # normalise and scale
    kinase_activity_kinex_df = quantile_normalize(kinase_activity_kinex_df)
    kinase_activity_kinex_df = scale_01(kinase_activity_kinex_df)
    # melt
    kinase_activity_kinex_df = pd.melt(
        kinase_activity_kinex_df.assign(index=kinase_activity_kinex_df.index),
        id_vars=['index']
    )
    kinase_activity_kinex_df.columns = ['Kinase','Experiment','Kinase activity change']
    # remove NAs due to missing phosphosite sequence
    kinase_activity_kinex_df = kinase_activity_kinex_df.dropna()
    # make index
    kinase_activity_kinex_df.index = kinase_activity_kinex_df.apply(
        lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
        axis=1
    )
    kinase_activity_kinex_df.to_csv(f'{out_prefix}cptac_kinase_activity_kinex.tsv', header=True, index=True, sep='\t')
    ##########
  
    # Take only instances for which a kinase activity could be computed by all methods
    intersection_index = list(set(kinase_activity_phosx_df.index).intersection(set(kinase_activity_kinex_df.index)).intersection(set(kinase_activity_gsea_df.index)))
    kinase_activity_phosx_df = kinase_activity_phosx_df.loc[intersection_index,]
    kinase_activity_kinex_df = kinase_activity_kinex_df.loc[intersection_index,]
    kinase_activity_gsea_df = kinase_activity_gsea_df.loc[intersection_index,]
   
    # metadata - ground truth kinase regulation
    metadata_df = pd.read_csv(metadata_tsv, sep='\t', index_col=None, header=None)
    metadata_df.columns = ['Experiment','Kinase','Regulation']
    metadata_df.index = metadata_df.apply(
        lambda x:  f'{x["Experiment"]}__{x["Kinase"]}',
        axis=1
    )
    ##########
   
    # make sets of indexes for positive and negative examples in upregulation and downregulation, separately 
    upregulation_df = metadata_df.loc[metadata_df['Regulation']==1]
    downregulation_df = metadata_df.loc[metadata_df['Regulation']==-1]
    
    upregulation_df.to_csv(f'{out_prefix}cptac_upregulated.tsv', header=True, index=True, sep='\t')
    downregulation_df.to_csv(f'{out_prefix}cptac_downregulated.tsv', header=True, index=True, sep='\t')
    
    metadata_kinases_set = set(metadata_df['Kinase'].unique())
    metadata_experiments_set = set(metadata_df['Experiment'].unique())
    possible_indexes_list = list(product(metadata_kinases_set, metadata_experiments_set))
    possible_indexes_set = set([f'{i[1]}__{i[0]}' for i in possible_indexes_list])
    
    upregulation_negative_indexes_list = list(possible_indexes_set.difference(set(upregulation_df.index)))
    downregulation_negative_indexes_list = list(possible_indexes_set.difference(set(downregulation_df.index)))
    ##########
    
    
    # top 5% score true positives performance
    phosx_top5percentScore, phosx_top5percentScore_total = top5percentScoreFraction(kinase_activity_phosx_df, upregulation_df, metadata_experiments_set)
    kinex_top5percentScore, kinex_top5percentScore_total = top5percentScoreFraction(kinase_activity_kinex_df, upregulation_df, metadata_experiments_set)
    gsea_top5percentScore, gsea_top5percentScore_total = top5percentScoreFraction(kinase_activity_gsea_df, upregulation_df, metadata_experiments_set)
    phosx_bottom5percentScore, phosx_bottom5percentScore_total = bottom5percentScoreFraction(kinase_activity_phosx_df, downregulation_df, metadata_experiments_set)
    kinex_bottom5percentScore, kinex_bottom5percentScore_total = bottom5percentScoreFraction(kinase_activity_kinex_df, downregulation_df, metadata_experiments_set)
    gsea_bottom5percentScore, gsea_bottom5percentScore_total = bottom5percentScoreFraction(kinase_activity_gsea_df, downregulation_df, metadata_experiments_set)
    ##########
    
    """
    phosx_top5percentile, phosx_top5percentile_total, phosx_top5percentile/phosx_top5percentile_total
    gsea_top5percentile, gsea_top5percentile_total, gsea_top5percentile/gsea_top5percentile_total
    kinex_top5percentile, kinex_top5percentile_total, kinex_top5percentile/kinex_top5percentile_total
    phosx_bottom5percentile, phosx_bottom5percentile_total, phosx_bottom5percentile/phosx_bottom5percentile_total
    gsea_bottom5percentile, gsea_bottom5percentile_total, gsea_bottom5percentile/gsea_bottom5percentile_total
    kinex_bottom5percentile, kinex_bottom5percentile_total, kinex_bottom5percentile/kinex_bottom5percentile_total
    """
    
    phosx_score_total = phosx_top5percentScore_total + phosx_bottom5percentScore_total
    kinex_score_total = kinex_top5percentScore_total + kinex_bottom5percentScore_total
    gsea_score_total = gsea_top5percentScore_total + gsea_bottom5percentScore_total
    
    tp_percentage_df = pd.DataFrame({
        'Upregulation':
            [
                phosx_top5percentScore/phosx_score_total,
                kinex_top5percentScore/kinex_score_total,
                gsea_top5percentScore/gsea_score_total,
            ],
        'Downregulation':
            [
                phosx_bottom5percentScore/phosx_score_total,
                kinex_bottom5percentScore/kinex_score_total,
                gsea_bottom5percentScore/gsea_score_total
            ],
        'Method': 
            [
                'PhosX',
                'Kinex',
                'GSEApy'
            ],
    })
    
    plt.clf()
    ax = tp_percentage_df.set_index('Method').plot(kind='bar', stacked=True, color=['red','blue'], alpha=0.5, figsize=(2.5, 3))
    plt.xticks(rotation=0)
    plt.legend(frameon=False)
    plt.ylim([0.0, 1])
    ax.set_ylabel("Fraction of regulated kinases")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}topNpercentScore_stacked_barplot.pdf')
    
    # top 5% kinases true positives performance
    phosx_top5percentKinases, phosx_top5percentKinases_total = top5percentKinasesFraction(kinase_activity_phosx_df, upregulation_df, metadata_experiments_set)
    kinex_top5percentKinases, kinex_top5percentKinases_total = top5percentKinasesFraction(kinase_activity_kinex_df, upregulation_df, metadata_experiments_set)
    gsea_top5percentKinases, gsea_top5percentKinases_total = top5percentKinasesFraction(kinase_activity_gsea_df, upregulation_df, metadata_experiments_set)
    phosx_bottom5percentKinases, phosx_bottom5percentKinases_total = bottom5percentKinasesFraction(kinase_activity_phosx_df, downregulation_df, metadata_experiments_set)
    kinex_bottom5percentKinases, kinex_bottom5percentKinases_total = bottom5percentKinasesFraction(kinase_activity_kinex_df, downregulation_df, metadata_experiments_set)
    gsea_bottom5percentKinases, gsea_bottom5percentKinases_total = bottom5percentKinasesFraction(kinase_activity_gsea_df, downregulation_df, metadata_experiments_set)
    ##########
    
    """
    phosx_top5percentKinases, phosx_top5percentKinases_total, phosx_top5percentKinases/phosx_top5percentKinases_total
    kinex_top5percentKinases, kinex_top5percentKinases_total, kinex_top5percentKinases/kinex_top5percentKinases_total
    gsea_top5percentKinases, gsea_top5percentKinases_total, gsea_top5percentKinases/gsea_top5percentKinases_total
    phosx_bottom5percentKinases, phosx_bottom5percentKinases_total, phosx_bottom5percentKinases/phosx_bottom5percentKinases_total
    kinex_bottom5percentKinases, kinex_bottom5percentKinases_total, kinex_bottom5percentKinases/kinex_bottom5percentKinases_total
    gsea_bottom5percentKinases, gsea_bottom5percentKinases_total, gsea_bottom5percentKinases/gsea_bottom5percentKinases_total
    """
     
    phosx_kinases_total = phosx_top5percentKinases_total + phosx_bottom5percentKinases_total
    kinex_kinases_total = kinex_top5percentKinases_total + kinex_bottom5percentKinases_total
    gsea_kinases_total = gsea_top5percentKinases_total + gsea_bottom5percentKinases_total
    
    tp_percentage_df = pd.DataFrame({
        'Upregulation':
            [
                phosx_top5percentKinases/phosx_kinases_total,
                kinex_top5percentKinases/kinex_kinases_total,
                gsea_top5percentKinases/gsea_kinases_total,
            ],
        'Downregulation':
            [
                phosx_bottom5percentKinases/phosx_kinases_total,
                kinex_bottom5percentKinases/kinex_kinases_total,
                gsea_bottom5percentKinases/gsea_kinases_total
            ],
        'Method': 
            [
                'PhosX',
                'Kinex',
                'GSEApy'
            ],
    })
    
    plt.clf()
    ax = tp_percentage_df.set_index('Method').plot(kind='bar', stacked=True, color=['red','blue'], alpha=0.5, figsize=(2.5, 3))
    plt.xticks(rotation=0)
    plt.legend(frameon=False)
    plt.ylim([0.0, 1])
    ax.set_ylabel("Fraction of regulated kinases")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}topNpercentKinases_stacked_barplot.pdf')
    ##########
    
     
    violinplots_upreg_df = pd.DataFrame()
    violinplots_downreg_df = pd.DataFrame()
    violinplots_joined_df = pd.DataFrame()
    
    # PhosX evaluation
    upreg_auc_list = []
    upreg_apr_list = []
    downreg_auc_list = []
    downreg_apr_list = []
    upregulation_phosx_positive_indexes_list = list(set(kinase_activity_phosx_df.index).intersection(set(upregulation_df.index)))
    upregulation_phosx_negative_indexes_list = [x for x in upregulation_negative_indexes_list if x in kinase_activity_phosx_df.index]
    downregulation_phosx_positive_indexes_list = list(set(kinase_activity_phosx_df.index).intersection(set(downregulation_df.index)))
    downregulation_phosx_negative_indexes_list = [x for x in downregulation_negative_indexes_list if x in kinase_activity_phosx_df.index]
    for i in range(100):
        # randomly sample neg ex in equal number as the pos ex, for up- and down- regulation separately
        upreg_neg_idx_list = sample(upregulation_phosx_negative_indexes_list, len(upregulation_phosx_positive_indexes_list))
        downreg_neg_idx_list = sample(downregulation_phosx_negative_indexes_list, len(downregulation_phosx_positive_indexes_list))
        
        # build the evaluation datasets for up- and down- regulation separately
        upreg_df = kinase_activity_phosx_df.loc[upregulation_phosx_positive_indexes_list+upreg_neg_idx_list]
        downreg_df = kinase_activity_phosx_df.loc[downregulation_phosx_positive_indexes_list+downreg_neg_idx_list]
        
        upreg_df = upreg_df.merge(upregulation_df, how='left')
        upreg_df['Regulation'].loc[upreg_df['Regulation'].isna()] = 0
        
        upreg_df.index = upreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
            axis=1
        )
        
        downreg_df = downreg_df.merge(downregulation_df, how='left')
        downreg_df['Regulation'].loc[downreg_df['Regulation'].isna()] = 0
        downreg_df['Regulation'].loc[downreg_df['Regulation'] == -1] = 1
        
        downreg_df.index = downreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
            axis=1
        )
        
        # upregulation
        fpr, tpr, thresholds = roc_curve(upreg_df['Regulation'], upreg_df['Kinase activity change'])
        upreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(upreg_df['Regulation'], upreg_df['Kinase activity change'])
        upreg_apr_list.append(auc(recall, precision))
        
        #downregulation
        fpr, tpr, thresholds = roc_curve(downreg_df['Regulation'], -downreg_df['Kinase activity change'])
        downreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(downreg_df['Regulation'], -downreg_df['Kinase activity change'])
        downreg_apr_list.append(auc(recall, precision))
        
    n_upreg_phosx = len(upreg_df.loc[upreg_df["Regulation"]==1])
    n_downreg_phosx = len(downreg_df.loc[downreg_df["Regulation"]==1])
    
    data = pd.DataFrame({'AUROC':upreg_auc_list, 'AUPR':upreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['PhosX' for i in range(len(data))]
    
    violinplots_upreg_df = pd.concat([violinplots_upreg_df, data])
    
    data = pd.DataFrame({'AUROC':downreg_auc_list, 'AUPR':downreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['PhosX' for i in range(len(data))]
    
    violinplots_downreg_df = pd.concat([violinplots_downreg_df, data])
    
    data = pd.DataFrame({'AUROC':upreg_auc_list+downreg_auc_list, 'AUPR':upreg_apr_list+downreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['PhosX' for i in range(len(data))]
    
    violinplots_joined_df = pd.concat([violinplots_joined_df, data])
    ##########
    
    # Kinex evaluation
    upreg_auc_list = []
    upreg_apr_list = []
    downreg_auc_list = []
    downreg_apr_list = []
    upregulation_kinex_positive_indexes_list = list(set(kinase_activity_kinex_df.index).intersection(set(upregulation_df.index)))
    upregulation_kinex_negative_indexes_list = [x for x in upregulation_negative_indexes_list if x in kinase_activity_kinex_df.index]
    downregulation_kinex_positive_indexes_list = list(set(kinase_activity_kinex_df.index).intersection(set(downregulation_df.index)))
    downregulation_kinex_negative_indexes_list = [x for x in downregulation_negative_indexes_list if x in kinase_activity_kinex_df.index]
    for i in range(100):
        # randomly sample neg ex in equal number as the pos ex, for up- and down- regulation separately
        upreg_neg_idx_list = sample(upregulation_kinex_negative_indexes_list, len(upregulation_kinex_positive_indexes_list))
        downreg_neg_idx_list = sample(downregulation_kinex_negative_indexes_list, len(downregulation_kinex_positive_indexes_list))
        
        # build the evaluation datasets for up- and down- regulation separately
        upreg_df = kinase_activity_kinex_df.loc[upregulation_kinex_positive_indexes_list+upreg_neg_idx_list]
        downreg_df = kinase_activity_kinex_df.loc[downregulation_kinex_positive_indexes_list+downreg_neg_idx_list]
        
        upreg_df = upreg_df.merge(upregulation_df, how='left')
        upreg_df['Regulation'].loc[upreg_df['Regulation'].isna()] = 0
        
        upreg_df.index = upreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
            axis=1
        )
        
        downreg_df = downreg_df.merge(downregulation_df, how='left')
        downreg_df['Regulation'].loc[downreg_df['Regulation'].isna()] = 0
        downreg_df['Regulation'].loc[downreg_df['Regulation'] == -1] = 1
        
        downreg_df.index = downreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
            axis=1
        )
        
        # upregulation
        fpr, tpr, thresholds = roc_curve(upreg_df['Regulation'], upreg_df['Kinase activity change'])
        upreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(upreg_df['Regulation'], upreg_df['Kinase activity change'])
        upreg_apr_list.append(auc(recall, precision))
        
        #downregulation
        fpr, tpr, thresholds = roc_curve(downreg_df['Regulation'], -downreg_df['Kinase activity change'])
        downreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(downreg_df['Regulation'], -downreg_df['Kinase activity change'])
        downreg_apr_list.append(auc(recall, precision))
    
    n_upreg_kinex = len(upreg_df.loc[upreg_df["Regulation"]==1])
    n_downreg_kinex = len(downreg_df.loc[downreg_df["Regulation"]==1])
    
    data = pd.DataFrame({'AUROC':upreg_auc_list, 'AUPR':upreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['Kinex' for i in range(len(data))]
    
    violinplots_upreg_df = pd.concat([violinplots_upreg_df, data])
    
    data = pd.DataFrame({'AUROC':downreg_auc_list, 'AUPR':downreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['Kinex' for i in range(len(data))]
    
    violinplots_downreg_df = pd.concat([violinplots_downreg_df, data])
    
    data = pd.DataFrame({'AUROC':upreg_auc_list+downreg_auc_list, 'AUPR':upreg_apr_list+downreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['Kinex' for i in range(len(data))]
    
    violinplots_joined_df = pd.concat([violinplots_joined_df, data])
    ##########
    
      
    # GSEApy evaluation
    upreg_auc_list = []
    upreg_apr_list = []
    downreg_auc_list = []
    downreg_apr_list = []
    upregulation_gsea_positive_indexes_list = list(set(kinase_activity_gsea_df.index).intersection(set(upregulation_df.index)))
    upregulation_gsea_negative_indexes_list = [x for x in upregulation_negative_indexes_list if x in kinase_activity_gsea_df.index]
    downregulation_gsea_positive_indexes_list = list(set(kinase_activity_gsea_df.index).intersection(set(downregulation_df.index)))
    downregulation_gsea_negative_indexes_list = [x for x in downregulation_negative_indexes_list if x in kinase_activity_gsea_df.index]
    for i in range(100):
        # randomly sample neg ex in equal number as the pos ex, for up- and down- regulation separately
        upreg_neg_idx_list = sample(upregulation_gsea_negative_indexes_list, len(upregulation_gsea_positive_indexes_list))
        downreg_neg_idx_list = sample(downregulation_gsea_negative_indexes_list, len(downregulation_gsea_positive_indexes_list))
        
        # build the evaluation datasets for up- and down- regulation separately
        upreg_df = kinase_activity_gsea_df.loc[upregulation_gsea_positive_indexes_list+upreg_neg_idx_list]
        downreg_df = kinase_activity_gsea_df.loc[downregulation_gsea_positive_indexes_list+downreg_neg_idx_list]
        
        upreg_df = upreg_df.merge(upregulation_df, how='left')
        upreg_df['Regulation'].loc[upreg_df['Regulation'].isna()] = 0
        
        upreg_df.index = upreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
            axis=1
        )
        
        downreg_df = downreg_df.merge(downregulation_df, how='left')
        downreg_df['Regulation'].loc[downreg_df['Regulation'].isna()] = 0
        downreg_df['Regulation'].loc[downreg_df['Regulation'] == -1] = 1
        
        downreg_df.index = downreg_df.apply(
            lambda x: f'{x["Experiment"]}__{x["Kinase"]}',
            axis=1
        )
        
        # upregulation
        fpr, tpr, thresholds = roc_curve(upreg_df['Regulation'], upreg_df['Kinase activity change'])
        upreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(upreg_df['Regulation'], upreg_df['Kinase activity change'])
        upreg_apr_list.append(auc(recall, precision))
        
        #downregulation
        fpr, tpr, thresholds = roc_curve(downreg_df['Regulation'], -downreg_df['Kinase activity change'])
        downreg_auc_list.append(auc(fpr, tpr))
        precision, recall, thresholds = precision_recall_curve(downreg_df['Regulation'], -downreg_df['Kinase activity change'])
        downreg_apr_list.append(auc(recall, precision))
        
    n_upreg_gseapy = len(upreg_df.loc[upreg_df["Regulation"]==1])
    n_downreg_gseapy = len(downreg_df.loc[downreg_df["Regulation"]==1])
    
    data = pd.DataFrame({'AUROC':upreg_auc_list, 'AUPR':upreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['GSEApy' for i in range(len(data))]
    
    violinplots_upreg_df = pd.concat([violinplots_upreg_df, data])
    
    data = pd.DataFrame({'AUROC':downreg_auc_list, 'AUPR':downreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['GSEApy' for i in range(len(data))]
    
    violinplots_downreg_df = pd.concat([violinplots_downreg_df, data])
    
    data = pd.DataFrame({'AUROC':upreg_auc_list+downreg_auc_list, 'AUPR':upreg_apr_list+downreg_apr_list})
    data = data.melt()
    data.columns = ['Metric', 'Value']
    data['Method'] = ['GSEApy' for i in range(len(data))]
    
    violinplots_joined_df = pd.concat([violinplots_joined_df, data])
    ##########
    
    plt.clf()
    plt.figure(figsize=[3.5,4])
    ax = sns.violinplot(
        data=violinplots_upreg_df,
        x='Metric',
        y='Value',
        hue='Method',
        cut=0
    )
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    plt.axhline(y=0.5, color='grey', linestyle='--', linewidth=2)
    ax.set_title(f'Positive upregulated examples\nPhosX: {n_upreg_phosx}; Kinex: {n_upreg_kinex}; GSEApy: {n_upreg_gseapy}')
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_upreg_auc_prc_violinplots_w_title.pdf')
    
    plt.clf()
    plt.figure(figsize=[3.5,4])
    ax = sns.violinplot(
        data=violinplots_downreg_df,
        x='Metric',
        y='Value',
        hue='Method',
        cut=0
    )
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    plt.axhline(y=0.5, color='grey', linestyle='--', linewidth=2)
    ax.set_title(f'Positive downregulated examples\nPhosX: {n_downreg_phosx}; Kinex: {n_downreg_kinex}; GSEApy: {n_downreg_gseapy}')
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_downreg_auc_prc_violinplots_w_title.pdf')
    
    plt.clf()
    plt.figure(figsize=[3.5,4])
    ax = sns.violinplot(
        data=violinplots_joined_df,
        x='Metric',
        y='Value',
        hue='Method',
        cut=0
    )
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    plt.axhline(y=0.5, color='grey', linestyle='--', linewidth=2)
    ax.set_title(f'Joint examples\nPhosX: {n_downreg_phosx + n_upreg_phosx}; Kinex: {n_downreg_kinex + n_upreg_kinex}; GSEApy: {n_downreg_gseapy + n_upreg_gseapy}')
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_joined_auc_prc_violinplots_w_title.pdf')
    
    plt.clf()
    plt.figure(figsize=[3.5,3])
    ax = sns.violinplot(
        data=violinplots_upreg_df,
        x='Metric',
        y='Value',
        hue='Method',
        cut=0
    )
    plt.ylim([0.0, 1.0])
    plt.axhline(y=0.5, color='grey', linestyle='--', linewidth=2)
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_upreg_auc_prc_violinplots.pdf')
    
    plt.clf()
    plt.figure(figsize=[3.5,3])
    ax = sns.violinplot(
        data=violinplots_downreg_df,
        x='Metric',
        y='Value',
        hue='Method',
        cut=0
    )
    plt.ylim([0.0, 1.0])
    plt.axhline(y=0.5, color='grey', linestyle='--', linewidth=2)
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_downreg_auc_prc_violinplots.pdf')
    
    plt.clf()
    plt.figure(figsize=[3.5,3])
    ax = sns.violinplot(
        data=violinplots_joined_df,
        x='Metric',
        y='Value',
        hue='Method',
        cut=0
    )
    plt.ylim([0.0, 1.0])
    plt.axhline(y=0.5, color='grey', linestyle='--', linewidth=2)
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_joined_auc_prc_violinplots.pdf')
    ##########    

    # Upregulated and Downregulated subset comparisons
    kinase_activity_phosx_df = kinase_activity_phosx_df.rename(columns={'Kinase activity change':'PhosX Activity Score'})
    kinase_activity_gsea_df = kinase_activity_gsea_df.rename(columns={'Kinase activity change':'GSEApy Activity Score'})
    kinase_activity_kinex_df = kinase_activity_kinex_df.rename(columns={'Kinase activity change':'Kinex Activity Score'})
    
    regulated_df = metadata_df.merge(kinase_activity_phosx_df, on=['Experiment','Kinase'],how='left').merge(kinase_activity_gsea_df, on=['Experiment','Kinase'],how='left').merge(kinase_activity_kinex_df, on=['Experiment','Kinase'],how='left')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    plt.ylim([0,30])
    sns.histplot(regulated_df['PhosX Activity Score'].loc[regulated_df['Regulation']==1], color='red', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Upregulated')
    sns.histplot(regulated_df['PhosX Activity Score'].loc[regulated_df['Regulation']==-1], color='blue', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Downregulated')
    plt.legend(loc='upper left', frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_PhosX_score_regulated_kinases.pdf')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    plt.ylim([0,30])
    sns.histplot(regulated_df['GSEApy Activity Score'].loc[regulated_df['Regulation']==1], color='red', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Upregulated')
    sns.histplot(regulated_df['GSEApy Activity Score'].loc[regulated_df['Regulation']==-1], color='blue', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Downregulated')
    plt.legend(loc='upper left', frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_GSEApy_score_regulated_kinases.pdf')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    plt.ylim([0,30])
    sns.histplot(regulated_df['Kinex Activity Score'].loc[regulated_df['Regulation']==1], color='red', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Upregulated')
    sns.histplot(regulated_df['Kinex Activity Score'].loc[regulated_df['Regulation']==-1], color='blue', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Downregulated')
    plt.legend(loc='upper left', frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_Kinex_score_regulated_kinases.pdf')
    
   
    regulated_df = regulated_df.dropna()
    
    n_upregulated = len(regulated_df.loc[regulated_df["Regulation"]==1])
    n_downregulated = len(regulated_df.loc[regulated_df["Regulation"]==-1])
    
    phosx_fpr, phosx_tpr, phosx_roc_auc = compute_roc(regulated_df['Regulation'], regulated_df['PhosX Activity Score'])
    gsea_fpr, gsea_tpr, gsea_roc_auc = compute_roc(regulated_df['Regulation'], regulated_df['GSEApy Activity Score'])
    kinex_fpr, kinex_tpr, kinex_roc_auc = compute_roc(regulated_df['Regulation'], regulated_df['Kinex Activity Score'])
    
    class_imbalance = n_upregulated / (n_upregulated + n_downregulated)
    if class_imbalance < 0.5:
        class_imbalance = 1 / class_imbalance
    
    plt.clf()
    plt.figure(figsize=(4, 4))
    #plt.plot(phosx_fpr, phosx_tpr, linestyle='-', color='grey', lw=2, label=f'PhosX (AUC = {phosx_roc_auc:.2f})')
    #plt.plot(gsea_fpr, gsea_tpr, linestyle='--', color='grey', lw=2, label=f'GSEApy (AUC = {gsea_roc_auc:.2f})')
    #plt.plot(kinex_fpr, kinex_tpr, linestyle=':', color='grey', lw=2, label=f'Kinex (AUC = {kinex_roc_auc:.2f})')
    plt.plot(phosx_fpr, phosx_tpr, lw=2, label=f'PhosX (AUC = {phosx_roc_auc:.2f})')
    plt.plot(gsea_fpr, gsea_tpr, lw=2, label=f'GSEApy (AUC = {gsea_roc_auc:.2f})')
    plt.plot(kinex_fpr, kinex_tpr, lw=2, label=f'Kinex (AUC = {kinex_roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='black', lw=0.5, linestyle='-')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.legend(loc='lower right', frameon=False)
    plt.title(f'Upregulated examples: {n_upregulated}\nDownregulated examples: {n_downregulated}')
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_regulated_kinases_ROC_w_title.pdf')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    #plt.plot(phosx_fpr, phosx_tpr, linestyle='-', color='grey', lw=2, label=f'PhosX (AUC = {phosx_roc_auc:.2f})')
    #plt.plot(gsea_fpr, gsea_tpr, linestyle='--', color='grey', lw=2, label=f'GSEApy (AUC = {gsea_roc_auc:.2f})')
    #plt.plot(kinex_fpr, kinex_tpr, linestyle=':', color='grey', lw=2, label=f'Kinex (AUC = {kinex_roc_auc:.2f})')
    plt.plot(phosx_fpr, phosx_tpr, lw=2, label=f'PhosX (AUC = {phosx_roc_auc:.2f})')
    plt.plot(gsea_fpr, gsea_tpr, lw=2, label=f'GSEApy (AUC = {gsea_roc_auc:.2f})')
    plt.plot(kinex_fpr, kinex_tpr, lw=2, label=f'Kinex (AUC = {kinex_roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='black', lw=0.5, linestyle='-')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.legend(loc='lower right', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_regulated_kinases_ROC.pdf')

    phosx_recall, phosx_precision, phosx_pr_auc = compute_pr(regulated_df['Regulation'], regulated_df['PhosX Activity Score'])
    gsea_recall, gsea_precision, gsea_pr_auc = compute_pr(regulated_df['Regulation'], regulated_df['GSEApy Activity Score'])
    kinex_recall, kinex_precision, kinex_pr_auc = compute_pr(regulated_df['Regulation'], regulated_df['Kinex Activity Score'])

    plt.clf()
    plt.figure(figsize=(4, 4))
    #plt.plot(phosx_recall, phosx_precision, linestyle='-', color='grey', lw=2, label=f'PhosX (AUC = {phosx_pr_auc:.2f})')
    #plt.plot(gsea_recall, gsea_precision, linestyle='--', color='grey', lw=2, label=f'GSEApy (AUC = {gsea_pr_auc:.2f})')
    #plt.plot(kinex_recall, kinex_precision, linestyle=':', color='grey', lw=2, label=f'Kinex (AUC = {kinex_pr_auc:.2f})')
    plt.plot(phosx_recall, phosx_precision, lw=2, label=f'PhosX (AUC = {phosx_pr_auc:.2f})')
    plt.plot(gsea_recall, gsea_precision, lw=2, label=f'GSEApy (AUC = {gsea_pr_auc:.2f})')
    plt.plot(kinex_recall, kinex_precision, lw=2, label=f'Kinex (AUC = {kinex_pr_auc:.2f})')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.axhline(y=class_imbalance, lw=0.5, color='black', linestyle='-', )
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='lower left', frameon=False)
    plt.title(f'Upregulated examples: {n_upregulated}\nDownregulated examples: {n_downregulated}')
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_regulated_kinases_PR_w_title.pdf')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    #plt.plot(phosx_recall, phosx_precision, linestyle='-', color='grey', lw=2, label=f'PhosX (AUC = {phosx_pr_auc:.2f})')
    #plt.plot(gsea_recall, gsea_precision, linestyle='--', color='grey', lw=2, label=f'GSEApy (AUC = {gsea_pr_auc:.2f})')
    #plt.plot(kinex_recall, kinex_precision, linestyle=':', color='grey', lw=2, label=f'Kinex (AUC = {kinex_pr_auc:.2f})')
    plt.plot(phosx_recall, phosx_precision, lw=2, label=f'PhosX (AUC = {phosx_pr_auc:.2f})')
    plt.plot(gsea_recall, gsea_precision, lw=2, label=f'GSEApy (AUC = {gsea_pr_auc:.2f})')
    plt.plot(kinex_recall, kinex_precision, lw=2, label=f'Kinex (AUC = {kinex_pr_auc:.2f})')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.axhline(y=class_imbalance, lw=0.5, color='black', linestyle='-', )
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='lower left', frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_regulated_kinases_PR.pdf')
    
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    plt.ylim([0,30])
    sns.histplot(regulated_df['PhosX Activity Score'].loc[regulated_df['Regulation']==1], color='red', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Upregulated')
    sns.histplot(regulated_df['PhosX Activity Score'].loc[regulated_df['Regulation']==-1], color='blue', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Downregulated')
    plt.legend(loc='upper left', frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_PhosX_score_regulated_kinases_intersection.pdf')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    plt.ylim([0,30])
    sns.histplot(regulated_df['GSEApy Activity Score'].loc[regulated_df['Regulation']==1], color='red', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Upregulated')
    sns.histplot(regulated_df['GSEApy Activity Score'].loc[regulated_df['Regulation']==-1], color='blue', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Downregulated')
    plt.legend(loc='upper left', frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_GSEApy_score_regulated_kinases_intersection.pdf')
    
    plt.clf()
    plt.figure(figsize=(3.5, 3))
    plt.ylim([0,30])
    sns.histplot(regulated_df['Kinex Activity Score'].loc[regulated_df['Regulation']==1], color='red', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Upregulated')
    sns.histplot(regulated_df['Kinex Activity Score'].loc[regulated_df['Regulation']==-1], color='blue', alpha=0.5, binrange=[0,1], binwidth=0.05, legend='Downregulated')
    plt.legend(loc='upper left', frameon=False, labels=["Upregulated", "Downregulated"])
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{out_prefix}cptac_Kinex_score_regulated_kinases_intersection.pdf')
 
    
if __name__ == '__main__':
    main()