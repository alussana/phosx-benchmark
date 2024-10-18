#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import pickle
import os
from kstar import helpers, mapping, calculate


def max_abs_value(series: pd.Series):
    # Drop NaN values and get the absolute values
    abs_series = series.dropna().abs()
    # If the series is empty after dropping NaN, return None
    if abs_series.empty:
        return None
    # Return the original element with the maximum absolute value
    return series.loc[abs_series.idxmax()]


def main():

    rnk = sys.argv[1]
    odir = sys.argv[2]
    logName = sys.argv[3]
    n_proc = int(sys.argv[4])
    network_dir = sys.argv[5]
    fc_threshold = float(sys.argv[6])

    df = pd.read_csv(rnk, sep="\t", header=None)
    df.columns = ["query_accession", "mod_sites", "peptide", "data:log2fc"]

    if not os.path.exists(f"{odir}/MAPPED_DATA"):
        os.mkdir(f"{odir}/MAPPED_DATA")

    # intialize logger
    mapping_log = helpers.get_logger(
        f"mapping_{logName}", f"{odir}/MAPPED_DATA/mapping_{logName}.log"
    )

    # map dataset and record process in the logger
    mapDict = {"peptide": "peptide", "accession_id": "query_accession"}
    exp_mapper = mapping.ExperimentMapper(
        experiment=df, columns=mapDict, logger=mapping_log
    )

    exp_mapper.experiment.to_csv(
        f"{odir}/MAPPED_DATA/{logName}_mapped.tsv", sep="\t", index=False
    )

    ##########################################
    # compute activities for Y kinases (upreg)

    # determine thresholds
    logName_new = logName + "_Y_upereg"
    # intialize log file
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )

    # load mapped data, if necessary
    experiment = pd.read_csv(
        f"{odir}/MAPPED_DATA/{logName}_mapped.tsv", sep="\t", index_col=0
    )
    # experiment = exp_mapper.experiment

    # Define parameters
    data_columns = None
    agg = "mean"
    greater = True  # for up-regulation
    threshold = fc_threshold
    evidence_size = 100  # overrides threshold

    # intialize KinaseActivity class object
    # kinact = calculate.KinaseActivity(experiment, activity_log,data_columns = data_columns, phospho_type='Y')
    # convert evidence into binary evidence based on the provided threshold
    # evidence_sizes = kinact.test_threshold(agg = agg, threshold = threshold,  greater = greater, return_evidence_sizes = True)

    # Calculate Statistical Enrichment (Hypergeometric p-values)
    phospho_types = ["Y"]
    # Load the pickles containing the 50 pruned networks for tyrosine kinases
    networks = {}
    networks["Y"] = pickle.load(open(f"{network_dir}/kstar_network_Y.p", "rb"))
    # Create activity log: if already did this, ignore.
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )
    kinact_dict = calculate.enrichment_analysis(
        experiment,
        activity_log,
        networks,
        phospho_types=phospho_types,
        agg=agg,
        threshold=threshold,
        evidence_size=evidence_size,
        greater=greater,
        PROCESSES=1,
    )
    # Generate random datasets, run kinase activity on random datasets, normalize original analysis
    # Set the number of random experiments
    num_random_experiments = 150
    # Generate random experiments
    calculate.randomized_analysis(
        kinact_dict, activity_log, num_random_experiments, PROCESSES=n_proc
    )
    # Calculate Mann Whitney Significance
    calculate.Mann_Whitney_analysis(
        kinact_dict, activity_log, number_sig_trials=100, PROCESSES=n_proc
    )
    # export results to disk
    calculate.save_kstar_slim(kinact_dict, logName, odir)
    kinact_dict["Y"].activities_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_activities.tsv", sep="\t"
    )
    kinact_dict["Y"].fpr_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_fpr.tsv", sep="\t"
    )
    # replace zeroes with tiny numbers
    kinact_dict["Y"].fpr_mann_whitney.loc[
        kinact_dict["Y"].fpr_mann_whitney["data:log2fc"] == 0
    ] = np.nextafter(np.float32(0), np.float32(1))
    activity_Y_upreg = -np.log2(kinact_dict["Y"].fpr_mann_whitney)
    activity_Y_upreg.columns = ["activity_Y_upreg"]
    # print(activity_Y_upreg.to_csv(sep="\t", header=False, index=False))

    ##########################################
    # compute activities for S/T kinases (upreg)

    # determine thresholds
    logName_new = logName + "_ST_upereg"
    # intialize log file
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )

    # load mapped data, if necessary
    experiment = pd.read_csv(
        f"{odir}/MAPPED_DATA/{logName}_mapped.tsv", sep="\t", index_col=0
    )
    # experiment = exp_mapper.experiment

    # Define parameters
    data_columns = None
    agg = "mean"
    greater = True  # for up-regulation
    threshold = fc_threshold
    evidence_size = 100  # overrides threshold

    # intialize KinaseActivity class object
    # kinact = calculate.KinaseActivity(experiment, activity_log,data_columns = data_columns, phospho_type='ST')
    # convert evidence into binary evidence based on the provided threshold
    # evidence_sizes = kinact.test_threshold(agg = agg, threshold = threshold,  greater = greater, return_evidence_sizes = True)

    # Calculate Statistical Enrichment (Hypergeometric p-values)
    phospho_types = ["ST"]
    # Load the pickles containing the 50 pruned networks for tyrosine kinases
    networks = {}
    networks["ST"] = pickle.load(open(f"{network_dir}/kstar_network_ST.p", "rb"))
    # Create activity log: if already did this, ignore.
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )
    kinact_dict = calculate.enrichment_analysis(
        experiment,
        activity_log,
        networks,
        phospho_types=phospho_types,
        agg=agg,
        threshold=threshold,
        evidence_size=evidence_size,
        greater=greater,
        PROCESSES=1,
    )
    # Generate random datasets, run kinase activity on random datasets, normalize original analysis
    # Set the number of random experiments
    num_random_experiments = 150
    # Generate random experiments
    calculate.randomized_analysis(
        kinact_dict, activity_log, num_random_experiments, PROCESSES=n_proc
    )
    # Calculate Mann Whitney Significance
    calculate.Mann_Whitney_analysis(
        kinact_dict, activity_log, number_sig_trials=100, PROCESSES=n_proc
    )
    # export results to disk
    calculate.save_kstar_slim(kinact_dict, logName_new, odir)
    kinact_dict["ST"].activities_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_activities.tsv", sep="\t"
    )
    kinact_dict["ST"].fpr_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_fpr.tsv", sep="\t"
    )
    # replace zeroes with tiny numbers
    kinact_dict["ST"].fpr_mann_whitney.loc[
        kinact_dict["ST"].fpr_mann_whitney["data:log2fc"] == 0
    ] = np.nextafter(np.float32(0), np.float32(1))
    activity_ST_upreg = -np.log2(kinact_dict["ST"].fpr_mann_whitney)
    activity_ST_upreg.columns = ["activity_ST_upreg"]
    # print(activity_ST_upreg.to_csv(sep="\t", header=False, index=False))

    ############################################
    # compute activities for Y kinases (downreg)

    # determine thresholds
    logName_new = logName + "_Y_downreg"
    # intialize log file
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )

    # load mapped data, if necessary
    experiment = pd.read_csv(
        f"{odir}/MAPPED_DATA/{logName}_mapped.tsv", sep="\t", index_col=0
    )
    # experiment = exp_mapper.experiment

    # Define parameters
    data_columns = None
    agg = "mean"
    greater = False  # for up-regulation
    threshold = -fc_threshold
    evidence_size = 100  # overrides threshold

    # intialize KinaseActivity class object
    # kinact = calculate.KinaseActivity(experiment, activity_log,data_columns = data_columns, phospho_type='Y')
    # convert evidence into binary evidence based on the provided threshold
    # evidence_sizes = kinact.test_threshold(agg = agg, threshold = threshold,  greater = greater, return_evidence_sizes = True)

    # Calculate Statistical Enrichment (Hypergeometric p-values)
    phospho_types = ["Y"]
    # Load the pickles containing the 50 pruned networks for tyrosine kinases
    networks = {}
    networks["Y"] = pickle.load(open(f"{network_dir}/kstar_network_Y.p", "rb"))
    # Create activity log: if already did this, ignore.
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )
    kinact_dict = calculate.enrichment_analysis(
        experiment,
        activity_log,
        networks,
        phospho_types=phospho_types,
        agg=agg,
        threshold=threshold,
        evidence_size=evidence_size,
        greater=greater,
        PROCESSES=1,
    )
    # Generate random datasets, run kinase activity on random datasets, normalize original analysis
    # Set the number of random experiments
    num_random_experiments = 150
    # Generate random experiments
    calculate.randomized_analysis(
        kinact_dict, activity_log, num_random_experiments, PROCESSES=n_proc
    )
    # Calculate Mann Whitney Significance
    calculate.Mann_Whitney_analysis(
        kinact_dict, activity_log, number_sig_trials=100, PROCESSES=n_proc
    )
    # export results to disk
    calculate.save_kstar_slim(kinact_dict, logName, odir)
    kinact_dict["Y"].activities_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_activities.tsv", sep="\t"
    )
    kinact_dict["Y"].fpr_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_fpr.tsv", sep="\t"
    )
    # replace zeroes with tiny numbers
    kinact_dict["Y"].fpr_mann_whitney.loc[
        kinact_dict["Y"].fpr_mann_whitney["data:log2fc"] == 0
    ] = np.nextafter(np.float32(0), np.float32(1))
    activity_Y_downreg = np.log2(kinact_dict["Y"].fpr_mann_whitney)
    activity_Y_downreg.columns = ["activity_Y_downreg"]
    # print(activity_Y_downreg.to_csv(sep="\t", header=False, index=False))

    ##############################################
    # compute activities for S/T kinases (downreg)

    # determine thresholds
    logName_new = logName + "_ST_downreg"
    # intialize log file
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )

    # load mapped data, if necessary
    experiment = pd.read_csv(
        f"{odir}/MAPPED_DATA/{logName}_mapped.tsv", sep="\t", index_col=0
    )
    # experiment = exp_mapper.experiment

    # Define parameters
    data_columns = None
    agg = "mean"
    greater = False  # for up-regulation
    threshold = -fc_threshold
    evidence_size = 100  # overrides threshold

    # intialize KinaseActivity class object
    # kinact = calculate.KinaseActivity(experiment, activity_log,data_columns = data_columns, phospho_type='ST')
    # convert evidence into binary evidence based on the provided threshold
    # evidence_sizes = kinact.test_threshold(agg = agg, threshold = threshold,  greater = greater, return_evidence_sizes = True)

    # Calculate Statistical Enrichment (Hypergeometric p-values)
    phospho_types = ["ST"]
    # Load the pickles containing the 50 pruned networks for tyrosine kinases
    networks = {}
    networks["ST"] = pickle.load(open(f"{network_dir}/kstar_network_ST.p", "rb"))
    # Create activity log: if already did this, ignore.
    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")
    activity_log = helpers.get_logger(
        f"activity_{logName_new}", f"{odir}/RESULTS/activity_{logName_new}.log"
    )
    kinact_dict = calculate.enrichment_analysis(
        experiment,
        activity_log,
        networks,
        phospho_types=phospho_types,
        agg=agg,
        threshold=threshold,
        evidence_size=evidence_size,
        greater=greater,
        PROCESSES=1,
    )
    # Generate random datasets, run kinase activity on random datasets, normalize original analysis
    # Set the number of random experiments
    num_random_experiments = 150
    # Generate random experiments
    calculate.randomized_analysis(
        kinact_dict, activity_log, num_random_experiments, PROCESSES=n_proc
    )
    # Calculate Mann Whitney Significance
    calculate.Mann_Whitney_analysis(
        kinact_dict, activity_log, number_sig_trials=100, PROCESSES=n_proc
    )
    # export results to disk
    calculate.save_kstar_slim(kinact_dict, logName_new, odir)
    kinact_dict["ST"].activities_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_activities.tsv", sep="\t"
    )
    kinact_dict["ST"].fpr_mann_whitney.to_csv(
        f"{logName_new}_mann_whitney_fpr.tsv", sep="\t"
    )
    # replace zeroes with tiny numbers
    kinact_dict["ST"].fpr_mann_whitney.loc[
        kinact_dict["ST"].fpr_mann_whitney["data:log2fc"] == 0
    ] = np.nextafter(np.float32(0), np.float32(1))
    activity_ST_downreg = np.log2(kinact_dict["ST"].fpr_mann_whitney)
    activity_ST_downreg.columns = ["activity_ST_downreg"]
    # print(activity_ST_downreg.to_csv(sep="\t", header=False, index=False))

    activity_merged_df = pd.concat(
        [activity_ST_upreg, activity_ST_downreg, activity_Y_upreg, activity_Y_downreg],
        axis=1,
        join="outer",
    )

    activity_df = activity_merged_df.apply(max_abs_value, axis=1)

    print(activity_df.to_csv(sep="\t", header=False, index=True))


if __name__ == "__main__":
    main()
