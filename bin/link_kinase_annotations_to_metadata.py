#!/usr/bin/env python3

import pandas as pd
import numpy as np
import h5py
import sys


"""
dict_tsv = "input/gene_synomym2gene_name_dict.tsv"
metadata_tsv = 'input/metadata.tsv'
kinase_metadata_annotated_h5 = 'input/kinase_metadata.h5'
"""


def hdf5_to_dict(hdf5_file_path):
    """
    Convert an HDF5 file's contents into a nested dictionary recursively,
    with special handling for string datasets.

    Args:
        hdf5_file_path (str): Path to the HDF5 file

    Returns:
        dict: Nested dictionary containing the HDF5 file's structure and data
    """
    def _convert_group(h5_group):
        result = {}
        for key, item in h5_group.items():
            # If item is a group, recursively convert it
            if isinstance(item, h5py.Group):
                result[key] = _convert_group(item)
            # If item is a dataset
            elif isinstance(item, h5py.Dataset):
                # Check if the dataset contains string type
                if h5py.check_string_dtype(item.dtype):
                    # Handle scalar string
                    if item.shape == ():
                        value = item[()]
                        result[key] = (
                            value.decode("utf-8")
                            if isinstance(value, bytes)
                            else str(value)
                        )
                    # Handle array of strings
                    else:
                        value = item[:]
                        result[key] = [
                            v.decode("utf-8") if isinstance(v, bytes) else str(v)
                            for v in value
                        ]
                else:
                    # Handle non-string scalar datasets
                    if item.shape == ():
                        result[key] = item[()]
                    # Handle non-string array datasets
                    else:
                        result[key] = item[:]
                    # Convert numpy types to Python native types if possible
                    if hasattr(result[key], "tolist"):
                        result[key] = result[key].tolist()
        return result
    # Open the HDF5 file and convert its contents
    try:
        with h5py.File(hdf5_file_path, "r") as f:
            return _convert_group(f)
    except Exception as e:
        raise Exception(f"Error reading HDF5 file: {str(e)}")


def map_kinase_annotations_to_metadata(kinase_metadata_dict, meta_df, dict_df, default_specificity=np.nan, default_family=np.nan):
    meta_df["Family"] = default_family
    meta_df["Family"] = meta_df["Family"].astype("object")
    meta_df["Specificity"] = default_specificity
    meta_df["Specificity"] = meta_df["Specificity"].astype("object")
    kinase_names = list(meta_df["Kinase"].unique())
    meta_kinase_names = list(kinase_metadata_dict["specificity"].keys())
    for kinase_name in kinase_names:
        meta_kinase_name = None
        if kinase_name in meta_kinase_names:
            meta_kinase_name = kinase_name
        else:
            tr_list = list(dict_df.loc[dict_df["Gene name"]==kinase_name, "Gene synonym"].values)
            tr_list = tr_list + list(dict_df.loc[dict_df["Gene synonym"]==kinase_name, "Gene name"].values)
            for tr_str in tr_list:
                if tr_str in meta_kinase_names:
                    meta_kinase_name = tr_str
                    break
        if meta_kinase_name is not None:
            meta_df.loc[meta_df["Kinase"]==kinase_name, "Specificity"] = kinase_metadata_dict["specificity"][meta_kinase_name]
            meta_df.loc[meta_df["Kinase"]==kinase_name, "Family"] = kinase_metadata_dict["family"][meta_kinase_name]
    return meta_df


def main():
    metadata_tsv = sys.argv[1]
    dict_tsv = sys.argv[2]
    kinase_metadata_annotated_h5 = sys.argv[3]

    dict_df = pd.read_csv(dict_tsv, sep="\t", header=None)
    dict_df.columns = ["Gene synonym", "UniProt AC", "Gene name"]

    metadata_df = pd.read_csv(metadata_tsv, sep="\t", header=None)
    metadata_df.columns = ["Experiment", "Kinase", "Regulation"]

    kinase_metadata_dict = hdf5_to_dict(kinase_metadata_annotated_h5)

    metadata_annotated_df = map_kinase_annotations_to_metadata(
        kinase_metadata_dict,
        metadata_df,
        dict_df
    )

    print(metadata_annotated_df.to_csv(sep="\t", index=False))
    

if __name__ == "__main__":
    main()
