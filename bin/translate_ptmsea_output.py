#!/usr/bin/env python3

import sys
import pandas as pd


alias_map_txt_gz = sys.argv[1]
alias_map_df = pd.read_csv(alias_map_txt_gz, sep="\t")

Ensembl_HGNC_symbol_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_HGNC_symbol", "alias"
].unique()

Ensembl_HGNC_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_HGNC", "alias"
].unique()

Ensembl_UniProt_list = alias_map_df.loc[
    alias_map_df["source"] == "Ensembl_UniProt", "alias"
].unique()


def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word,][to_col - 1]
    if len(trs) == 0:
        return None
    else:
        return trs.values[0]


def find_translation(
    key,
):
    # The following mappings would collide if not for the ad-hoc flow control that follows:
    # Ser/Thr:
    # HGK     MAP4K4
    # NIK     MAP4K4
    # PDHK1   PDK1
    # PDK1    PDK1
    # Tyr:
    # EPHA3   EPHA3
    # ETK     EPHA3
    if key == "NIK":
        return "MAP3K14"
    elif key == "PDK1":
        return "PDPK1"
    elif key == "ETK":
        return "BMX"
    key_parts = key.split("_")
    key = key_parts[0]
    if (
        key in Ensembl_HGNC_symbol_list
        or key in Ensembl_HGNC_list
        or key in Ensembl_UniProt_list
    ):
        if len(key_parts) > 1:
            new_key = f"{key}_TYR"
        else:
            new_key = key
        return new_key
    key_alias_df = alias_map_df.loc[alias_map_df["alias"] == key]
    string_id_array = key_alias_df["#string_protein_id"].unique()
    if len(string_id_array) == 0:
        return key
    else:
        for string_id_str in string_id_array:
            try:
                map_df = alias_map_df.loc[
                    alias_map_df["#string_protein_id"] == string_id_str
                ]
                new_key = map_df.loc[map_df["source"] == "Ensembl_HGNC_symbol", "alias"]
                if len(new_key) == 0:
                    new_key = map_df.loc[map_df["source"] == "Ensembl_HGNC", "alias"]
                if len(new_key) == 0:
                    new_key = map_df.loc[map_df["source"] == "Ensembl_UniProt", "alias"]
                if len(new_key) == 0:
                    next
                new_key = new_key.values[0]
                if len(key_parts) > 1:
                    new_key = f"{new_key}_TYR"
                break
            except Exception as e:
                print(f"{key}: {e}")
                next
        return new_key


def main():

    target_file = sys.argv[2]
    field_to_translate = int(sys.argv[3])

    print("\tActivity Score")
    with open(target_file, "r") as target_fh:
        for line in target_fh:
            fields = line.removesuffix("\n").split("\t")
            tr_word = find_translation(fields[field_to_translate - 1])
            fields[field_to_translate - 1] = tr_word
            print("\t".join([str(field) for field in fields]))


if __name__ == "__main__":
    main()
