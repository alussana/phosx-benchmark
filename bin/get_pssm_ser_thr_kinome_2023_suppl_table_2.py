#!/usr/bin/env python3

import numpy as np
import pandas as pd

def get_kinase_pssm(kinase_str:str, matrix_df:pd.DataFrame, densitometry_df:pd.DataFrame):
    aa_list = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T','s','t','y']
    pos_list = [-5, -4, -3, -2, -1, 1, 2, 3, 4]
    kinase_series = matrix_df.loc[kinase_str,]
    densitometry_series = densitometry_df.loc[kinase_str,]
    #field_list = list(kinase_series.index)
    #pssm_values_list = np.array([[0 for aa in aa_list] for pos in pos_list], dtype=float)
    pssm_values_list = np.empty([len(pos_list), len(aa_list)], dtype = float)
    for i in range(len(pos_list)):
        pos = pos_list[i]
        for j in range(len(aa_list)):
            aa = aa_list[j]
            pssm_values_list[i][j] = kinase_series[f'{pos}{aa}']
    pssm = pd.DataFrame(pssm_values_list, index=pos_list, columns=aa_list)
    S_0 = np.empty(len(pos_list), dtype = float)
    T_0 = np.empty(len(pos_list), dtype = float)
    for i in range(len(pos_list)):
        pos = pos_list[i]
        S_0[i] = densitometry_series[f'{pos}S']
        T_0[i] = densitometry_series[f'{pos}T']
    S_S = S_0.sum()
    S_T = T_0.sum()
    S_ctrl = 0.75 * S_S - 0.25 * S_T
    T_ctrl = 0.75 * S_T - 0.25 * S_S
    S_0 = S_ctrl / max(S_ctrl, T_ctrl)
    T_0 = T_ctrl / max(S_ctrl, T_ctrl)
    pos_0_df = pd.DataFrame(index=[0], columns=aa_list)
    pos_0_df['S'] = S_0
    pos_0_df['T'] = T_0
    pos_0_df = pos_0_df.fillna(0)
    pssm = pd.concat([pssm, pos_0_df])
    pssm = pssm.sort_index()    
    #pssm = pssm.apply(lambda x: x / x.sum(), axis=1)
    return pssm

def plot_sequence_logo(pssm_df:pd.DataFrame, file_path:str='logo.svg', ylabel_str:str='', title_str:str='') -> None:
    import logomaker
    import matplotlib.pyplot as plt
    central_residue_df = pd.DataFrame([pssm_df.loc[0,]])
    pssm_df = pssm_df.drop(0, axis=0)
    pssm_median_shifted_df = pssm_df.apply(lambda x: np.log2(x / x.median()), axis=1)
    max_pos_height_float = pssm_median_shifted_df[pssm_median_shifted_df > 0].T.sum().max()
    central_residue_df = central_residue_df * max_pos_height_float / central_residue_df.T.sum().values[0]
    pssm_median_shifted_df = pd.concat([pssm_median_shifted_df, central_residue_df], axis=0)
    pssm_median_shifted_df = pssm_median_shifted_df.sort_index()
    plt.clf()
    # create Logo object
    pssm_logo = logomaker.Logo(
        pssm_median_shifted_df,
        shade_below=0,
        fade_below=0,
        width=0.9,
        color_scheme='chemistry',
        figsize=[4,2],
        flip_below=False,
        stack_order='big_on_top',
        center_values=False,
        show_spines=True,
        #font_name='Arial Rounded MT Bold')
    )
    # style using Logo methods
    pssm_logo.style_spines(visible=False)
    pssm_logo.style_spines(spines=['left', 'bottom'], visible=True)
    pssm_logo.style_xticks(rotation=0, fmt='%d', anchor=0)
    # style using Axes methods
    pssm_logo.ax.set_ylabel(ylabel_str, labelpad=-1)
    pssm_logo.ax.set_title(title_str)
    pssm_logo.ax.xaxis.set_ticks_position('none')
    pssm_logo.ax.xaxis.set_tick_params(pad=-1)
    plt.savefig(file_path)
    plt.close()
    
def save_dict_to_hdf5(dictionary:dict, output_file:str):
    """
    Save a Python dictionary to an HDF5 file.

    Parameters:
        pssm_dict (dict): The Python dictionary to be saved.
        output_file (str): The name of the HDF5 file to be created.

    Returns:
        None
    """
    import h5py
    try:
        with h5py.File(output_file, 'w') as h5file:
            # Traverse the dictionary and save each key-value pair as a dataset in the HDF5 file
            for key, value in dictionary.items():
                h5file.create_dataset(str(key), data=value)
        print("HDF5 file successfully created and data saved.")
    except Exception as e:
        print(f"An error occurred while saving the data: {e}")
        
def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word, ][to_col - 1]
    if len(trs) == 0:
        return(None)
    else:
        return(trs.values[0])

def main():
    norm_scaled_pssm_df = pd.read_excel(
        'input/kinome_2023_suppl_table_2.xlsx',
        sheet_name='ser_thr_all_norm_scaled_matrice',
        index_col=0
    )

    densitometry_df = pd.read_excel(
        'input/kinome_2023_suppl_table_2.xlsx',
        sheet_name='ser_thr_all_raw_matrices',
        index_col=0
    )
    
    """
    norm_pssm_df = pd.read_excel(
        '41586_2022_5575_MOESM4_ESM.xlsx',
        sheet_name='ser_thr_all_norm_matrices',
        index_col=0
    )
    """
    
    dict_file = 'input/gene_synonym_2_gene_name_dict.tsv'
    word_dict = pd.read_csv(dict_file, sep='\t', header=None)
    
    # remove PAK1 --> PKN1 translation 
    # remove PDHK1 --> PDK1 translation
    # see get_ser_thr_kinase_pssm() process in modules/pssm.nf
    word_dict.loc[word_dict[0]=='PAK1'] = None
    word_dict.loc[word_dict[0]=='PDHK1'] = None
    
    # translate gene synonyms into gene names for kinases in norm_scaled_pssm_df.index
    tr_index = []
    for i in range(len(norm_scaled_pssm_df.index)):
        word = norm_scaled_pssm_df.index[i]
        tr_word = translateWord(word, word_dict, from_col=1, to_col=3)
        if tr_word == None:
            tr_word = word
        tr_index.append(tr_word)
    norm_scaled_pssm_df.index = tr_index
    
    # translate gene synonyms into gene names for kinases in densitometry_df.index
    tr_index = []
    for i in range(len(densitometry_df.index)):
        word = densitometry_df.index[i]
        tr_word = translateWord(word, word_dict, from_col=1, to_col=3)
        if tr_word == None:
            tr_word = word
        tr_index.append(tr_word)
    densitometry_df.index = tr_index
    
    kinase_list = list(norm_scaled_pssm_df.index)
    
    pssm_dict = {}
    for kinase_str in kinase_list:
        pssm_dict[kinase_str] = get_kinase_pssm(kinase_str, norm_scaled_pssm_df, densitometry_df)
        plot_sequence_logo(
            pssm_dict[kinase_str],
            file_path=f'pssm/logos/{kinase_str}.svg',
            ylabel_str='log2(Ratio to Median)',
            title_str=kinase_str
        )
        pssm_dict[kinase_str].to_csv(
            path_or_buf=f'pssm/{kinase_str}.tsv',
            sep='\t',
            index=True,
            header=True
        )
    save_dict_to_hdf5(pssm_dict, 'pssm/pssm_dict.h5')

if __name__ == '__main__':
    main()