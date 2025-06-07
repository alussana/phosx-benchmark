#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns
import sys


plt.rcParams['axes.titlesize'] = 8        
plt.rcParams['axes.labelsize'] = 6        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6 


def plot_kin_freq(df, bar_pdf):
    # Count occurrences of each kinase
    counts = df["Kinase"].value_counts()

    # Create a new dataframe for plotting
    plot_data = pd.DataFrame({
        "Kinase": counts.index,
        "Count": counts.values
    })

    # Define a color palette
    palette = sns.color_palette("viridis", len(plot_data))

    # Create the barplot
    plt.clf()
    plt.figure(figsize=(5, 2.5))
    ax = sns.barplot(x="Kinase", y="Count", data=plot_data, palette=palette)

    # Customize the plot
    plt.xlabel(f"Kinase reported in PhosphoSitePlus (n={len(counts)})")
    plt.ylabel("Count of annotated target phosphosites")
    plt.xticks([])
    sns.despine()

    # Save the plot
    plt.tight_layout()
    plt.savefig(bar_pdf)


def plot_kin_freq_v2(df, bar_pdf):
    # Count the frequency of each kinase
    kinase_counts = df['Kinase'].value_counts()

    # Create the plot
    plt.figure(figsize=(12, 6))

    # Bar plot with rotated labels
    kinase_counts.plot(kind='bar')
    plt.title('Frequency of Different Kinases')
    plt.xlabel('Kinase')
    plt.ylabel('Frequency')
    plt.xticks([])
    plt.tight_layout()

    # Save the plot
    plt.savefig(bar_pdf)


def main():
    kin_phos_tsv = sys.argv[1]
    kin_sub_tsv = sys.argv[2]
    kin_freq_pdf = sys.argv[3]
    phos_freq_pdf = sys.argv[4]
    sub_freq_pdf = sys.argv[5]

    df = pd.read_csv(kin_phos_tsv, sep="\t", header=None)
    df.columns = ["Kinase", "Phosphosite"]

    plot_kin_freq(df, kin_freq_pdf)

    #df = pd.read_csv(kin_sub_tsv, sep="\t")
    #df.columns = ["Kinase", "Substrate"]


if __name__ == "__main__":
    main()
