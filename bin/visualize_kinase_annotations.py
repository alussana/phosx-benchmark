#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys


"""
metadata_tsv = 'input/metadata.tsv'
family_bar_pdf = 'datasets/hernandez2017/kinase_family_regulation_counts.pdf'
"""


plt.rcParams['axes.titlesize'] = 8        
plt.rcParams['axes.labelsize'] = 6        
plt.rcParams['xtick.labelsize'] = 6       
plt.rcParams['ytick.labelsize'] = 6       
plt.rcParams['legend.fontsize'] = 6       
plt.rcParams['figure.titlesize'] = 8   
plt.rcParams['font.size'] = 6 


kinase_family_colors = {
    "AGC": "#00BA38",
    "CMGC": "#00C19F",
    "CAMK": "#619CFF",
    "CK1": "#DB72FB",
    "STE": "#FF61C3",
    "TK": "#D39200",
    "TKL": "#93AA00",
    "Atypical": "#00B9E3",
    "Other": "#F8766D",
    # "NA": "#BEBEBE",
}
kinase_families = list(kinase_family_colors.keys())


kinase_specificity_colors = {
    "SerThr": "#ADD8E6",
    "Tyr": "#E6DAA6",
    "Dual": "#000000",
}
kinase_specificities = list(kinase_specificity_colors.keys())


def plot_kinase_family_regulation_counts(df, x_labels, bar_pdf):
    # Group by Family and Regulation, count occurrences, and unstack to get counts for Regulation=1 and -1
    counts = df.groupby(["Family", "Regulation"]).size().unstack(fill_value=0)

    # Reindex to include all families from kinase_families, filling missing with 0
    counts = counts.reindex(x_labels, fill_value=0)

    # Create a new dataframe for plotting
    # Positive counts for Regulation=1, negative counts for Regulation=-1
    plot_data = pd.DataFrame(
        {
            "Family": counts.index,
            "Activation": counts[1],
            "Inhibition": -counts[-1],  # Negative for left side of x-axis
        }
    )

    # Melt the dataframe to long format for seaborn
    plot_data_melted = pd.melt(
        plot_data,
        id_vars="Family",
        value_vars=["Activation", "Inhibition"],
        var_name="Regulation",
        value_name="Count",
    )

    # Define a color palette
    palette = {"Activation": "#FF000080", "Inhibition": "#0000FF80"}

    legend_elements = [
        Patch(facecolor=palette["Activation"], label="Activation"),
        Patch(facecolor=palette["Inhibition"], label="Inhibition"),
    ]

    # Create the barplot
    plt.clf()
    plt.figure(figsize=(4.0, 2.5))
    ax = sns.barplot(
        x="Family",
        y="Count",
        hue="Regulation",
        data=plot_data_melted,
        palette=palette,
        order=x_labels,
    )

    # Set alpha for all bars after creation
    for patch in ax.patches:
        patch.set_alpha(0.5)

    # Customize the plot
    plt.xlabel("Kinase family")
    plt.ylabel("Count")
    plt.axhline(0, color="black", linewidth=0.5)  # Add vertical line at x=0
    # plt.legend(title='Regulation', labels=['Activation', 'Inhibition'], frameon=False)
    plt.legend(handles=legend_elements, title="Regulation", frameon=False)
    # plt.xticks(rotation=45)
    plt.tight_layout()
    sns.despine()

    # Save the plot
    plt.savefig(bar_pdf)


def plot_kinase_specificity_regulation_counts(df, x_labels, bar_pdf):
    # Group by Family and Regulation, count occurrences, and unstack to get counts for Regulation=1 and -1
    counts = df.groupby(["Specificity", "Regulation"]).size().unstack(fill_value=0)

    # Reindex to include all families from kinase_families, filling missing with 0
    counts = counts.reindex(x_labels, fill_value=0)

    # Create a new dataframe for plotting
    # Positive counts for Regulation=1, negative counts for Regulation=-1
    plot_data = pd.DataFrame(
        {
            "Specificity": counts.index,
            "Activation": counts[1],
            "Inhibition": -counts[-1],  # Negative for left side of x-axis
        }
    )

    # Melt the dataframe to long format for seaborn
    plot_data_melted = pd.melt(
        plot_data,
        id_vars="Specificity",
        value_vars=["Activation", "Inhibition"],
        var_name="Regulation",
        value_name="Count",
    )

    # Define a color palette
    palette = {"Activation": "#FF000080", "Inhibition": "#0000FF80"}


    legend_elements = [
        Patch(facecolor=palette["Activation"], label="Activation"),
        Patch(facecolor=palette["Inhibition"], label="Inhibition"),
    ]

    # Create the barplot
    plt.clf()
    plt.figure(figsize=(1.75, 2.5))
    ax = sns.barplot(
        x="Specificity",
        y="Count",
        hue="Regulation",
        data=plot_data_melted,
        palette=palette,
        order=x_labels,
    )

    # Set alpha for all bars after creation
    for patch in ax.patches:
        patch.set_alpha(0.5)

    # Customize the plot
    plt.xlabel("Kinase specificity")
    plt.ylabel("Count")
    plt.axhline(0, color="black", linewidth=0.5)  # Add vertical line at x=0
    # plt.legend(title='Regulation', labels=['Activation', 'Inhibition'], frameon=False)
    plt.legend(handles=legend_elements, title="Regulation", frameon=False)
    # plt.xticks(rotation=45)
    plt.tight_layout()
    sns.despine()

    # Save the plot
    plt.savefig(bar_pdf)


def main():
    metadata_tsv = sys.argv[1]
    family_bar_pdf = sys.argv[2]
    specificity_bar_pdf = sys.argv[3]

    df = pd.read_csv(metadata_tsv, sep="\t")

    plot_kinase_family_regulation_counts(df, kinase_families, family_bar_pdf)
    
    plot_kinase_specificity_regulation_counts(
        df, kinase_specificities, specificity_bar_pdf
    )


if __name__ == "__main__":
    main()
