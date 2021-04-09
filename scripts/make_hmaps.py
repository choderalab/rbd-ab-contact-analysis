import argparse
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="create heat maps from count data")
parser.add_argument(
    "-data_file",
    dest="data_file",
    type=str,
    help="the pickle file of interaction counts",
)
parser.add_argument(
    "-save_path",
    dest="save_path",
    type=str,
    help="the path to store the generatred plots",
)
parser.add_argument(
    "-ab_name", dest="ab_name", type=str, help="the name of the antibody, e.g. S309"
)
parser.add_argument(
    "-cutoff", dest="cutoff", type=str, help="the cut off used to generate the "
)
parser.add_argument(
    "-custom_residues_file",
    dest="custom_residues_file",
    type=str,
    help="yaml file to use for custom RBD residues order",
)

args = parser.parse_args()


def load_yaml(file):
    import yaml

    with open(file) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    return data


def plot_hmap(
    df,
    ylabel,
    xlabel,
    n_frames,
    save_path,
    percent=True,
    figsize=(12, 12),
    vmax_val=False,
    vmin_val=False,
    plot_name="",
    cmap="rocket",
    cbar_title="",
):

    """
    Plot a heatmap

    Parameters
    ----------

    df : pandas.DataFrame
        A 2D pandas dataframe.
    ylabel : str
        The ylabel title.
    xlabel : str
        The xlabel title.
    n_frames : int
        The total number of frames
    vmax_val : int, float
        The maximum value to use for the colourbar.
    vmin_val : int, float
        The minimum value to use for the colourbar.
    plot_name : str
        The name of the plot, useful for when save=True.
    save : bool
        Set to True to save the plot as a .svg file.
    cmap : str
        The seaborn colourmap to use
    cbar_title : str
        The colourbar title.

    """

    if percent:

        df = (df / n_frames) * 100
        df = df.round(1)

    fig, ax = plt.subplots(figsize=figsize)
    cbar_ax = fig.add_axes([0.92, 0.3, 0.03, 0.4])

    g1 = sns.heatmap(
        data=df,
        annot=True,
        cbar=True,
        cmap=cmap,
        ax=ax,
        vmin=vmin_val,
        vmax=vmax_val,
        cbar_ax=cbar_ax,
        cbar_kws={"label": f"{cbar_title}"},
        fmt="g",
    )

    g1.set_ylabel(ylabel, size=20)
    g1.set_xlabel(xlabel, size=20)
    g1.yaxis.labelpad = 20
    g1.xaxis.labelpad = 20

    # sort the colour bar
    cbar_ax.yaxis.label.set_size(20)
    cbar_ax.yaxis.labelpad = 10

    # rotate the tick labels
    for ax in [g1]:
        tl = ax.get_xticklabels()
        ax.set_xticklabels(tl, rotation=90, size=10)
        tly = ax.get_yticklabels()
        ax.set_yticklabels(tly, rotation=0)

    ax.tick_params(axis="both", which="major", labelsize=12)
    cbar_ax.tick_params(axis="both", which="major", labelsize=12)

    save_name = os.path.join(save_path, plot_name)
    plt.savefig(save_name + ".svg", dpi=300, bbox_inches="tight")


def plot_collapsed_hmap(
    df,
    plot_name,
    xlabel,
    total_ints,
    save_path,
    transpose=True,
    cbar_title="",
    figsize=(12, 12),
    custom_list=False,
    make_percent=False,
    n_frames=None,
    average="frames",
    vmin=0,
    vmax=25,
    annot=False,
):

    """
    Plot a collapsed heatmap

    Parameters
    ----------

    df : pandas.DataFrame
        A 2D pandas dataframe.
    ylabel : str
        The ylabel title.
    xlabel : str
        The xlabel title.
    transpose : bool
        Whether or not to use the transpose of data, this corresponds to the two protein chains
    total_ints : int
        The total number of interactions over all (selected) residues
    n_frames : int
        The total number of frames
    vmax_val : int, float
        The maximum value to use for the colourbar.
    vmin_val : int, float
        The minimum value to use for the colourbar.
    plot_name : str
        The name of the plot, useful for when save=True.
    save : bool
        Set to True to save the plot as a .svg file.
    cmap : str
        The seaborn colourmap to use
    cbar_title : str
        The colourbar title.

    """

    fig, ax = plt.subplots(figsize=figsize)

    if transpose:
        df = df.T

    # sum columns, returns df with one column
    df_sum = df.sum()

    df2 = pd.DataFrame(df_sum)

    # average over all frames
    if average == "frames":
        df2 = (df2 / n_frames) * 100

    # average over interactions
    elif average == "ints":
        # average over user supplied counted interactions
        if total_ints:
            df2 = (df2 / total_ints) * 100
        # average over all counted interactions
        else:
            df2 = (df2 / df2.sum()) * 100

    # transpose to get one row and N columns
    if custom_list:
        # use custom residue list if supplied
        df_plot = df2.T[custom_list]
    else:
        df_plot = df2.T

    # plot the counts
    g1 = sns.heatmap(
        data=df_plot.round(2),
        linewidths=1,
        linecolor="black",
        square=True,
        cmap="Blues",
        ax=ax,
        cbar_kws={"shrink": 0.25, "label": f"{cbar_title}"},
        vmin=vmin,
        vmax=vmax,
        annot=annot,
    )

    cbar = g1.collections[0].colorbar
    cbar.set_label(label=cbar_title, size=12, labelpad=10)

    g1.tick_params(left=False)
    g1.set_yticklabels("", size=12)
    g1.set_xticklabels(g1.get_xticklabels(), rotation=90, size=12)

    ax.set_xlabel(xlabel, size=20)
    ax.set_ylabel("", size=20)
    ax.xaxis.labelpad = 20

    save_name = os.path.join(save_path, plot_name)
    plt.savefig(save_name + ".svg", dpi=300, bbox_inches="tight")

    # save the CSV
    print("saving CSV...")
    df_csv = df2.round(2)
    df_csv.to_csv(save_name + ".csv")


if __name__ == "__main__":

    # gather info
    data_file = args.data_file
    save_path = args.save_path
    ab_name = args.ab_name
    cutoff = args.cutoff
    custom_residues_file = args.custom_residues_file

    with open(data_file, "rb") as data:
        df = pickle.load(data)

    # get system info from loaded df
    ab_chain_key = df.attrs["ab_chain_key"]
    n_micros = int(np.round(df.attrs["microseconds"], 0))
    n_frames = df.attrs["n_frames"]
    mutant = df.attrs["mutant"]
    ab_residues = list(df.index)
    rbd_residues = list(df.columns)

    # specify glycan residue names
    glycan_strings = ("UYB ", "4YB ", "VMB ", "0fA ", "2MA ", "0LB ", "0YB ")

    # separate protein and glycan residues
    ab_protein_residues = [i for i in ab_residues if not i.startswith(glycan_strings)]
    rbd_protein_residues = [i for i in rbd_residues if not i.startswith(glycan_strings)]
    rbd_glycan_residues = [i for i in rbd_residues if i.startswith(glycan_strings)]

    # Take transpose, this sets RBD residues on the y-axis, Ab residues on the x-axis
    df_T = df.T
    ab_rbd_protein_ints = df_T.loc[rbd_protein_residues, ab_protein_residues]
    ab_protein_rbd_glycan_ints = df_T.loc[rbd_glycan_residues, ab_protein_residues]

    # count all interactions that make up close contacts (protein + glycans)
    total_ints = (
        ab_rbd_protein_ints.sum().sum() + ab_protein_rbd_glycan_ints.sum().sum()
    )

    print("making large heat maps...")

    # large protein heat map
    plot_hmap(
        df=ab_rbd_protein_ints.round(1),
        ylabel=f"RBD ({mutant}) interface residues",
        xlabel=f"{ab_name} interface residues",
        n_frames=n_frames,
        percent=True,
        plot_name=f"RBD_{mutant}_{ab_name}_chain{ab_chain_key}_hmap_percentage_{n_micros}_us",
        vmin_val=0,
        vmax_val=100,
        cmap="mako",
        cbar_title=f"Percentage interaction < {cutoff} $\AA$ (%)",
        save_path=save_path,
        figsize=(26, 12),
    )

    # large glycan heat map
    plot_hmap(
        df=ab_protein_rbd_glycan_ints.round(1),
        ylabel=f"RBD ({mutant}) glycan residues",
        xlabel=f"{ab_name} interface residues",
        n_frames=n_frames,
        percent=True,
        plot_name=f"RBDg_{mutant}_{ab_name}_chain{ab_chain_key}_hmap_percentage_{n_micros}_us",
        vmin_val=0,
        vmax_val=100,
        cmap="mako",
        cbar_title=f"Percentage interaction < {cutoff} $\AA$ (%)",
        save_path=save_path,
        figsize=(22, 12),
    )

    print("making collapsed heatmaps...")

    # collapsed RBD glycan interactions
    if custom_residues_file:
        # use a custom order of RBD residues
        custom_residues = load_yaml(custom_residues_file)[mutant]

        plot_collapsed_hmap(
            df=ab_rbd_protein_ints,
            xlabel=f"RBD ({mutant}) interface residue",
            plot_name=f"RBD_{mutant}_interactions_lt{cutoff}_w_{ab_name}_chain{ab_chain_key}_{n_micros}_us_collapsed",
            save_path=save_path,
            cbar_title=f"Interaction with {ab_name} < {cutoff} Å (%)",
            figsize=(18, 12),
            custom_list=custom_residues,
            average="ints",
            total_ints=total_ints,
            annot=True,
        )

    else:
        plot_collapsed_hmap(
            df=ab_rbd_protein_ints,
            xlabel=f"RBD ({mutant}) interface residue",
            plot_name=f"RBD_{mutant}_interactions_lt{cutoff}_w_{ab_name}_chain{ab_chain_key}_{n_micros}_us_collapsed",
            save_path=save_path,
            cbar_title=f"Interaction with {ab_name} < {cutoff} Å (%)",
            figsize=(18, 12),
            average="ints",
            total_ints=total_ints,
            annot=True,
        )

    # collapsed RBD glycan interactions
    plot_collapsed_hmap(
        df=ab_protein_rbd_glycan_ints,
        xlabel=f"RBD ({mutant}) glycan residue",
        plot_name=f"RBDg_{mutant}_interactions_lt{cutoff}_w_{ab_name}_chain{ab_chain_key}_{n_micros}_us_collapsed",
        save_path=save_path,
        cbar_title=f"Interaction with {ab_name} < {cutoff} Å (%)",
        figsize=(18, 12),
        average="ints",
        total_ints=total_ints,
        annot=True,
    )
