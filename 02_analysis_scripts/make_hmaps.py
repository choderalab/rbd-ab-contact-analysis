import argparse
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser(description="create heat maps from count data")
parser.add_argument(
    "-data",
    dest="data",
    nargs="+",
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
parser.add_argument(
    "-percent",
    dest="percent",
    default=True,
    type=bool,
    help="set to False to plot raw counts"
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
    cbar_title="",
    figsize=(12, 12),
    custom_list=False,
    n_frames=None,
    average="ints",
    vmin=0,
    vmax=25,
    annot=False,
    percent=True,
):

    """
    Plot a collapsed heatmap

    Parameters
    ----------

    df : pandas.DataFrame
        A 2D pandas dataframe.
    plot_name : str
        the name of the plot.
    xlabel : str
        The xlabel title.
    total_ints : int
        The total number of interactions over all (selected) residues.
    save_path : str
        The path to save the plot.
    cbar_title : str
        The colourbar title.
    figsize : (int, int)
        The figure size.
    custom_list : .yaml file
        A yaml file containing the residues to plot.
    n_frames : int
        The total number of frames
    average : str
        Average interactions over total interactions (default) or via custom interactions ("ints")
    vmax_val : int, float
        The maximum value to use for the colourbar.
    vmin_val : int, float
        The minimum value to use for the colourbar.
    annot : bool
        Annotate the heatmap plots with values.
    percent : bool
        Plot with percentages



    """

    fig, ax = plt.subplots(figsize=figsize)

    # use custom residue list if supplied
    if custom_list:
        # transpose to get one row and N columns
        df = df.T[custom_list]
    else:
        df = df.T

    # sum columns, returns df with one column
    df_sum = df.sum()

    df2 = pd.DataFrame(df_sum)

    # average over all frames
    if average == "frames":
        print("Averaging counts over frames")
        df_avg = (df2 / n_frames)

    # average over interactions
    elif average == "ints":
        print("Averaging over interactions")
        # average over user supplied counted interactions
        if total_ints:
            df_avg = (df2 / total_ints)
        # average over all counted interactions
        else:
            df_avg = (df2 / df2.sum())
    
    else:
        print("---> INFO: not averaging counts")

    if percent:
        print("Creating percentages...")
        df_plot = df_avg.T * 100
    else:
        df_plot = df_avg.T

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
    df_plot.to_csv(save_name + ".csv")


if __name__ == "__main__":

    # gather info
    data = args.data
    save_path = args.save_path
    ab_name = args.ab_name
    cutoff = args.cutoff
    custom_residues_file = args.custom_residues_file
    percent = args.percent

    # get system info from first loaded df, assume the same for both files if provided
    with open(data[0], "rb") as d:
            df_info = pickle.load(d)

    n_micros = int(np.round(df_info.attrs["microseconds"], 0))
    n_frames = df_info.attrs["n_frames"]
    mutant = df_info.attrs["mutant"]

    # load the data, either one or two chains
    if len(data) > 1:
        with open(data[0], "rb") as data_0:
            df0 = pickle.load(data_0)
        with open(data[1], "rb") as data_1:
            df1 = pickle.load(data_1)

        # combine dataframes (i.e. combine both chains)
        df = pd.concat([df0, df1])

        ab_chain_key = "HL" # heavy and light Ab chains

    else:

        with open(data[0], "rb") as data:
            df = pickle.load(data)

        ab_chain_key = df.attrs["ab_chain_key"]

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

    print(df_T)
    print(ab_rbd_protein_ints)
    print(ab_protein_rbd_glycan_ints)

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
        percent=percent,
        plot_name=f"RBD_{mutant}_{ab_name}_chain{ab_chain_key}_hmap_percentage_{n_micros}_us",
        vmin_val=0,
        vmax_val=100,
        cmap="mako",
        cbar_title=f"Interaction (%)",
        save_path=save_path,
        figsize=(26, 12),
    )

    # large glycan heat map
    plot_hmap(
        df=ab_protein_rbd_glycan_ints.round(1),
        ylabel=f"RBD ({mutant}) glycan residues",
        xlabel=f"{ab_name} interface residues",
        n_frames=n_frames,
        percent=percent,
        plot_name=f"RBDg_{mutant}_{ab_name}_chain{ab_chain_key}_hmap_percentage_{n_micros}_us",
        vmin_val=0,
        vmax_val=100,
        cmap="mako",
        cbar_title=f"Interaction (%)",
        save_path=save_path,
        figsize=(22, 12),
    )

    print("making collapsed heatmaps...")

    # collapsed RBD glycan interactions
    if custom_residues_file:
        # use a custom order of RBD residues
        custom_residues = load_yaml(custom_residues_file)[mutant]

        ab_rbd_protein_ints_custom = df_T.loc[custom_residues, ab_protein_residues]
        print(ab_rbd_protein_ints_custom)

        # count all interactions that make up close contacts (protein + glycans) of custom protein subset
        total_ints_custom = (
            ab_rbd_protein_ints_custom.sum().sum() + ab_protein_rbd_glycan_ints.sum().sum()
        )

        plot_collapsed_hmap(
            df=ab_rbd_protein_ints,
            xlabel=f"RBD ({mutant}) interface residue",
            plot_name=f"RBD_{mutant}_interactions_w_{ab_name}_chain{ab_chain_key}_{n_micros}_us_collapsed",
            save_path=save_path,
            cbar_title=f"Interaction with {ab_name} (%)",
            figsize=(18, 12),
            custom_list=custom_residues,
            average="ints",
            total_ints=total_ints_custom,
            annot=True,
        )

    else:
        plot_collapsed_hmap(
            df=ab_rbd_protein_ints,
            xlabel=f"RBD ({mutant}) interface residue",
            plot_name=f"RBD_{mutant}_interactions_w_{ab_name}_chain{ab_chain_key}_{n_micros}_us_collapsed",
            save_path=save_path,
            cbar_title=f"Interaction with {ab_name} (%)",
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
        cbar_title=f"Interaction with {ab_name} (%)",
        figsize=(18, 12),
        average="ints",
        total_ints=total_ints,
        annot=True,
    )
