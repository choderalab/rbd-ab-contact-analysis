import argparse
import glob
import multiprocessing as mp
import os
import pickle
import time
from functools import partial
from multiprocessing import Manager, Pool
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

parser = argparse.ArgumentParser(
    description="count interactions at the Ab:RBD interface"
)
parser.add_argument(
    "-path_name",
    dest="path_name",
    type=str,
    help="the path to pickle files containing distance arrays",
)
parser.add_argument(
    "-metadata_file",
    dest="metadata_file",
    type=str,
    help="the name of the pickle file to read for metadata",
)
parser.add_argument(
    "-file_slice_val",
    dest="file_slice_val",
    type=int,
    default=1,
    help="the number of files to skip",
)
parser.add_argument(
    "-frame_slice_val",
    dest="frame_slice_val",
    type=int,
    default=1,
    help="the number of frames to skip in a trajectory",
)
parser.add_argument(
    "-save_collected_data",
    dest="save_collected_data",
    type=bool,
    default=True,
    help="specify whether to save the data",
)
parser.add_argument(
    "-distance_cutoff",
    dest="cutoff",
    # type=float,
    # default=3.5,
    help="the cutoff distance to use to count an interaction",
)
parser.add_argument(
    "-ab_chain_key",
    dest="ab_chain_key",
    type=int,
    default=1,
    help="the antibody chain to analyse, e.g. for S309 this can either be 1 or 2",
)
parser.add_argument(
    "-project_number",
    dest="project_number",
    type=int,
    default=1,
    help="the F@H project number",
)
parser.add_argument(
    "-mutant",
    dest="mutant",
    type=str,
    help="the mutant, e.g. WT, P337A etc",
)
parser.add_argument(
    "-n_runs",
    dest="n_runs",
    type=int,
    default=1,
    help="the number of runs to analyse, default=3.5",
)

args = parser.parse_args()


def load_pickle(
    file=None,
    metadata=False,
    file_slice_val=1,
    frame_slice_val=1,
    save=False,
):

    """
    Loads data from pickle file(s) and returns a dictionary containing results
    Parameters
    ----------
    path_name : str
        The path to the .pkl file(s)
    metadata : bool
        Request that only metadata is returned
    file_slice_val : int
        The number of files to skip over, default = 1
    frame_slice_val : int
        The number of frames to skip over, default = 1
    Returns
    -------
    data : dict
        Dictionary of {"data": dict} or {"metadata": dict}
    """

    if metadata:
        data = {"metadata": None}

        with open(file, "rb") as raw_data:
            data["metadata"] = pickle.load(raw_data)

    else:
        data = {"data": {}}

        with open(file, "rb") as raw_data:

            loaded_data = pickle.load(raw_data)

        frames_to_read = loaded_data[::frame_slice_val, :, :]

        data["data"] = frames_to_read

    return data


def make_results_array(
    ab_resids,
    ab_resnames,
    rbd_resids,
    rbd_resnames,
):

    """
    Creates an empty xarray to be populated
    Parameters
    ----------
    ab_resids : np.ndarray
        1D np array of antibody residue IDs
    ab_resname_list : list
        1D np array of antibody residue names
    rbd_resid_list : list
        1D np array of RBD residue IDs
    rbd_resname_list : list
        1D np array of RBD residue names
    Returns
    -------
    results_array : xarray
        An empty xarray to be populated
    """

    # combine single letter AA code with resid
    aa_dict = {
        "ACE": "ACE ",  # C-terminal cap
        "NME": "NME ",  # N-terminal cap
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "NLN": "N",  # glycosylated ASN
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        "UYB": "UYB ",  # glycan
        "4YB": "4YB ",  # glycan
        "2MA": "2MA ",  # glycan
        "0LB": "0LB ",  # glycan
        "VMB": "VMB ",  # glycan
        "0fA": "0fA ",  # glycan
        "0YB": "0YB ",  # glycan
    }

    # get a dict of residues as keys in the form, for example, {D45 : 0, A46 : 0, R47 : 0}
    ab_res_resids_d = {}
    rbd_res_resids_d = {}

    for i, j in zip(ab_resnames, ab_resids):

        entry = aa_dict[i] + str(j)
        ab_res_resids_d[entry] = 0  # placeholder

    for i, j in zip(rbd_resnames, rbd_resids):

        entry = aa_dict[i] + str(j)
        rbd_res_resids_d[entry] = 0  # placeholder

    # create an xarray of zeros to store counting data
    results_array = xr.DataArray(
        np.zeros((len(ab_res_resids_d), len(rbd_res_resids_d))),
        dims=("ab_residues", "rbd_residues"),
        coords={
            "ab_residues": list(ab_res_resids_d.keys()),
            "rbd_residues": list(rbd_res_resids_d.keys()),
        },
    )

    return results_array


def create_alt_resids(resids, glycan_start):

    """
    Creates alternative residue numbering (000) for glycosylated protein chains
    Parameters
    ----------
    resids : np.ndarray
        1D np array of antibody residue IDs
    glycan_start : int
        The index in a residue ID list where the glycan residue begins
    Returns
    -------
    resids_alt : np.array
        A 1D np.array containing new residue numberings
    """

    # add "000" to resid list to prevent duplication with other residues
    resids_alt = []
    for i, res in enumerate(resids):
        if i < glycan_start:  # index where glycans start
            resids_alt.append(res)
        else:
            resids_alt.append(int(str(res) + "000"))

    return np.array(resids_alt)


def get_result(result):
    """
    Process results packet from Pool process
    """
    global_list.append(result)


def count_percent_int(
    results_array,
    clone_data,
    ab_resids,
    rbd_resids,
    res_resid_ab_list,
    res_resid_rbd_list,
    user_cutoff,
):

    """
    Count interactions of residue pairs based on minimum distances
    Parameters
    ----------
    results_array : xarray
        An empty array generated from make_results_array
    clone_data : xarray
        An xarray loaded from a pickle file containing distance arrays for a particular trajectory
    ab_resids : np array
        A 1D numpy array containing antibody residue IDs (for each atom)
    rbd_resids : np array
        A 1D numpy array containing rbd residue IDs (for each atom), corrected for any glycans present
    res_resid_ab_list : list
        A list of antibody residues in the form X1, Y2, Z3 etc
    res_resid_rbd_list : list
        A list of RBD residues in the form X1, Y2, Z3 etc
    user_cutoff : float
        Optional, the cutoff at which to count an interaction defined by the user.
    Returns
    -------
    results_array : xarray
        A populated 2D array with counts
    """

    hydrophobic_residues = ["ALA", "VAL", "ILE", "LEU", "MET","PHE","TYR","TRP"]

    for frame in clone_data:
        # loop over the residues in the first residue list
        for resi in res_resid_ab_list:
            # create the index based on the residue number (e.g. 45 for D45)
            try:
                resn = int(resi[1:])
            except:
                # handle non-standard residues like ACE, NME
                if resi[:3] == "ACE" or resi[:3] == "NME":
                    resn = int(resi[3:])

            # locate ab resid indices
            maskn = ab_resids == resn

            for resj in res_resid_rbd_list:

                try:
                    resm = int(resj[1:])
                except:
                    if resj[:3] == "ACE" or resj[:3] == "NME":
                        resm = int(resj[3:])
                    else:
                        # handle modified glycan resids
                        resm = int(str(f"{resj[3:]}000"))

                # locate rbd resid indices
                maskm = rbd_resids == resm

                # create the mask for the 2D subset
                mm, nn = np.meshgrid(maskm, maskn)  # x, y
                mask = np.logical_and(mm, nn)

                # create subset using a mask of the 2D distance array
                # get the minimum of this 2D subset (i.e. of all distances between residue X and Y atoms)
                # this can be done since the N x M array at the current frame has the same atom order
                # as in ab_resids / rbd_resids

                if user_cutoff:
                    cutoff = user_cutoff

                else:
                    if resi or resj in hydrophobic_residues:
                        cutoff = 4.5
                    else:
                        cutoff = 3.5

                min_val = np.min(frame.values[mask]) <= cutoff

                # locate corresponding residues in results_array (xarray)
                if min_val:
                    results_array.loc[dict(ab_residues=resi, rbd_residues=resj)] += 1

    return results_array


if __name__ == "__main__":

    # Define everything
    path_name = args.path_name
    file_slice_val = args.file_slice_val
    frame_slice_val = args.frame_slice_val
    metadata_file = args.metadata_file
    save_collected_data = args.save_collected_data
    cutoff = args.cutoff
    ab_chain_key = args.ab_chain_key
    project_number = args.project_number
    mutant = args.mutant
    n_runs = args.n_runs

    print("Loading metadata...")
    metadata = load_pickle(
        file=metadata_file,
        metadata=True,
    )

    print(f"Cutoff: {cutoff}")

    print("Finding protein residues...")

    # get resids and resnames, this contains all the atoms of each residue
    # (i.e. multiple resid occurances)
    rbd_resids = metadata["metadata"]["RBD_residues"]
    rbd_resnames = metadata["metadata"]["RBD_resnames"]

    if ab_chain_key == 1:
        ab_resids = metadata["metadata"]["Ab_chain_1_residues"]
        ab_resnames = metadata["metadata"]["Ab_chain_1_resnames"]

    elif ab_chain_key == 2:
        ab_resids = metadata["metadata"]["Ab_chain_2_residues"]
        ab_resnames = metadata["metadata"]["Ab_chain_2_resnames"]

    else:
        print("You need to specify a chain key")
        os.sys.exit()

    print("Finding RBD glycan residues...")
    # Get the location of where the glycan residues begin
    # this is a hack since the RBD chain starts at resid 1
    for i, r in enumerate(rbd_resids):
        if i > 20 and r == 1:
            break

    rbd_glycan_start = i

    print("Creating an empty array to populate with interactions...")

    results_array = make_results_array(
        ab_resids=ab_resids,
        ab_resnames=ab_resnames,
        rbd_resids=rbd_resids,
        rbd_resnames=rbd_resnames,
    )

    print("Sorting RBD glycan residues...")
    # Create alternative RBD resids to handle glycans
    rbd_resids_alt = create_alt_resids(rbd_resids, rbd_glycan_start)

    # get the names of the residues in the form X1, Y2, Z3, etc from the empty results array
    res_resid_ab_list = list(results_array.coords["ab_residues"].values)
    res_resid_rbd_list = list(results_array.coords["rbd_residues"].values)

    print("Setting up pool...")
    manager = Manager()
    global_list = manager.list()

    # use all available cores, otherwise specify the number you want as an argument
    pool = Pool(mp.cpu_count())

    ts = time.time()
    print("Starting loop...")

    # get list of clones to analyse
    # TODO remove hard coded clone number?
    clones_to_analyse = []
    for clone_number in range(0, 5000, file_slice_val):

        for run_number in range(n_runs):

            file = os.path.join(
                path_name,
                f"PROJ{project_number}_{mutant}_Abchain{ab_chain_key}_distances_run_{str(run_number)}_clone{clone_number}.pkl",
            )

            if os.path.isfile(file):
                clones_to_analyse.append(file)

    # Loop over each clone
    n_frames = 0

    for clone_file in clones_to_analyse:

        loaded_data = load_pickle(
            file=clone_file,
            frame_slice_val=frame_slice_val,
        )

        n_frames += len(loaded_data["data"])

        pool.apply_async(
            count_percent_int,
            args=(
                results_array,  # empty results array
                loaded_data["data"],  # distance arrays for clone
                ab_resids,  # antibody residue IDs (for each atom)
                rbd_resids_alt,  # rbd residue IDs (for each atom) with glycans "corrected"
                res_resid_ab_list,  # antibody residues [X1, Y2, Z3]
                res_resid_rbd_list,  # rbd residues [X1, Y2, Z3]
                cutoff,
            ),
            callback=get_result,
        )

    pool.close()
    pool.join()
    print("Time in parallel:", time.time() - ts)
    print(f"{len(global_list)} clones analysed")
    print(f"{n_frames} frames / {n_frames / 1000} microseconds (aggregated) analysed")

    print("Concatenating results...")
    # concatenate all count arrays across all clones
    concat_array = xr.concat(global_list, "frames")
    concat_array_sum = concat_array.sum("frames")

    print("Creating final dataframe...")
    df = concat_array_sum.to_pandas()
    df.attrs["n_frames"] = n_frames
    df.attrs["microseconds"] = n_frames / 1000 # each frame = 1 ns
    df.attrs["ab_resids"] = ab_resids
    df.attrs["ab_resnames"] = ab_resnames
    df.attrs["rbd_resids"] = rbd_resids
    df.attrs["rbd_resnames"] = rbd_resnames
    df.attrs["project_number"] = project_number
    df.attrs["mutant"] = mutant
    df.attrs["ab_chain_key"] = ab_chain_key

    if not cutoff:
        pickle_name = f"PROJ{project_number}_count_interactions_angstroms_chain{ab_chain_key}_results.pkl"

    else:
        pickle_name = f"PROJ{project_number}_count_interactions_below_{cutoff}_angstroms_chain{ab_chain_key}_results.pkl"
    
    df.to_pickle(pickle_name)

    print(f"Results saved as: {pickle_name}")
