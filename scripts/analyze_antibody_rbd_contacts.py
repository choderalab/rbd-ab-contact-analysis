import argparse
import multiprocessing as mp
import os
import pickle
import time
from functools import partial
from multiprocessing import Manager, Pool
from pathlib import Path

import MDAnalysis as mda
import MDAnalysis.transformations as trans
import numpy as np
import pandas as pd
import xarray as xr
from MDAnalysis.analysis import align, distances, rms

parser = argparse.ArgumentParser(
    description="Analyse contacts at the Ab:RBD interface."
)
parser.add_argument(
    "-pdb_file",
    dest="pdb_file",
    type=str,
    help="the equilibrated structure, with canonical residue numbering",
)
parser.add_argument(
    "-traj_path",
    dest="traj_path",
    type=str,
    help="the location of the FAH trajectories i.e. /path/to/traj/PROJXXXXX/RUNX",
)
parser.add_argument(
    "-stride_clones",
    dest="stride_clones",
    type=int,
    help="the stride length. This will be the number of CLONES skipped during analysis",
    default=1,
)
parser.add_argument(
    "-fah_project_code",
    dest="project_code",
    type=str,
    help="the FAH project code, numbers only",
)
parser.add_argument(
    "-split",
    dest="split_value",
    type=int,
    help="the number used to split the list of trajectories into smaller chunks for pooling",
    default=100,
)
parser.add_argument(
    "-n_clones",
    dest="n_clones",
    type=int,
    help="the total number of clones to analyse",
    default=10000,
)
parser.add_argument(
    "-mutant",
    dest="mutant",
    type=str,
    help="the mutant to analyse e.g. WT",
)
parser.add_argument(
    "-clone_file",
    dest="clone_file",
    help="a .txt file containing the CLONES to analyse on each line. Only these CLONE numbers will be analysed if parsed.",
    default=None,
)
parser.add_argument(
    "-stride_frames",
    dest="stride_frames",
    type=int,
    help="the number of frames to stride in a trajectory",
    default=1,
)
parser.add_argument(
    "-cutoff",
    dest="cutoff",
    type=int,
    help="the distance cutoff to use to select interface residues",
    default=10,
)
parser.add_argument(
    "-chain_info_file",
    dest="chain_info_file",
    help="A YAML file containing chain information",
)
parser.add_argument(
    "-run_number",
    dest="run_number",
    default=0,
    help="The F@H run number, default = 0",
)
parser.add_argument(
    "-calphas",
    dest="calphas",
    default=False,
    help="If True, only calpha atoms are used in the calculation",
)


args = parser.parse_args()


def get_chain_info(chain_info_file):

    """
    Returns a dictionary containing chain information for an Ab:RBD (+ glycan) system.

    Parameters
    ----------
    chain_info_file : str
        The path to a yaml file containing chain information in the format

        ``
        antibody_chains:
          chain_1: C
          chain_2: D
        rbd_chain:
          chain_1: A
        glycan_chain:
          chain_1: B
        ``

    Returns
    -------
    loaded_chain_info : dict
        A dictionary containing chain information
    """

    import yaml

    with open(chain_info_file) as f:
        loaded_chain_info = yaml.load(f, Loader=yaml.FullLoader)

    return loaded_chain_info


def get_traj_list(clones, stride_clones=1):

    """
    Returns a list of trajectory file paths.

    Parameters
    ----------
    clones : list
        A list of clone numbers to analyse
    stride_clones : int, optional
        Specifies the number of clones to stride over. Default = 1

    Returns
    -------
    traj_list : list
        A list of trajectory file paths
    """

    traj_list = []

    for clone_number in clones[::stride_clones]:

        xtc_file = f"{traj_path}run{run_number}-clone{clone_number}.xtc"
        xtc_file_path = Path(xtc_file)

        if xtc_file_path.is_file():
            traj_list.append(xtc_file)

    return traj_list


def get_selection(pdb_file, chain_info, cutoff, calphas):

    """
    Returns two dictionaries of results containing metadata and selection strings.

    Parameters
    ----------
    pdb_file : str
        The path to the PDB file to read.
    chain_info : dict
        A dictionary containing chain information, see get_chain_info().
    cutoff : float
        The cutoff distance (in Ångströms) in which to consider a residue in contact with a different protein chain.
    Returns
    -------
    metadata : dict
        A dictionary of metadata
    selections : dict
        A dictionary of reisdue selections
    """

    # set the reference to be the equilibrated structure
    ref = mda.Universe(pdb_file)
    ref.trajectory[0]

    # set up the selections using the reference
    ab_chain_1 = chain_info["antibody_chains"]["chain_1"]
    ab_chain_2 = chain_info["antibody_chains"]["chain_2"]
    rbd_chain = chain_info["rbd_chain"]["chain_1"]
    glycan_chain = chain_info["glycan_chain"]["chain_1"]

    # set selection string based on chain IDs and cutoff
    ab_chain_1_interface_string = f"(segid {ab_chain_1} and around {cutoff} (segid {rbd_chain} or segid {glycan_chain})) and not name H*"
    ab_chain_2_interface_string = f"(segid {ab_chain_2} and around {cutoff} (segid {rbd_chain} or segid {glycan_chain})) and not name H*"
    rbd_interface_string = f"((segid {rbd_chain} and around {cutoff} (segid {ab_chain_1} or segid {ab_chain_2})) or segid {glycan_chain}) and not name H*"

    # get the initial selection (only contains some atoms of each residue)
    ab_chain_1_sel = ref.select_atoms(ab_chain_1_interface_string)
    ab_chain_2_sel = ref.select_atoms(ab_chain_2_interface_string)
    rbd_sel = ref.select_atoms(rbd_interface_string)

    # set string of residue numbers to be used in selection
    # e.g. '10 11 12'
    ab_chain_1_interface_resids = " ".join(map(str, set(ab_chain_1_sel.resids)))
    ab_chain_2_interface_resids = " ".join(map(str, set(ab_chain_2_sel.resids)))
    rbd_interface_resids = " ".join(map(str, set(rbd_sel.resids)))

    # create final selection strings
    if calphas:
        atom_filter = "and name CA"
    else:
        # only count heavy atoms
        atom_filter = "and not name H*"

    ab_chain_1_interface_sel_string = (
        f"segid {ab_chain_1} and (resid {ab_chain_1_interface_resids} {atom_filter})"
    )

    ab_chain_2_interface_sel_string = (
        f"segid {ab_chain_2} and (resid {ab_chain_2_interface_resids} {atom_filter})"
    )

    # include glycans in RBD selection, enforce using only heavy atoms for glycans
    rbd_sel_string = f"(segid {rbd_chain} and resid {rbd_interface_resids} {atom_filter}) or (segid {glycan_chain} and not name H*)"

    # create final selections, these now contain all atoms (with chosen filter) of each residue
    ab_chain_1_atoms = ref.select_atoms(ab_chain_1_interface_sel_string)
    ab_chain_2_atoms = ref.select_atoms(ab_chain_2_interface_sel_string)
    rbd_atoms = ref.select_atoms(rbd_sel_string)

    # make the metadata results
    # NOTE the residues, resnames, and names returned in metadata are used in other analysis scripts
    metadata = {
        "Ab_chain_1_residues": ab_chain_1_atoms.resids,
        "Ab_chain_2_residues": ab_chain_2_atoms.resids,
        "RBD_residues": rbd_atoms.resids,
        "Ab_chain_1_resnames": ab_chain_1_atoms.resnames,
        "Ab_chain_2_resnames": ab_chain_2_atoms.resnames,
        "RBD_resnames": rbd_atoms.resnames,
        "Ab_chain_1_names": ab_chain_1_atoms.atoms.names,
        "Ab_chain_2_names": ab_chain_2_atoms.atoms.names,
        "RBD_names": rbd_atoms.atoms.names,
    }

    selections = {
        "ab_chain_1_interface_sel_string": ab_chain_1_interface_sel_string,
        "ab_chain_2_interface_sel_string": ab_chain_2_interface_sel_string,
        "rbd_sel_string": rbd_sel_string,
    }

    return metadata, selections


def populate_dict(
    traj,
    pdb_file,
    mutant_sel,
    project_code,
    frames_to_stride,
    selection_strings,
    clone_id,
    equil_time=20,
    discard_time=50,
):

    """
    Returns a dictionary of results containing calculated distances and metadata.

    Parameters
    ----------
    traj : str
        The path to the trajectory to be analysed.
    pdb_file : str
        The path to the PDB file to read.
    mutant_sel : str
        The name of the mutant to analyse e.g. 'WT'.
    project_code : int
        The number of the Folding@Home project code for the system e.g. 17317.
    frames_to_stride : int
        The number of frames to stride in each trajectory.
    selection_strings : dict
        A dictionary containing the selection strings to use, see get_selection().
    clone_id : int
        The clone number associated with the supplied trajectory
    equil_time : int
        The equilibration time (in ns) used for the system, default = 20
    discard_time : int
        The amount of time (in ns) to discard before analysing contacts, default = 50

    Returns
    -------
    results_dict : list
        A list of [xr.DataArray(), xr.DataArray()], each array = N x M x n_frames
    """

    results_dict = {}

    print("--> Analysing trajectory: ", traj)

    # set the universe using the current traj
    u = mda.Universe(pdb_file, traj)

    # set empty lists to store xarrays in the current trajectory
    results_array_1_store, results_array_2_store = [], []

    # loop over frames in the current trajectory
    for ts in u.trajectory[::frames_to_stride]:

        time_ns = ts.time / 1000

        # discard first X ns of a trajectory
        if time_ns < discard_time + equil_time:
            continue

        print(f"--> Current frame: {ts.frame}")
        print(f"--> Current time: {time_ns} ns")

        # make interface selections at current time step
        ab_chain_1_sel_ts = u.select_atoms(
            selection_strings["ab_chain_1_interface_sel_string"]
        )
        ab_chain_2_sel_ts = u.select_atoms(
            selection_strings["ab_chain_2_interface_sel_string"]
        )
        rbd_sel_ts = u.select_atoms(selection_strings["rbd_sel_string"])

        # calculate the distances between residues at the Ab:RBD interface
        dist_array_chain_1 = distances.distance_array(
            ab_chain_1_sel_ts.atoms.positions,
            rbd_sel_ts.atoms.positions,
        )

        dist_array_chain_2 = distances.distance_array(
            ab_chain_2_sel_ts.atoms.positions,
            rbd_sel_ts.atoms.positions,
        )

        results_array_abchain_1 = xr.DataArray(
            dist_array_chain_1,
            dims=("ab_chain_1_sel_atoms", "rbd_sel_atoms"),
            coords={
                "ab_chain_1_sel_atoms": list(ab_chain_1_sel_ts.atoms.names),
                "rbd_sel_atoms": list(rbd_sel_ts.atoms.names),
            },
            attrs={
                "chain_num": 1,
                "clone": clone_id,
            },
        )

        results_array_1_store.append(results_array_abchain_1)

        results_array_abchain_2 = xr.DataArray(
            dist_array_chain_2,
            dims=("ab_chain_2_sel_atoms", "rbd_sel_atoms"),
            coords={
                "ab_chain_2_sel_atoms": list(ab_chain_2_sel_ts.atoms.names),
                "rbd_sel_atoms": list(rbd_sel_ts.atoms.names),
            },
            attrs={
                "chain_num": 2,
                "clone": clone_id,
            },
        )

        results_array_2_store.append(results_array_abchain_2)

    # concatenate xarrays in results lists
    results_array_abchain_1_concat = xr.concat(results_array_1_store, "frames")
    results_array_abchain_2_concat = xr.concat(results_array_2_store, "frames")

    return [results_array_abchain_1_concat, results_array_abchain_2_concat]


def write_result(
    result, i, project_code, mutant, n_clones, stride_clones, stride_frames
):

    """
    Writes a result packet from Pool.

    Parameters
    ----------
    result : list
        A returned list from populate_dict() of [xr.DataArray, xr.DataArray]
    i : int
        A value to label the returned result
    project_code : str
        The FAH project code
    mutant : str
        The system, typically this is WT (wild-type)
    n_clones : int
        The number of clones
    stride_clones : int
        The number of clones strided / skipped
    stride_frames : int
        The number of frames in a trajectory strided / skipped
    """

    for r in result:

        pickle_name = f"PROJ{project_code}_{mutant}_Abchain{r.chain_num}_distances_clone{r.clone}.pkl"

        with open(pickle_name, "wb") as output:
            pickle.dump(r, output)


if __name__ == "__main__":

    # Define everything
    pdb_file = args.pdb_file
    mutant = args.mutant
    project_code = args.project_code
    traj_path = args.traj_path
    n_clones = args.n_clones
    stride_clones = args.stride_clones
    clone_file = args.clone_file
    stride_frames = args.stride_frames
    cutoff = args.cutoff
    chain_info_file = args.chain_info_file
    run_number = args.run_number
    calphas = args.calphas

    # Print for a sanity check
    print(f"--> Using PDB file: {pdb_file}")
    print(f"--> Analysing {mutant} system")
    print(f"--> Using project code: {project_code}")
    print(f"--> Analysing run number: {run_number}")

    print("\n--> Processing chain information")

    chain_info = get_chain_info(chain_info_file)

    if calphas:
        print("\n-->Using C-alpha atoms for distance calculation")

    print(f"--> RBD chain: {chain_info['rbd_chain']}")
    print(f"--> Glycan chain: {chain_info['glycan_chain']}")
    print(f"--> Antibody chains: {chain_info['antibody_chains']}")

    # read in the trajs as a list
    print("\n--> Creating list of trajectories")
    print(f"--> Number of clones to stride: {stride_clones}")
    print(f"--> Number of frames to stride in a trajectory: {stride_frames}")
    print(f"--> Distance cutoff: {cutoff} Å")

    if clone_file is not None:  # using a input file of specific clones

        with open(clone_file) as f:
            temp_file = f.readlines()
            clones_to_analyse = [line.rstrip("\n") for line in temp_file]

        traj_list = get_traj_list(
            clones=clones_to_analyse,
            stride_clones=stride_clones,
        )

    else:  # running without an input file of clones

        clones_to_analyse = [x for x in range(n_clones)]
        traj_list = get_traj_list(
            clones=clones_to_analyse,
            stride_clones=stride_clones,
        )

    print(f"--> Number of trajectory files to analyse: {len(traj_list)}")

    print("\n--> Staring pool...")

    manager = Manager()
    shared_list = manager.list()

    # use all available cores, otherwise specify the number you want as an argument
    pool = Pool(mp.cpu_count())

    # get metadata and selection strings
    metadata_dict, selection_strings = get_selection(
        pdb_file=pdb_file, chain_info=chain_info, cutoff=cutoff, calphas=calphas
    )

    # save metadata file
    metadata_pickle_name = f"PROJ{project_code}_{mutant}_selection_metadata.pkl"
    with open(metadata_pickle_name, "wb") as output:
        pickle.dump(metadata_dict, output)

    ts = time.time()

    # loop over trajectories
    for traj in traj_list:

        # assumes path/to/file/run0-clone0.xtc
        # TODO clean this up?
        clone_id = int(
            traj.split("/")[-1].split(".")[0].split("-")[1].split("clone")[1]
        )

        # use partial to specify associated clone_id with saved pkl
        write_result_i = partial(
            write_result,
            i=clone_id,
            project_code=project_code,
            mutant=mutant,
            n_clones=n_clones,
            stride_clones=stride_clones,
            stride_frames=stride_frames,
        )

        pool.apply_async(
            populate_dict,
            args=(
                traj,
                pdb_file,
                mutant,
                project_code,
                stride_frames,
                selection_strings,
            ),
            callback=write_result_i,
        )

    pool.close()
    pool.join()
    print("Time in parallel:", time.time() - ts)
    print("Finished!")
