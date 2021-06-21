import pickle
import pandas as pd
import argparse


parser = argparse.ArgumentParser(
    description="Combine dataframes from two antibody chains after counting interactions."
)

parser.add_argument(
    "-path",
    dest="path",
    help="path to pickle files containing interaction counts, structured as /path/to/files/runX, where X is an int"
)

parser.add_argument(
    "-run_number",
    dest="run_number",
    help="the F@H run number"
)

parser.add_argument(
    "-project_number",
    dest="project_number",
    help="the F@H project number"
)


args = parser.parse_args()

def combine_dataframes(data_path, project_number, chain):

    chain_df_list = []
    for i in range(0, 5):

        try:
            with open(f"{data_path}run{run_number}/PROJ{project_number}_count_interactions_angstroms_chain{str(chain)}_results.pkl", "rb") as f:
                df = pickle.load(f)

            chain_df_list.append(df)
        else:
            continue

    # get the attributes
    chain_attrs = chain_df_list[0].attrs

    total_frames = 0
    total_micros = 0
    for df in chain_df_list:
        total_frames += df.attrs["n_frames"]
        total_micros += df.attrs["microseconds"]

    df_out = chain_df_list[0]

    for df in chain_df_list[1:]:

        df_add = df_out.add(df)
        df_out = df_add

    df_out.attrs = chain_attrs
    df_out.attrs["n_frames"] = total_frames
    df_out.attrs["microseconds"] = total_micros

    df_out.to_pickle(f"PROJ{project_number}_count_interactions_angstroms_chain{str(chain)}_results_allruns.pkl")

if __name__ == "__main__":

    # loop over chains 1 and 2
    for i in range(1,3):
        combine_dataframes(
            data_path=argspath,
            run_number=args.run_number,
            project_number=args.project_number,
            chain=i
        )
