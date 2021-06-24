import pickle
import pandas as pd

def combine_dataframes(chain):


    chain_df_list = []
    for i in range(0, 5):

        with open(f"../run{i}/PROJ17341_count_interactions_angstroms_chain{str(chain)}_results.pkl", "rb") as f:
            df = pickle.load(f)

        chain_df_list.append(df)

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

    df_out.to_pickle(f"PROJ17341_count_interactions_angstroms_chain{str(chain)}_results_allruns.pkl")

if __name__ == "__main__":

    combine_dataframes(1)
    combine_dataframes(2)
