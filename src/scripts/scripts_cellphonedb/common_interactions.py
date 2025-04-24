# This script will take the significant interations from cellphonedb output and its gonna be removing the interations from the uninjered/sham file 
# which will only be left with the injured new ones 
#
# Goncalo Alves Msc Thesis 2025

import pandas as pd
import scanpy as sc
import numpy as np

def load_and_simplify(control_path, injured_15_path, injured_60_path):
    """
    Loads the p-values tables from the three conditions and keeps only 
    the 'id_cp_interaction' and the cluster interaction columns (those with '|').
    Returns three simplified DataFrames.
    """

    control_df = simplify_table(control_path)
    injured_15_df = simplify_table(injured_15_path)
    injured_60_df = simplify_table(injured_60_path)

    return control_df, injured_15_df, injured_60_df


def simplify_table(file_path):
      
    df = pd.read_csv(file_path, sep='\t')
    # Keep only 'id_cp_interaction' and columns with '|' (which identify cluster interactions)
    columns_to_keep = ['id_cp_interaction'] + [col for col in df.columns if '|' in col]
    simplified_df = df[columns_to_keep]

    return simplified_df

def export_to_excel(df, output_path):
    """
    Exports the given DataFrame to an Excel file.
    
    Parameters:
    - df: pd.DataFrame
    - output_path: str, full path to save the Excel file
    """
    try:
        df.to_excel(output_path, index=False)
        print(f"DataFrame exported successfully to {output_path}")
    except Exception as e:
        print(f"Failed to export DataFrame: {e}")


def df_to_significant_dict(df, threshold=0.05):
    """
    Converts a simplified CellPhoneDB p-values dataframe into a dictionary.
    Each key is an 'id_cp_interaction', and its value is a list of cluster interaction
    columns where the p-value is below the threshold (default: 0.05).
    """
    sig_dict = {}
    cluster_columns = [col for col in df.columns if col != 'id_cp_interaction']

    for _, row in df.iterrows():
        interaction_id = row['id_cp_interaction']
        significant_clusters = [
            col for col in cluster_columns if row[col] < threshold
        ]
        if significant_clusters:
            sig_dict[interaction_id] = significant_clusters

    return sig_dict


def filter_injured_by_control(control_dict, injured_dict, verbose=True):
    """
    Removes cluster interactions from the injured_dict that are also present 
    in the control_dict for the same interaction ID.
    Returns a new filtered dictionary.

    If verbose=True, it prints step-by-step details for debugging.
    """
    filtered_dict = {}

    for interaction_id, injured_clusters in injured_dict.items():
        control_clusters = control_dict.get(interaction_id, [])

        if verbose:
            print(f"\n--- Checking interaction ID: {interaction_id} ---")
            print(f"Injured clusters: {injured_clusters}")
            if control_clusters:
                print(f"Control clusters: {control_clusters}")
            else:
                print("No matching interaction in control (new interaction).")

        # Compare and filter
        unique_clusters = [cl for cl in injured_clusters if cl not in control_clusters]

        if verbose:
            removed = set(injured_clusters) - set(unique_clusters)
            if removed:
                print(f"Removed clusters (shared with control): {list(removed)}")
            if unique_clusters:
                print(f"Remaining (injury-specific) clusters: {unique_clusters}")
            else:
                print("No injury-specific clusters remain after filtering.")

        if unique_clusters:
            filtered_dict[interaction_id] = unique_clusters

    return filtered_dict


def transform_dict():
  """
    This function will remove the unwanted rows from each df 
    (we only want the interations that are only significant and unique in the injured condictions)
    
    Ent a minha ideia é tendo em conta o facto que o txt dos pvalues tem todas as interacoes que tem e n tem pvalues decentes
    a ideia seria primeiro remover todas as rows que n tem pelo menos um pvalues menor que 0.05.
    dps podiamos isolar as interacoes itself, ou seja temos uma row com um id, e dps remover todas as colunas que tem o pvalue superior a 0.05.
    dps o proximo step seria criar um df com cada id de interacao por row e as colunas seriam os conjuntos de clusters que sao detetados as interacoes.
   e dps ias ver os ids e clusters que se verificam no uninjured e ver quais os que se repetem os outros e se sim retirar
"""

# Main execution block
if __name__ == "__main__":
    # Load data
    
    control = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_uninjured_nona.txt"
    injured_15 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_injured_15_nona.txt"
    injured_60 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_injured_60_nona.txt"

    # Load and simplify
    control_df, injured_15_df, injured_60_df = load_and_simplify(control, injured_15, injured_60)
    
    # # Check shapes
    # print("Control shape:7", control_df.shape)
    # print("Injured 15 shape:", injured_15_df.shape)
    # print("Injured 60 shape:", injured_60_df.shape)

    # print("Control shape:", control_df)
    # print("Injured 15 shape:", injured_15_df)
    # print("Injured 60 shape:", injured_60_df)

    # Convert to dictionaries of significant interactions
    control_dict = df_to_significant_dict(control_df)
    injured_15_dict = df_to_significant_dict(injured_15_df)
    injured_60_dict = df_to_significant_dict(injured_60_df)

    # Filter injury-specific interactions
    filtered_15_dict = filter_injured_by_control(control_dict, injured_15_dict, verbose=False)
    filtered_60_dict = filter_injured_by_control(control_dict, injured_60_dict, verbose=False)

    from pprint import pprint
    print("\nExample from Injured 15:")
    pprint(list(injured_15_dict.items())[:3])

    # export_to_excel(control_df, "control_simplified.xlsx")
    # export_to_excel(injured_15_df, "injured_15_simplified.xlsx")
    # export_to_excel(injured_60_df, "injured_60_simplified.xlsx")


