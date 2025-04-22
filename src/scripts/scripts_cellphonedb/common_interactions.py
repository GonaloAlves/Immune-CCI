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


def remove_rows():
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
    
    control = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_means_final_merged_uninjured_nona.txt"
    injured_15 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_injured_15_nona.txt"
    injured_60 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_injured_60_nona.txt"

    # Load and simplify
    control_df, injured_15_df, injured_60_df = load_and_simplify(control, injured_15, injured_60)
    
    # Optional: Check shapes
    print("Control shape:7", control_df.shape)
    print("Injured 15 shape:", injured_15_df.shape)
    print("Injured 60 shape:", injured_60_df.shape)

    print("Control shape:", control_df)
    print("Injured 15 shape:", injured_15_df)
    print("Injured 60 shape:", injured_60_df)


