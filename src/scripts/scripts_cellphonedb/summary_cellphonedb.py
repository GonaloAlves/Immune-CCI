# Perform CellPhoneDB analysis - LeonorSaude10x
#
# 
#
# Daniel Ribeiro, Gonçalo Alves 2025
import os
import pandas as pd

cpdb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database"
cellphonedb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list"
cellphonedb_dir_out = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary"


datset_names = {'Neu': 'Neuron',
                'MeV': 'Meningeal_Vascular',
                'Imm': 'Immune'
                }

def simple_uniprot_to_gene_name(df_simple: pd.DataFrame,      # gene_input.csv dataframe
                                uniprot: str) -> str:
    mask = df_simple.loc[:, "uniprot"] == uniprot
    subset = df_simple[mask]
    # Get fisrt row
    row = subset.iloc[0, :]
    
    return row["gene_name"]


def complex_uniprot_to_gene_name(df_complex: pd.DataFrame,      # complex_input.csv dataframe
                                 df_simple: pd.DataFrame,       # gene_input.csv dataframe
                                 complex_name: str) -> list[str]:
    result = []
    mask = df_complex.loc[:, "complex_name"] == complex_name
    subset = df_complex[mask]
    row = subset.iloc[0, :]
    # Get all possilbe Uniprot ids
    row = row.loc[["uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4", "uniprot_5"]]
    uniprot = []
    uniprot.extend(row[pd.notna(row)].values.tolist())
    for uni in uniprot:
        result.append(simple_uniprot_to_gene_name(df_simple, uni))

    return result


def genes_from_partner(df_complex: pd.DataFrame,      # complex_input.csv dataframe
                       df_simple: pd.DataFrame,       # gene_input.csv dataframe
                       partner: str):
    result = []
    if partner.find("simple:") != -1:
        result.append(simple_uniprot_to_gene_name(df_simple, partner.split(':')[1]))
    elif partner.find("complex:") != -1:
        result.extend(complex_uniprot_to_gene_name(df_complex, df_simple, partner.split(':')[1]))
    
    
    return result

def get_genes_from_id_cp_interaction(interaction_id: str, df_stat: pd.DataFrame, df_complex: pd.DataFrame, df_simple: pd.DataFrame) -> list[str]:
    
    # Find row in df_stat for this interaction
    row = df_stat[df_stat['id_cp_interaction'] == interaction_id]
    if row.empty:
        return []

    partner_a = row['partner_a'].values[0]
    partner_b = row['partner_b'].values[0]
    
    genes = []
    # Handle partner_a
    genes += genes_from_partner(df_complex, df_simple, partner_a)
    # Handle partner_b
    genes += genes_from_partner(df_complex, df_simple, partner_b)

    return sorted(set(genes))  # optional: remove duplicates


def enrich_simplified_excel_with_genes(simplified_excel_path: str, output_excel_path: str, cpdb_gene_input_path: str, cpdb_complex_input_path: str, cpdb_stat_path: str):

    # Load CellPhoneDB inputs
    df_simple = pd.read_csv(cpdb_gene_input_path, sep=',')
    df_complex = pd.read_csv(cpdb_complex_input_path, sep=',')
    df_stat = pd.read_csv(cpdb_stat_path, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})

    # Load all sheets
    xl = pd.read_excel(simplified_excel_path, sheet_name=None)

    updated_sheets = {}
    for sheet_name, df in xl.items():
        print(f"Processing {sheet_name}...")
        genes_list = []

        for interaction_id in df["id_cp_interaction"]:
            genes = get_genes_from_id_cp_interaction(interaction_id, df_stat, df_complex, df_simple)
            genes_list.append(", ".join(genes))  # Join as string for Excel

        df["genes"] = genes_list
        updated_sheets[sheet_name] = df

    # Write to new Excel
    with pd.ExcelWriter(output_excel_path, engine='xlsxwriter') as writer:
        for sheet_name, df in updated_sheets.items():
            df.to_excel(writer, index=False, sheet_name=sheet_name)

    print("✅ Done writing enriched Excel.")



def alphabetic_order_int(filepath: str, output_path: str):
    """
    Combine all sheets into one DataFrame, keep only relevant columns,
    and sort interactions alphabetically. Includes cluster_pair in the first column.

    Args:
        filepath (str): Path to input Excel file with multiple sheets.
        output_path (str): Path to save the output Excel file.
    """
    # Load all sheets
    xls = pd.read_excel(filepath, sheet_name=None)
    all_data = []

    for sheet_name, df in xls.items():
        if df.empty:
            continue

        df = df.copy()

        # Build cluster_pair column
        df["cluster_pair"] = df.apply(
            lambda row: f"{sheet_name}|{row['partner']}" if row["direction"] == "sent"
            else f"{row['partner']}|{sheet_name}",
            axis=1
        )

        # Drop unneeded columns
        df.drop(columns=["direction", "partner"], inplace=True)

        # Move cluster_pair to first column
        cols = ["cluster_pair"] + [col for col in df.columns if col != "cluster_pair"]
        df = df[cols]

        all_data.append(df)

    # Merge all sheets and sort
    merged_df = pd.concat(all_data, ignore_index=True)
    merged_df = merged_df.sort_values(by="interacting_pair", ascending=True)

    # Save result
    merged_df.to_excel(output_path, index=False)
    print(f"Processed file saved to: {output_path}")


def alphabetic_order_cluster(filepath: str, output_path: str):
    """
    Combine all sheets into one DataFrame, keep only relevant columns,
    and sort interactions by cluster_pair and then by interacting_pair alphabetically.

    Args:
        filepath (str): Path to input Excel file with multiple sheets.
        output_path (str): Path to save the output Excel file.
    """
    # Load all sheets
    xls = pd.read_excel(filepath, sheet_name=None)
    all_data = []

    for sheet_name, df in xls.items():
        if df.empty:
            continue

        df = df.copy()

        # Build cluster_pair column
        df["cluster_pair"] = df.apply(
            lambda row: f"{sheet_name}|{row['partner']}" if row["direction"] == "sent"
            else f"{row['partner']}|{sheet_name}",
            axis=1
        )

        # Drop unneeded columns
        df.drop(columns=["direction", "partner"], inplace=True)

        # Move cluster_pair to first column
        cols = ["cluster_pair"] + [col for col in df.columns if col != "cluster_pair"]
        df = df[cols]

        all_data.append(df)

    # Merge all sheets
    merged_df = pd.concat(all_data, ignore_index=True)

    # Sort by cluster_pair first, then interacting_pair
    merged_df = merged_df.sort_values(by=["cluster_pair", "interacting_pair"], ascending=True)

    # Save result
    merged_df.to_excel(output_path, index=False)
    print(f"✅ Processed file saved to: {output_path}")


def start(n_proc: int = None) -> None:
    import pandas as pd

##
    # enrich_simplified_excel_with_genes(
    # simplified_excel_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list/simplified_significant_interactions_15.xlsx",
    # output_excel_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/significant_interactions_15_with_genes.xlsx",
    # cpdb_gene_input_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database/v4.1.0/gene_input.csv",
    # cpdb_complex_input_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database/v4.1.0/complex_input.csv",
    # cpdb_stat_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_significant_means_final_merged_injured_15_nona.txt"
    # )
    # enrich_simplified_excel_with_genes(
    # simplified_excel_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list/simplified_significant_interactions_60.xlsx",
    # output_excel_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/significant_interactions_60_with_genes.xlsx",
    # cpdb_gene_input_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database/v4.1.0/gene_input.csv",
    # cpdb_complex_input_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database/v4.1.0/complex_input.csv",
    # cpdb_stat_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_significant_means_final_merged_injured_60_nona.txt"
    # )

    alphabetic_order_int(
    filepath="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/simple_significant_interactions_15_with_genes.xlsx",
    output_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/merged_sorted_interactions_15.xlsx"
    )

    alphabetic_order_int(
    filepath="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/simple_significant_interactions_60_with_genes.xlsx",
    output_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/merged_sorted_interactions_60.xlsx"
    )

    alphabetic_order_cluster(
    filepath="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/simple_significant_interactions_15_with_genes.xlsx",
    output_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/merged_interactions_15.xlsx"
    )

    alphabetic_order_cluster(
    filepath="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/simple_significant_interactions_60_with_genes.xlsx",
    output_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/merged_interactions_60.xlsx"
    )




    

# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=os.cpu_count() - 1)
    
        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
    