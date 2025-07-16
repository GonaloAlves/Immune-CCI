# This script will take the significant interations from cellphonedb output and its gonna be removing the interations from the uninjered/sham file 
# which will only be left with the injured new ones 
#
# Goncalo Alves Msc Thesis 2025

import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict
import os



cellphonedb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/heatmaps/cutted"


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



def build_cluster_interaction_matrix(filtered_dict):
    """
    Builds a symmetric cluster-cluster interaction matrix from the interaction dictionary.
    """
    interaction_counts = defaultdict(lambda: defaultdict(int))
    all_clusters = set()

    for interaction_id, cluster_pairs in filtered_dict.items():
        for pair in cluster_pairs:
            clusterA, clusterB = pair.split('|')
            all_clusters.update([clusterA, clusterB])

            # Symmetric matrix: increment both directions, unless it's a self-interaction
            interaction_counts[clusterA][clusterB] += 1
            if clusterA != clusterB:
                interaction_counts[clusterB][clusterA] += 1

    # Create a sorted list of unique clusters for the index/columns
    sorted_clusters = sorted(all_clusters)
    matrix = pd.DataFrame(0, index=sorted_clusters, columns=sorted_clusters)

    for clusterA in sorted_clusters:
        for clusterB in sorted_clusters:
            matrix.at[clusterA, clusterB] = interaction_counts[clusterA][clusterB]

    return matrix



def test_heatmap(category: str = None, remove_clusters: list = [], matrix: pd.DataFrame = None , vmin: int = None, vmax: int = None):
 
    custom_order = [
        "Imm.M0Like.0", "Imm.M0Like.1", "Imm.M0Like.2", "Imm.MHCII.0", "Imm.Interferon.0","Imm.DAM.0", "Imm.DAM.1", "Imm.PVM.0", "Imm.Proliferative.0",
        "Neu.CSFcN.0", "Neu.Epend.0", "MeV.Endothelial.0", "MeV.Endothelial.1","MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Epithelial.0", "MeV.SMC.0", "MeV.Pericytes.0",
        "MeV.VLMC.0", "MeV.VLMC.1", "MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3", "MeV.FibLaminin.0", 
        "MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2" ,"MeV.Fib.3", "MeV.FibProlif.0"
    ]

    # Check if all your custom labels exist in the matrix
    missing = set(custom_order) - set(matrix.index)
    if missing:
        print("Warning: These cluster names are not in the matrix:", missing)

    # Remove clusters you want to exclude from the matrix
    if remove_clusters:
        print(f"Removing the following clusters: {remove_clusters}")
        matrix = matrix.drop(index=remove_clusters, columns=remove_clusters, errors='ignore')

    # Reorder the matrix rows and columns (if custom_order still applies)
    ordered_matrix = matrix.loc[custom_order, custom_order]
    print(f"Ordered Matrix:\n{ordered_matrix}")

    #show only one part
    mask = np.tril(np.ones(ordered_matrix.shape, dtype=bool), k=-1) 

    # print("######")
    # print(ordered_matrix.shape)
    # print("######")

    plt.figure(figsize=(60, 60))

    # Define your custom gradient colors: Blue → Light Beige → Purplish-Red
    custom_colors = ["#2166ac", "#ffead0", "#b2182b"]  # Replace with exact hex if needed

    # Create the colormap
    custom_cmap = LinearSegmentedColormap.from_list("custom_bluered", custom_colors, N=256)

    ax = sns.heatmap(
        ordered_matrix,
        mask = mask,
        annot=True,
        fmt=".0f",
        cmap=custom_cmap,
        linewidths=0.2,
        linecolor='gray',
        square=True,
        cbar_kws={"shrink": 0.8},
        xticklabels=True,
        yticklabels=True,
        annot_kws={"size": 30},
        vmin=vmin,
        vmax=vmax
    )

    # Move x-tick labels to top
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=40, rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=40, rotation=0)
    ax.yaxis.set_ticks_position('right')

    # Add title at the bottom manually
    plt.figtext(0.5, 0.01, f"Number of Significant Interactions in {category}", ha='center', fontsize=60)

    # Adjust colorbar ticksM0Like.1
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    cbar.ax.set_position([0.85, 0.2, 0.5, 0.3])  # [left, bottom, width, height]

    # Save the figure
    output_path = f"{cellphonedb_dir}/manual_filtered_heatmap_{category}_70.png"
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()


def export_to_excel_inverted(interaction_dict, output_path):
    """
    Inverts the interaction dictionary and exports it to Excel with sorted columns.
    
    - Keys of the original dict are interaction IDs.
    - Values are lists of cluster pairs (e.g., 'clusterA|clusterB').

    The output Excel will have:
    - Columns = sorted unique cluster pairs
    - Rows = interaction IDs where that pair occurred
    """
    from collections import defaultdict
    import pandas as pd

    # Step 1: Invert the dictionary
    inverted = defaultdict(list)
    for interaction_id, pairs in interaction_dict.items():
        for pair in pairs:
            inverted[pair].append(interaction_id)

    # Step 2: Sort the cluster pairs alphabetically
    sorted_pairs = sorted(inverted.keys())

    # Step 3: Create a DataFrame with equal-length columns
    max_rows = max(len(inverted[pair]) for pair in sorted_pairs)
    inverted_df = pd.DataFrame({
        pair: inverted[pair] + [None] * (max_rows - len(inverted[pair])) for pair in sorted_pairs
    })

    # Step 4: Export to Excel
    try:
        inverted_df.to_excel(output_path, index=False)
        print(f"Inverted interaction matrix exported successfully to {output_path}")
    except Exception as e:
        print(f"Failed to export inverted DataFrame: {e}")

def export_detailed_excel_inverted(interaction_dict, pval_df_path, output_path):
    """
    Export an Excel file where:
    - Each cluster|cluster pair gets two columns:
        * One for interaction IDs
        * One for corresponding interacting_pair names
    - Row i of column A is an ID
    - Row i of column B is its corresponding name
    """

    import pandas as pd
    from collections import defaultdict

    # Load p-value DataFrame and map ID to interacting_pair
    full_df = pd.read_csv(pval_df_path, sep='\t')
    interaction_pair_map = full_df.set_index("id_cp_interaction")["interacting_pair"].to_dict()

    # Invert interaction_dict to group by cluster pairs
    from collections import defaultdict
    inverted = defaultdict(list)
    for interaction_id, pairs in interaction_dict.items():
        for pair in pairs:
            inverted[pair].append(interaction_id)

    # Sort cluster pairs alphabetically
    sorted_pairs = sorted(inverted.keys())

    # Prepare output structure
    output_data = {}

    for pair in sorted_pairs:
        ids = inverted[pair]
        id_col = []
        name_col = []

        for inter_id in ids:
            id_col.append(inter_id)
            name_col.append(interaction_pair_map.get(inter_id, "NA"))

        # Add columns for IDs and their corresponding names
        output_data[pair] = id_col
        output_data[f"{pair}_names"] = name_col

    # Normalize lengths
    max_len = max(len(col) for col in output_data.values())
    for k in output_data:
        output_data[k] += [None] * (max_len - len(output_data[k]))

    # Create DataFrame and save
    df_out = pd.DataFrame(output_data)
    df_out.to_excel(output_path, index=False)
    print(f"✅ Exported enriched Excel to {output_path}") 
    
def eexport_detailed_excel_inverted(interaction_dict, pval_df_path, output_path):
    """
    Export an Excel where:
    - Each column is a sorted cluster-cluster pair
    - Each column contains alternating rows: [interaction_id, interacting_pair]
    """
    import pandas as pd
    from collections import defaultdict

    # Load full p-value dataframe
    full_df = pd.read_csv(pval_df_path, sep='\t')

    # Map id_cp_interaction → interacting_pair
    interaction_pair_map = full_df.set_index("id_cp_interaction")["interacting_pair"].to_dict()

    # Invert the interaction_dict: cluster_pair → list of interaction_ids
    inverted = defaultdict(list)
    for interaction_id, pairs in interaction_dict.items():
        for pair in pairs:
            inverted[pair].append(interaction_id)

    # Sort the cluster pairs alphabetically
    sorted_pairs = sorted(inverted.keys())

    # Build the new table with alternating id/pair rows
    output_data = {}
    for pair in sorted_pairs:
        ids = inverted[pair]
        paired_list = []
        for inter_id in ids:
            paired_list.append(inter_id)
            paired_list.append(interaction_pair_map.get(inter_id, "NA"))
        output_data[pair] = paired_list

    # Normalize lengths
    max_len = max(len(col) for col in output_data.values())
    for k in output_data:
        output_data[k] += [None] * (max_len - len(output_data[k]))

    # Create DataFrame and export
    df_out = pd.DataFrame(output_data)
    df_out.to_excel(output_path, index=False)
    print(f"✅ Exported enriched Excel to {output_path}")


def build_directional_edge_list(filtered_dict):
    """
    Builds a directional interaction DataFrame (edge list) from the interaction dictionary.
    Each row represents a directed pair (from, to) with a count of interactions.
    
    Returns:
    - A DataFrame with columns: ['from', 'to', 'value']
    """
    from collections import defaultdict
    import pandas as pd

    interaction_counts = defaultdict(int)

    for interaction_id, cluster_pairs in filtered_dict.items():
        for pair in cluster_pairs:
            clusterA, clusterB = pair.split('|')
            interaction_counts[(clusterA, clusterB)] += 1

    # Convert to DataFrame
    edge_list = pd.DataFrame(
        [(a, b, count) for (a, b), count in interaction_counts.items()],
        columns=["from", "to", "value"]
    )

    return edge_list


def plot_interaction_distribution(edge_list_df, condition_label="", output_dir="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/histograms"):
    """
    Saves a histogram showing the distribution of interaction counts for a given condition.

    Parameters:
    - edge_list_df: pd.DataFrame with columns ['from', 'to', 'value']
    - condition_label: str, used in plot title and filename (e.g., "injured_15")
    - output_dir: str, directory where the histogram will be saved
    """
    if 'value' not in edge_list_df.columns:
        raise ValueError("Input DataFrame must contain a 'value' column.")
    
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    output_path = os.path.join(output_dir, f"interaction_distribution_{condition_label}.png")

    plt.figure(figsize=(10, 6))
    sns.histplot(edge_list_df['value'], bins=30, kde=True, color='skyblue', edgecolor='black')
    plt.title(f"Distribution of Interaction Counts - {condition_label}", fontsize=16)
    plt.xlabel("Number of Significant Interactions", fontsize=14)
    plt.ylabel("Number of Cluster Pairs", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"✅ Histogram saved to {output_path}")


def reorder_edge_list_by_groups(edge_list_df, group_dict, output_edge_path, output_group_path):
    """
    Reorders the edge list based on a user-defined biological grouping dictionary
    and exports both the reordered edge list and group annotations.

    Parameters:
    - edge_list_df: pd.DataFrame with 'from', 'to', 'value'
    - group_dict: dict, e.g., {"Immune": ["Imm.M0Like.1", ...], "Mesenchymal": [...], ...}
    - output_edge_path: path to save the reordered edge list
    - output_group_path: path to save the cluster-to-group mapping
    """

    # Flatten and order all clusters
    ordered_clusters = [cl for group in group_dict.values() for cl in group]

    # Filter only relevant cluster pairs (i.e., both in group_dict)
    included_clusters = set(ordered_clusters)
    filtered_df = edge_list_df[
        edge_list_df['from'].isin(included_clusters) & edge_list_df['to'].isin(included_clusters)
    ]

    # Enforce ordering on 'from' and 'to' for consistent plotting
    filtered_df['from'] = pd.Categorical(filtered_df['from'], categories=ordered_clusters, ordered=True)
    filtered_df['to'] = pd.Categorical(filtered_df['to'], categories=ordered_clusters, ordered=True)

    # Optional: sort edge list by order in groups
    filtered_df = filtered_df.sort_values(by=['from', 'to'])

    # Save edge list
    filtered_df.to_csv(output_edge_path, index=False)
    print(f"✅ Reordered edge list saved to: {output_edge_path}")

    # Build and export cluster → group mapping
    cluster_to_group = {}
    for group_name, clusters in group_dict.items():
        for cl in clusters:
            cluster_to_group[cl] = group_name

    group_df = pd.DataFrame(list(cluster_to_group.items()), columns=["cluster", "group"])
    group_df.to_csv(output_group_path, index=False)
    print(f"✅ Cluster-group mapping saved to: {output_group_path}")


def plot_interaction_distribution_matplotlib(edge_list_df, condition_label="", output_dir="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/histograms"):
    """
    Saves a histogram using Matplotlib to show the distribution of interaction counts
    with correct bar alignment and axis labels.

    Parameters:
    - edge_list_df: pd.DataFrame with a 'value' column
    - condition_label: str, label for the plot
    - output_dir: str, path to save the figure
    """
    if 'value' not in edge_list_df.columns:
        raise ValueError("Input DataFrame must contain a 'value' column.")
    
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"interaction_distribution_{condition_label}.png")

    values = edge_list_df['value']
    
    plt.figure(figsize=(10, 6))
    counts, bins, patches = plt.hist(values, bins=range(int(values.min()), int(values.max()) + 2), color='skyblue', edgecolor='black', align='left')

    plt.title(f"Distribution of Interaction Counts - {condition_label}", fontsize=16)
    plt.xlabel("Number of Significant Interactions", fontsize=14)
    plt.ylabel("Number of Cluster Pairs", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xticks(bins[:-1])  # show each bin value on x-axis

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"✅ Histogram saved to {output_path}")


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
    #print(filtered_15_dict)
    filtered_60_dict = filter_injured_by_control(control_dict, injured_60_dict, verbose=False)
    #print(filtered_60_dict)
    
    # from pprint import pprint
    # print("\nExample from Injured 15:")
    # pprint(list(injured_15_dict.items())[:3])

    matrix_15 = build_cluster_interaction_matrix(filtered_15_dict)
    matrix_60 = build_cluster_interaction_matrix(filtered_60_dict)

    edge_list_15 = build_directional_edge_list(filtered_15_dict)
    edge_list_60 = build_directional_edge_list(filtered_60_dict)

    # Optional: Save to CSV or Excel
    edge_list_15.to_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list_injured_15.csv", index=False)
    edge_list_60.to_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list_injured_60.csv", index=False)

    print(matrix_15)
    print(matrix_60)

    # plot_interaction_distribution(edge_list_15, condition_label="Injured 15 min")
    # plot_interaction_distribution(edge_list_60, condition_label="Injured 60 min")

    remove_clusters = ["MeV.EndoUnknow.4" ,"MeV.ImmuneDoublets.0", "MeV.FibUnknown.6", "MeV.LowQuality.0"]
    test_heatmap(category="injured_15", matrix = matrix_15 , remove_clusters=remove_clusters, vmin = 0, vmax = 70)
    test_heatmap(category="injured_60", matrix = matrix_60 , remove_clusters=remove_clusters, vmin = 0, vmax = 70)

    # export_detailed_excel_inverted(
    #     interaction_dict=filtered_15_dict,
    #     pval_df_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_injured_15_nona.txt",
    #     output_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/test_injured_15_enriched.xlsx"
    # )

    # export_detailed_excel_inverted(
    #     interaction_dict=filtered_60_dict,
    #     pval_df_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_pvalues_final_merged_injured_60_nona.txt",
    #     output_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/test_injured_60_enriched.xlsx"
    # )

    # export_to_excel_inverted(filtered_15_dict, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/injured_15_inverted.xlsx")
    # export_to_excel_inverted(filtered_60_dict, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/injured_60_inverted.xlsx")

    # export_to_excel(control_df, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/control_simplified.xlsx")
    # export_to_excel(injured_15_df, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/injured_15_simplified.xlsx")
    # export_to_excel(injured_60_df, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/injured_60_simplified.xlsx")

    biological_groups = {
    "Imm_Resting": ["Imm.M0Like.0", "Imm.M0Like.1", "Imm.M0Like.2"], 
    "Imm_Other": ["Imm.MHCII.0", "Imm.PVM.0"], 
    "Imm_Injury": ["Imm.Interferon.0", "Imm.DAM.0", "Imm.DAM.1"],

    "Neu": ["Neu.CSFcN.0", "Neu.Epend.0"],

    "MeV_Vascular": ["MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Pericytes.0", "MeV.SMC.0"],
    "MeV_Epithelial": ["MeV.Epithelial.0"],
    "MeV_Fibroblast": ["MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2", "MeV.Fib.3", "MeV.Fib.4", "MeV.Fib.5"],
    "MeV_Fib_Col": ["MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"],
    "MeV_Fib_Lam": ["MeV.FibLaminin.0"],

    "VLMC": ["MeV.VLMC.0", "MeV.VLMC.1"],

    "Prolifs": ["Imm.Proliferative.0","MeV.FibProlif.0"]
    }

    # # Run it
    # reorder_edge_list_by_groups(
    #     edge_list_df=edge_list_15,
    #     group_dict=biological_groups,
    #     output_edge_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_15.csv",
    #     output_group_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/group_annotation_15.csv"
    # )
    # reorder_edge_list_by_groups(
    #     edge_list_df=edge_list_60,
    #     group_dict=biological_groups,
    #     output_edge_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_60.csv",
    #     output_group_path="/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/group_annotation_60.csv"
    # )

    # plot_interaction_distribution_matplotlib(edge_list_15, condition_label="Injured 15")
    # plot_interaction_distribution_matplotlib(edge_list_60, condition_label="Injured 60")



