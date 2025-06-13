import pandas as pd
from collections import defaultdict
import seaborn as sns
import numpy as np
import os


output_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list"


def calculate_cluster_percentages(edge_list_path, output_prefix):
    """
    Calculates percentage of interactions sent (from clusters) and received (to clusters),
    and saves them as Excel files.

    Parameters:
        edge_list_path (str): Path to the CSV file containing 'from', 'to', and 'value' columns.
        output_prefix (str): Prefix for the output Excel file names.
    """
    df = pd.read_csv(edge_list_path)

    # ----- Percentage Sent (from cluster)
    from_totals = df.groupby("from")["value"].sum()
    total_sent = from_totals.sum()
    from_percentages = (from_totals / total_sent * 100).round(2).reset_index()
    from_percentages.columns = ["Cluster", "Percentage_Sent"]

    # ----- Percentage Received (to cluster)
    to_totals = df.groupby("to")["value"].sum()
    total_received = to_totals.sum()
    to_percentages = (to_totals / total_received * 100).round(2).reset_index()
    to_percentages.columns = ["Cluster", "Percentage_Received"]

    # ----- Save both to Excel
    from_output_path = os.path.join(output_dir, f"{output_prefix}_sent_percentages.xlsx")
    to_output_path = os.path.join(output_dir, f"{output_prefix}_rec_percentages.xlsx")

    from_percentages.to_excel(from_output_path, index=False)
    to_percentages.to_excel(to_output_path, index=False)

    print(f"Saved sent percentages to: {from_output_path}")
    print(f"Saved received percentages to: {to_output_path}")


def calculate_combined_cluster_influence(edge_list_path, output_prefix):
    """
    Calculates total interaction influence (sent + received %) per cluster
    and saves to Excel.

    Parameters:
        edge_list_path (str): Path to the CSV file with 'from', 'to', and 'value' columns.
        output_prefix (str): Prefix for the output Excel file.
    """
    df = pd.read_csv(edge_list_path)

    # Calculate % Sent
    from_totals = df.groupby("from")["value"].sum()
    total_sent = from_totals.sum()
    from_percentages = (from_totals / total_sent * 100).round(2).reset_index()
    from_percentages.columns = ["Cluster", "Percentage_Sent"]

    # Calculate % Received
    to_totals = df.groupby("to")["value"].sum()
    total_received = to_totals.sum()
    to_percentages = (to_totals / total_received * 100).round(2).reset_index()
    to_percentages.columns = ["Cluster", "Percentage_Received"]

    # Merge and fill missing clusters (e.g., cluster that only sends or only receives)
    combined = pd.merge(from_percentages, to_percentages, on="Cluster", how="outer").fillna(0)
    combined["Total_Influence"] = (combined["Percentage_Sent"] + combined["Percentage_Received"]).round(2)

    # Save
    combined_output_path = os.path.join(output_dir, f"{output_prefix}_combined_influence.xlsx")
    combined.to_excel(combined_output_path, index=False)
    print(f"Saved combined influence to: {combined_output_path}")


# Main execution block
if __name__ == "__main__":
    # Input paths
    injured_15 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/filtered_pvalues/grouped_edge_list_15.csv"
    injured_60 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/filtered_pvalues/grouped_edge_list_60.csv"


    calculate_cluster_percentages(injured_15, "cluster_percentages_15")
    calculate_cluster_percentages(injured_60, "cluster_percentages_60")


    calculate_combined_cluster_influence(injured_15, "cluster_combined_15")
    calculate_combined_cluster_influence(injured_60, "cluster_combined_60")
