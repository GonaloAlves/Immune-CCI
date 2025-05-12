import pandas as pd
import os

# List of significant clusters
significant_clusters = [
    "Imm.M0Like.1", "Imm.DAM.0", "Imm.Interferon.0", "Imm.PVM.0", "Imm.DAM.1",
    "Neu.CSFcN.0", "Neu.Epend.0", "MeV.Pericytes.0", "MeV.Endothelial.0", "MeV.FibCollagen.1", "MeV.Fib.5", "MeV.Fib.4", 
    "MeV.FibCollagen.2", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.FibCollagen.3"
]

# Paths
base_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb"
output_dir = os.path.join(base_dir, "filtered_pvalues")
os.makedirs(output_dir, exist_ok=True)

# Input files
input_files = {
    "control": "statistical_analysis_pvalues_final_merged_uninjured_nona.txt",
    "injured_15": "statistical_analysis_pvalues_final_merged_injured_15_nona.txt",
    "injured_60": "statistical_analysis_pvalues_final_merged_injured_60_nona.txt"
}

def is_valid_column(col, mode):
    if "|" not in col:
        return False
    parts = col.split("|")
    if mode == "both":
        return all(p in significant_clusters for p in parts)
    elif mode == "one":
        return any(p in significant_clusters for p in parts)
    return False

# Process
for label, filename in input_files.items():
    df = pd.read_csv(os.path.join(base_dir, filename), sep="\t")

    meta_cols = [col for col in df.columns if "|" not in col]
    interact_cols = [col for col in df.columns if "|" in col]

    # Filters
    both_cols = [col for col in interact_cols if is_valid_column(col, mode="both")]
    one_cols = [col for col in interact_cols if is_valid_column(col, mode="one")]

    df_both = df[meta_cols + both_cols]
    df_one = df[meta_cols + one_cols]

    df_both.to_csv(os.path.join(output_dir, f"{label}_only_significant_clusters.txt"), sep="\t", index=False)
    df_one.to_csv(os.path.join(output_dir, f"{label}_atleast_one_significant_cluster.txt"), sep="\t", index=False)

    print(f"{label} done: {len(both_cols)} kept (both), {len(one_cols)} kept (one)")
