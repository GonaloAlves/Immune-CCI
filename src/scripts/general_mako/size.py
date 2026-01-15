# Import necessary packages
import os
import scanpy as sc
import pandas as pd


def load_data(file_path):
    print(f"\n🔹 Loading h5ad file: {file_path}")
    adata = sc.read_h5ad(file_path)
    print(f"Data loaded successfully! Shape: {adata.shape}")
    return adata


def remove_NA_cat(adata: sc.AnnData, key="leiden_fusion_old1", na_value="Imm.NA"):
    print("Removing NA cells category")
    mask = adata.obs[key] != na_value
    return adata[mask].copy()


def summarize_injury_control_per_cluster(
    adata,
    cluster_key="leiden_fusion",
    cluster_order=None
):
    """
    control = sham + uninjured
    injured = injured (day 15 + 60)
    """

    df = adata.obs[[cluster_key, "injury", "day"]].copy()

    # Define groups
    df["group"] = "other"
    df.loc[df["injury"].isin(["sham", "uninjured"]), "group"] = "control"
    df.loc[
        (df["injury"] == "injured") & (df["day"].isin([15, 60])),
        "group"
    ] = "injured"

    # Keep only relevant cells
    df = df[df["group"].isin(["control", "injured"])]

    # Per-cluster counts
    count_table = (
        pd.crosstab(df[cluster_key], df["group"])
        .reindex(columns=["control", "injured"], fill_value=0)
    )

    # Per-cluster totals & percentages
    count_table["total_cells"] = count_table.sum(axis=1)
    count_table["control_pct"] = (
        count_table["control"] / count_table["total_cells"] * 100
    ).round(2)
    count_table["injured_pct"] = (
        count_table["injured"] / count_table["total_cells"] * 100
    ).round(2)

    # Reorder clusters if requested
    if cluster_order is not None:
        count_table = count_table.reindex(
            [c for c in cluster_order if c in count_table.index]
        )

    # 🔹 ADD GLOBAL TOTAL ROW
    total_control = count_table["control"].sum()
    total_injured = count_table["injured"].sum()
    total_cells = total_control + total_injured

    total_row = pd.DataFrame(
        {
            "control": [total_control],
            "injured": [total_injured],
            "total_cells": [total_cells],
            "control_pct": [round(total_control / total_cells * 100, 2)],
            "injured_pct": [round(total_injured / total_cells * 100, 2)],
        },
        index=["TOTAL"]
    )

    count_table = pd.concat([count_table, total_row])

    return count_table

DATASETS = {
    "Immune": {
        "path": "/home/makowlg/Documents/Immune-CCI/h5ad_files/"
                "adata_final_Immune_raw_norm_ranked_copy_copy.h5ad",
        "cluster_order": [
            "Imm.M0Like.0", "Imm.M0Like.1",
            "Imm.M0Like.2", "Imm.MHCII.0",
            "Imm.Interferon.0", "Imm.DAM.0",
            "Imm.DAM.1", "Imm.PVM.0", "Imm.Proliferative.0"
        ],
        "remove_na": True
    },

    "Meningeal": {
        "path": "/home/makowlg/Documents/Immune-CCI/h5ad_files/"
                "adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad",
        "cluster_order": [
            "MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2",
            "MeV.Endothelial.3", "MeV.Epithelial.0",
            "MeV.SMC.0", "MeV.Pericytes.0",
            "MeV.VLMC.0", "MeV.VLMC.1",
            "MeV.FibCollagen.0", "MeV.FibCollagen.1",
            "MeV.FibCollagen.2", "MeV.FibCollagen.3",
            "MeV.FibLaminin.0",
            "MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2",
            "MeV.Fib.5", "MeV.Fib.3", "MeV.Fib.4",
            "MeV.FibProlif.0"
        ],
        "remove_na": False
    },

    "Neural": {
        "path": "/home/makowlg/Documents/Immune-CCI/h5ad_files/"
                "adata_final_Neu_CentralCanal_raw_norm_ranked_copy_copy.h5ad",
        "cluster_order": [
            "Neu.CSFcN.0", "Neu.Epend.0"
        ],
        "remove_na": False
    }
}


# Main execution block
if __name__ == "__main__":

    output_dir = "/home/makowlg/Documents/Immune-CCI/src/fractions_related"
    os.makedirs(output_dir, exist_ok=True)

    for name, cfg in DATASETS.items():
        print(f"\n===== Processing {name} dataset =====")

        adata = load_data(cfg["path"])

        if cfg.get("remove_na", False):
            adata = remove_NA_cat(adata)

        summary = summarize_injury_control_per_cluster(
            adata,
            cluster_key="leiden_fusion",
            cluster_order=cfg["cluster_order"]
        )

        print(summary)

        out_path = os.path.join(
            output_dir,
            f"{name.lower()}_injury_control_per_cluster.tsv"
        )
        summary.to_csv(out_path, sep="\t")

        print(f"✅ Saved: {out_path}")

