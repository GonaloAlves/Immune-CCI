# Import necessary packages
import os
import gc
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# This script will display the gene expression values of the canonical genes from Immune dataset



# Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    print("Loading h5ad file")
    return sc.read_h5ad(file_path)


def imm_remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    
    mask_NA = adata.obs['leiden_fusion'] != 'Imm.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2

def mev_remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    
    mask_NA = adata.obs['leiden_fusion'] != 'MeV.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2

def neu_remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    
    mask_NA = adata.obs['leiden_fusion'] != 'Neu.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2


# Load canonical genes
def load_canonical_from_dir(directory):
    """
    Load gene lists from all .txt files in a directory.

    Parameters:
    directory (str): Path to the directory containing .txt files.

    Returns:
    dict: Dictionary of gene lists with filenames (without extensions) as keys.
    """
    print(f"Loading canonical gene lists from: {directory}")
    gene_files = [f for f in os.listdir(directory) if f.endswith('.txt')]

    gene_dict = {}
    for gene_file in gene_files:
        file_path = os.path.join(directory, gene_file)
        gene_name = os.path.splitext(gene_file)[0]  # Use the filename without the extension
        with open(file_path, 'r') as f:
            gene_dict[gene_name] = f.read().splitlines()

    print(f"Loaded gene lists: {list(gene_dict.keys())}")
    return gene_dict

def create_dotplots_with_thresholds(adata, genes, thresholds, output_dir, user_order, order_txt, name, order_name=None):
    """
    Create and save dotplots for different pts thresholds, with and without dendrograms.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    genes (dict): Dictionary of canonical genes grouped by gene groups.
    thresholds (list of float): List of pts thresholds to generate dotplots for.
    output_dir (str): Directory to save the dotplots.

    Returns:
    None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    adata = imm_keep_only_selected_clusters(adata=adata, clusters_to_keep=user_order)

    # Ensure leiden_fusion is categorical and reorder it
    #adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype('category')
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.reorder_categories(user_order, ordered=True)

    # ✅ Format gene names right away
    genes = format_gene_names(genes)
    
    print(genes)

    for threshold in thresholds:
        print(f"\nProcessing pts threshold: {threshold}")

        # Extract DGE data
        gene_names, pts = extract_dge_data(adata)

        # Create cluster DataFrames with the current threshold
        cluster_dfs = create_cluster_dfs(gene_names, pts, pts_threshold=threshold)

        # Compare canonical genes with cluster-specific genes
        filtered_genes = compare_canonical(genes, cluster_dfs)

        # Aggregate filtered genes by gene group
        top_genes_names = top_gene_names(filtered_genes, genes)

        # Example user-defined gene group order
        user_gene_group_order = order_txt#use the order that is down

        # Reorder the dictionary based on user order
        top_genes_names = {key: top_genes_names[key] for key in user_gene_group_order}

        # Generate four different dotplots per threshold
        print(f"Generating dotplots for pts threshold: {threshold}")

        # (1) Scaled expression (Greys) without dendrogram
        dotplot_scaled_no_dendro = sc.pl.dotplot(
            adata,
            var_names=top_genes_names,
            groupby='leiden_fusion',
            cmap='Greys',
            colorbar_title='Scaled expression',
            use_raw=False,
            standard_scale='var',
            dendrogram=False,
            return_fig=True
        )

        # (3) Raw expression (Reds) without dendrogram
        dotplot_normal_no_dendro = sc.pl.dotplot(
            adata,
            var_names=top_genes_names,
            groupby='leiden_fusion',
            cmap='Reds',
            use_raw=False,
            dendrogram=False,
            return_fig=True
        )

        # Save dotplots with appropriate filenames
        output_scaled_no_dendro = os.path.join(output_dir, f"{name}_{order_name}_dotplot_scaled_{threshold}.png")
        output_normal_no_dendro = os.path.join(output_dir, f"{name}_{order_name}_dotplot_normal_{threshold}.png")

        dotplot_scaled_no_dendro.savefig(output_scaled_no_dendro, bbox_inches="tight")
        dotplot_normal_no_dendro.savefig(output_normal_no_dendro, bbox_inches="tight")

        
        plt.close()
        gc.collect()
        print(f"Saved dotplots for threshold {threshold}:")
        print(f"  - {output_scaled_no_dendro}")


def extract_dge_data(adata):
    """
    Extract the different type of data from the AnnData object like the differential gene expression values from each cell.

    Parameters:
    adata (AnnData): The AnnData object containing the DGE data.

    Returns:
    tuple: DataFrames for gene names and pts.
    """
    dge_fusion = adata.uns['rank_genes_groups_leiden_fusion']
    
    # Convert the extracted data into DataFrames
    gene_names = pd.DataFrame(dge_fusion['names'])
    pts = pd.DataFrame(dge_fusion['pts'])
    
    return gene_names, pts



# Step 3: Create cluster dataframes with filtered data
def create_cluster_dfs(gene_names, pts, pts_threshold):
    """
    Create a dictionary of dataframes per cluster, with the respective gene expression data.
    By default this function will not order the genes by fold change and the default minimun pts is 0

    Parameters:
    gene_names (DataFrame): DataFrame with gene names.
    pts (DataFrame): DataFrame with p values for each gene from each cell.
    pts_threshold (float): The minimum percentage of expression value to filter genes.

    Pts: Proportion of cells in each group that express a particular gene

    Returns:
    dict: Dictionary of DataFrames for each cluster.
    """
    cluster_dfs = {} # creates dictionary

    for i in range(len(gene_names.columns)):  # For to read each cluster
        gene_reindex = gene_names.iloc[:, i]  # Takes correct order of the genes index of the gene names
        pts_reindexed = pts.iloc[:, i].reindex(gene_reindex.values)  # Reindex the pts with the gene names index

        # Create dataframe for each cluster
        cluster_df = pd.DataFrame({
            'pts': pts_reindexed.values},
            index=gene_reindex.values
        )

        # Mask to remove mitochondrial genes (those starting with "mt-")
        mask_mt = ~gene_reindex.str.startswith("mt-")  # Create mask where gene names don't start with "mt-"
        filtered_gene_reindex = gene_reindex[mask_mt]  # Filter out mitochondrial genes in the gene names
        cluster_df = cluster_df.loc[filtered_gene_reindex]  # Apply mask to the DataFrame
        
        # Filter by 'pts' using the user-specified threshold
        mask_pts = (cluster_df['pts'] >= pts_threshold)
        filtered_df = cluster_df[mask_pts]


        # Assigns the cluster name as the name of the resolution atributed
        cluster_names = gene_names.columns[i]
        cluster_dfs[cluster_names] = filtered_df  # Store the respective clusters in the dictionary


    return cluster_dfs


# Step 6: Compare canonical genes
def compare_canonical(genes, cluster_dfs):
    new_dic = {}
    for group_name, canonical_genes in genes.items():
        print(f"Processing canonical gene group: {group_name}")
        group_results = {}
        for cluster_name, cluster_df in cluster_dfs.items():
            filtered_cluster = cluster_df[cluster_df.index.isin(canonical_genes)]
            if not filtered_cluster.empty:
                group_results[cluster_name] = filtered_cluster
        new_dic[group_name] = group_results
    return new_dic



def top_gene_names(filtered_genes, original_gene_dict):
    """
    Aggregate the remaining filtered genes for each gene group while maintaining the original order.

    Parameters:
    filtered_genes (dict): Dictionary containing filtered genes per group and cluster.
    original_gene_dict (dict): Dictionary of gene groups as they appear in the original .txt files.

    Returns:
    dict: Dictionary where the keys are gene group names (txt file titles),
          and the values are lists of remaining genes, maintaining their original order.
    """
    top_genes_names = {}

    for group_name, clusters in filtered_genes.items():
        # Combine gene names across all clusters for this group while maintaining order
        combined_genes = set()
        for cluster_name, df in clusters.items():
            combined_genes.update(df.index.tolist())  # Store unique genes from all clusters
        
        # Preserve original order of genes from the text file
        ordered_genes = [gene for gene in original_gene_dict[group_name] if gene in combined_genes]
        top_genes_names[group_name] = ordered_genes  # Maintain original order

    return top_genes_names

def imm_keep_only_selected_clusters(adata: sc.AnnData, clusters_to_keep: list):
    """
    Keep only cells whose 'leiden_fusion' is in the clusters_to_keep list.

    Parameters:
    adata (AnnData): The AnnData object.
    clusters_to_keep (list): List of cluster names to keep.

    Returns:
    AnnData: Filtered AnnData object containing only the selected clusters.
    """
    print("Keeping only selected clusters")
    
    mask = adata.obs['leiden_fusion'].isin(clusters_to_keep)
    filtered_adata = adata[mask].copy()
    
    return filtered_adata

def format_gene_names(genes_dict):
    """
    Format all gene names in the given dictionary so only the first letter is uppercase, rest lowercase.
    
    Parameters:
    genes_dict (dict): Dictionary like {'group1': ['CADM1', 'IL6R', ...], ...}
    
    Returns:
    dict: Same structure, but with gene names formatted.
    """
    formatted_genes = {}
    for group, gene_list in genes_dict.items():
        formatted_genes[group] = [gene.capitalize() for gene in gene_list]
    return formatted_genes


# Main execution block
if __name__ == "__main__":
    # Load data
    adataimm = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad")
    adatamev = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")
    adataneu = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Neu_CentralCanal_raw_norm_ranked_copy_copy.h5ad")
    adatamerged = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_merged_raw_norm_annot_nona_copy.h5ad")


    
    filtered_adataimm = imm_remove_NA_cat(adataimm)
    filtered_adatamev = mev_remove_NA_cat(adatamev)
    filtered_adataneu = neu_remove_NA_cat(adataneu)

    # Load canonical gene lists from a directory

    # IMM receive(receptors) genes
    can_dir_imm_rec_15_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_rec/15days_all"
    can_dir_imm_rec_15_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_rec/15days_seperate"
    can_dir_imm_rec_60_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_rec/60days_all"
    can_dir_imm_rec_60_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_rec/60days_seperate"

    # IMM send(ligands) genes
    can_dir_imm_send_15_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_send/15days_all"
    can_dir_imm_send_15_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_send/15days_seperate"
    can_dir_imm_send_60_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_send/60days_all"
    can_dir_imm_send_60_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Immune/Immune_genes/Immune_genes_send/60days_seperate"

    ####
    # MeV receive(receptors) genes
    can_dir_mev_rec_15_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_rec/15days_all"
    can_dir_mev_rec_15_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_rec/15days_seperate"
    can_dir_mev_rec_60_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_rec/60days_all"
    can_dir_mev_rec_60_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_rec/60days_seperate"

    # MeV send(ligands) genes
    can_dir_mev_send_15_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_send/15days_all"
    can_dir_mev_send_15_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_send/15days_seperate"
    can_dir_mev_send_60_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_send/60days_all"
    can_dir_mev_send_60_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal/MeV_genes/MeV_genes_send/60days_seperate"

    ####
    # Neu receive(receptors) genes
    can_dir_neu_rec_15_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_rec/15days_all"
    can_dir_neu_rec_15_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_rec/15days_seperate"
    can_dir_neu_rec_60_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_rec/60days_all"
    can_dir_neu_rec_60_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_rec/60days_seperate"

    # Neu send(ligands) genes
    can_dir_neu_send_15_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_send/15days_all"
    can_dir_neu_send_15_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_send/15days_seperate"
    can_dir_neu_send_60_all = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_send/60days_all"
    can_dir_neu_send_60_seperate = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Neuron/Neu_genes/Neu_genes_send/60days_seperate"

    can_dir_mev_send_custom = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/txt"

    can_dir_custom_15d_coll = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/txt/15d_coll"
    can_dir_custom_60d_coll = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/txt/60d_coll"
    can_dir_custom_15d_endo = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/txt/15d_endo"
    can_dir_custom_15d_epend = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/txt/15d_epend"
    can_dir_custom_60d_epend = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/txt/60d_epend"

    
    recgenesimm15all = load_canonical_from_dir(can_dir_imm_rec_15_all)
    recgenesimm15sep = load_canonical_from_dir(can_dir_imm_rec_15_seperate)
    recgenesimm60all = load_canonical_from_dir(can_dir_imm_rec_60_all)
    recgenesimm60spe = load_canonical_from_dir(can_dir_imm_rec_60_seperate)

    sendgenesimm15all = load_canonical_from_dir(can_dir_imm_send_15_all)
    sendgenesimm15sep = load_canonical_from_dir(can_dir_imm_send_15_seperate)
    sendgenesimm60all = load_canonical_from_dir(can_dir_imm_send_60_all)
    sendgenesimm60spe = load_canonical_from_dir(can_dir_imm_send_60_seperate)


    recgenesmev15all = load_canonical_from_dir(can_dir_mev_rec_15_all)
    recgenesmev15sep = load_canonical_from_dir(can_dir_mev_rec_15_seperate)
    recgenesmev60all = load_canonical_from_dir(can_dir_mev_rec_60_all)
    recgenesmev60spe = load_canonical_from_dir(can_dir_mev_rec_60_seperate)

    sendgenesmev15all = load_canonical_from_dir(can_dir_mev_send_15_all)
    sendgenesmev15sep = load_canonical_from_dir(can_dir_mev_send_15_seperate)
    sendgenesmev60all = load_canonical_from_dir(can_dir_mev_send_60_all)
    sendgenesmev60spe = load_canonical_from_dir(can_dir_mev_send_60_seperate)


    recgenesneu15all = load_canonical_from_dir(can_dir_neu_rec_15_all)
    recgenesneu15sep = load_canonical_from_dir(can_dir_neu_rec_15_seperate)
    recgenesneu60all = load_canonical_from_dir(can_dir_neu_rec_60_all)
    recgenesneu60spe = load_canonical_from_dir(can_dir_neu_rec_60_seperate)

    sendgenesneu15all = load_canonical_from_dir(can_dir_neu_send_15_all)
    sendgenesneu15sep = load_canonical_from_dir(can_dir_neu_send_15_seperate)
    sendgenesneu60all = load_canonical_from_dir(can_dir_neu_send_60_all)
    sendgenesneu60spe = load_canonical_from_dir(can_dir_neu_send_60_seperate)

    sendgenesneucustom = load_canonical_from_dir(can_dir_mev_send_custom)

    custom15d_coll = load_canonical_from_dir(can_dir_custom_15d_coll)
    custom60d_coll = load_canonical_from_dir(can_dir_custom_60d_coll)
    custom15d_endo = load_canonical_from_dir(can_dir_custom_15d_endo)
    custom15d_epend = load_canonical_from_dir(can_dir_custom_15d_epend)
    custom60d_epend= load_canonical_from_dir(can_dir_custom_60d_epend)



    # # Create dendogram ot the top genes
    # dendogram_sc(filtered_adata)

    # Define thresholds
    pts_thresholds = [0, 0.2]

####

    imm_custom_cluster_order = [
    "Imm.M0Like.0", "Imm.M0Like.1", 
    "Imm.M0Like.2", "Imm.MHCII.0" ,
    "Imm.Interferon.0", "Imm.DAM.0", 
    "Imm.DAM.1", "Imm.PVM.0", "Imm.Proliferative.0"]

    mev_custom_cluster_order = ["MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Epithelial.0",
                            "MeV.SMC.0", "MeV.Pericytes.0", "MeV.VLMC.0", "MeV.VLMC.1" , "MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3",
                            "MeV.FibLaminin.0", "MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2", "MeV.Fib.5", "MeV.Fib.3", "MeV.Fib.4", "MeV.FibProlif.0"]
    
    neu_custom_cluster_order = ["Neu.CSFcN.0", "Neu.Epend.0"]

    imm_cluster_order_sigs = ["Imm.M0Like.1", 
    "Imm.Interferon.0", "Imm.DAM.0", "Imm.DAM.1", 
    "Imm.PVM.0"]

    mev_cluster_order_sigs = ["MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2",
                              "MeV.Pericytes.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3",
                              "MeV.Fib.4", "MeV.Fib.5"]
    
    neu_cluster_order_sigs = ["Neu.CSFcN.0", "Neu.Epend.0"]

    merged_custom_cluster_order = ["Imm.M0Like.1", 
    "Imm.Interferon.0", "Imm.DAM.0", "Imm.DAM.1", "Imm.PVM.0",
    
    "MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2",
    "MeV.Pericytes.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3",
    "MeV.Fib.4", "MeV.Fib.5",
    
    "Neu.Epend.0"]

    rec_cluster_remove_imm_15 = [
    "Imm.M0Like.1", 
    "Imm.Interferon.0",
    "Imm.PVM.0"]

    send_cluster_remove_imm_15 = [ 
    "Imm.Interferon.0", "Imm.DAM.0", 
    "Imm.PVM.0"]

    rec_cluster_remove_imm_60 = [
    "Imm.M0Like.1", 
    "Imm.Interferon.0", 
    "Imm.DAM.1"]

    send_cluster_remove_imm_60 = [
    "Imm.M0Like.1", 
    "Imm.Interferon.0", "Imm.DAM.0", 
    "Imm.DAM.1"]

####

    rec_cluster_remove_mev_15 = [
    "MeV.Endothelial.0",
    "MeV.Endothelial.1", 
    "MeV.Endothelial.2",
    "MeV.Pericytes.0"]

    send_cluster_remove_mev_15 = [
    "MeV.Endothelial.2", "MeV.FibCollagen.1", 
    "MeV.FibCollagen.2", "MeV.FibCollagen.3"]

    rec_cluster_remove_mev_60 = ["MeV.Endothelial.0", "MeV.Pericytes.0"]

    send_cluster_remove_mev_60 = [
    "MeV.FibCollagen.1", 
    "MeV.FibCollagen.3"]

####

    rec_cluster_remove_neu_15 = [
    "Neu.Epend.0"]

    rec_cluster_remove_neu_60 = [
    "Neu.Epend.0"]

    send_cluster_remove_neu_15 = [
    "Neu.Epend.0"]

    send_cluster_remove_neu_60 = [
    "Neu.Epend.0"]

#### Customs orders

    send_collagen_order_15 = [
        "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"
    ]

    send_collagen_order_60 = [
        "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"
    ]

    send_custom_collagen_order = ["MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"]


    #### Name of the list of genes txts

    immune_genes_rec_15_all = ["receptors|Interferon", "receptors|M0.1", "receptors|PVM"]

    immune_genes_rec_15_seperate = ["FibColl.2|Interferon", "FibColl.2|M0.1", "FibColl.3|PVM"]
    
    immune_genes_rec_60_all = ["receptors|DAM.1", "receptors|Intreferon", "receptors|M0.1"]

    immune_genes_rec_60_seperate = ["FibColl.3|DAM.1", "FibColl.3|Intreferon", "FibColl.3|M0.1"]
##
    immune_genes_send_15_all = ["DAM.0|ligands", "Intreferon|ligands", "PVM|ligands"]

    immune_genes_send_15_seperate = ["DAM.0|Endo.2", "DAM.0|Epend", "Intreferon|Endo.1", "PVM|Epend"]
    
    immune_genes_send_60_all = ["DAM|ligands", "Intreferon|ligands", "M0.1|ligands"]

    immune_genes_send_60_seperate = ["DAM.0|Epend", "DAM.1|Epend", "Intreferon|Epend", "M0.1|Epend"]
###
    meningeal_genes_rec_15_all = ["receptors|Endo.0", "receptors|Endo.1", "receptors|Endo.2", "receptors|Pericytes"]

    meningeal_genes_rec_15_seperate = ["FibColl.3|Endo.2", "FibColl.1|Endo.0", "DAM.0|Endo.2", "Interferon|Endo.1", "FibColl.1|Pericytes"]
    
    meningeal_genes_rec_60_all = ["receptors|Endo.0", "receptors|Pericytes"]

    meningeal_genes_rec_60_seperate = ["FibColl.3|Endo.0", "FibColl.1|Pericytes"]
##
    meningeal_genes_send_15_all = ["Endo.2|ligands", "FibColl.1|ligands", "FibColl.2|ligands", "FibColl.3|ligands"]

    meningeal_genes_send_15_seperate = ["FibColl.1|Endo.0", "FibColl.1|Pericytes", "FibColl.2|Interferon", "FibColl.2|M0.1", "FibColl.3|Endo.2", "FibColl.3|PVM"]
    
    meningeal_genes_send_60_all = ["FibColl.1|ligands", "FibColl.3|ligands"]

    meningeal_genes_send_60_seperate = ["FibColl.1|Pericytes", "FibColl.3|DAM.1", "FibColl.3|Endo.0", "FibColl.3|Epend", "FibColl.3|Interferon", "FibColl.3|M0.1"]
###
    neu_genes_rec_15_all = ["receptors|Epend"]

    neu_genes_rec_15_seperate = ["DAM.0|Epend", "PVM|Epend"]
    
    neu_genes_rec_60_all = ["receptors|Epend"]

    neu_genes_rec_60_seperate = ("DAM.0|Epend", "DAM.1|Epend", "FibColl.3|Epend", "Interferon|Epend", "M0.1|Epend")
##
    neu_genes_send_15_all_none = []

    neu_genes_send_15_seperate_none = []
    
    neu_genes_send_60_all_none = []

    neu_genes_send_60_seperate_none = []


    mev_custome_genes_coll = ["Fibroblasts", "ECM_Laminin", "ECM_Collagen", "CPDB_BMP", "CPDB_Glutamate", "CPDB_Collagen"]

    d15_coll = ["FibColl.1|Endo.0_ligands", "FibColl.1|Endo.0_receptors", 
                "FibColl.1|Pericytes_ligands", "FibColl.1|Pericytes_receptors",
                "FibColl.2|Interferon_ligands", "FibColl.2|Interferon_receptors", 
                "FibColl.2|M0.1_ligands", "FibColl.2|M0.1_receptors",
                "FibColl.3|Endo.2_ligands", "FibColl.3|Endo.2_receptor",
                "FibColl.3|PVM_ligands", "FibColl.3|PVM_receptors"]
    
    d60_coll = ["FibColl.1|Pericytes_ligands", "FibColl.1|Pericytes_receptors",
                "FibColl.3|DAM.1_ligands", "FibColl.3|DAM.1_receptors",
                "FibColl.3|Endo.0_ligands", "FibColl.3|Endo.0_receptors",
                "FibColl.3|Epend_ligands", "FibColl.3|Epend_receptors",
                "FibColl.3|Interferon_ligands", "FibColl.3|Intreferon_receptors",
                "FibColl.3|M0.1_ligands", "FibColl.3|M0.1_receptors"]
    
    d15_endo = ["DAM.0|Endo.2_ligands", "DAM.0|Endo.2_receptor",
                "FibColl.1|Endo.0_ligands", "FibColl.1|Endo.0_receptor",
                "FibColl.3|Endo.2_ligands", "FibColl.3|Endo.2_receptor",
                "Intreferon|Endo.1_ligands", "Interferon|Endo.1_receptor"]
    
    d15_epend = ["DAM.0|Epend_ligands", "DAM.0|Epend_receptors", 
                 "PVM|Epend_ligands", "PVM|Epend_receptors"]
    
    d60_epend = ["DAM.0|Epend_ligands", "DAM.0|Epend_receptors",
                 "DAM.1|Epend_ligands", "DAM.1|Epend_receptors",
                 "Intreferon|Epend_ligands", "Intreferon|Epend_receptors", 
                 "M0.1|Epend_ligands", "M0.1|Epend_receptors"]

    output_dir_immune = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_immune/cpdb_genes"
    output_dir_meningeal = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_meningeal/cpdb_genes"
    output_dir_neu = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_neuron/cpdb_genes"

    output_dir_costum = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/output"

    output_15d_coll = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/output/15d_coll"
    output_60d_coll = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/output/60d_coll"
    output_15d_endo = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/output/15d_endo"
    output_15d_epend = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/output/15d_epend"
    output_60d_epend = "/home/makowlg/Documents/Immune-CCI/src/canonical/custom_cpdb/output/60d_epend"


    name1= "imm_rec_15_all"
    name2= "imm_rec_15_sep"
    name3= "imm_rec_60_all"
    name4= "imm_rec_60_sep"
    name5= "imm_send_15_all"
    name6= "imm_send_15_sep"
    name7= "imm_send_60_all"
    name8= "imm_send_60_sep"
    
    name9= "mev_rec_15_all"
    name10= "mev_rec_15_sep"
    name11= "mev_rec_60_all"
    name12= "mev_rec_60_sep"
    name13= "mev_send_15_all"
    name14= "mev_send_15_sep"
    name15= "mev_send_60_all"
    name16= "mev_send_60_sep"

    name17= "neu_rec_15_all"
    name18= "neu_rec_15_sep"
    name19= "neu_rec_60_all"
    name20= "neu_rec_60_sep"
    name21= "neu_send_15_all"
    name22= "neu_send_15_sep"
    name23= "neu_send_60_all"
    name24= "neu_send_60_sep"

    name25= "custom_colls"

    name26= "15days_collagen"
    name27= "60days_collagen"
    name28= "15days_endothelial"
    name29= "15days_ependymal"
    name30= "60days_ependymal"

    name999= "full"


   

    #print(filtered_adataimm2.obs['leiden_fusion'])

    # # Case1 (imm_rec_15_all)
    # print(name1)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_imm_15, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_15_all,
    #                                 name=name1)
    
    # # Case2 (imm_rec_15_sep)
    # print(name2)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_imm_15, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_15_seperate,
    #                                 name=name2)
    
    # # Case3 (imm_rec_60_all)
    # print(name3)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_imm_60, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_60_all,
    #                                 name=name3)
    
    # # Case4 (imm_rec_60_sep)
    # print(name4)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_imm_60, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_60_seperate,
    #                                 name=name4)


    # # Case5 (imm_send_15_all)
    # print(name5)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_imm_15, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_15_all,
    #                                 name=name5)
    
    # # Case6 (imm_send_15_sep)
    # print(name6)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_imm_15, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_15_seperate,
    #                                 name=name6)
    
    # # Case7 (imm_send_60_all)
    # print(name7)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_imm_60, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_60_all,
    #                                 name=name7)
    
    # # Case8 (imm_send_60_sep)
    # print(name8)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_imm_60, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_60_seperate,
    #                                 name=name8)

    # # Case9 (mev_rec_15_all)
    # print(name9)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_mev_15, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_15_all,
    #                                 name=name9)
    
    # # Case10 (mev_rec_15_sep)
    # print(name10)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_mev_15, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_15_seperate,
    #                                 name=name10)
    
    # # Case11 (mev_rec_60_all)
    # print(name11)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_mev_60, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_60_all,
    #                                 name=name11)
    
    # # Case12 (mev_rec_60_sep)
    # print(name12)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_mev_60, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_60_seperate,
    #                                 name=name12)


    # # Case13 (mev_send_15_all)
    # print(name13)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_mev_15, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_15_all,
    #                                 name=name13)
    
    # # Case14 (mev_send_15_sep)
    # print(name14)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_mev_15, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_15_seperate,
    #                                 name=name14)
    
    # # Case15 (mev_send_60_all)
    # print(name15)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_mev_60, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_60_all,
    #                                 name=name15)
    
    # # Case16 (mev_send_60_sep)
    # print(name16)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_cluster_remove_mev_60, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_60_seperate,
    #                                 name=name16)
    
    # # Case17 (neu_rec_15_all)
    # print(name17)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_neu_15, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_15_all,
    #                                 name=name17)
    
    # # Case18 (neu_rec_15_sep)
    # print(name18)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_neu_15, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_15_seperate,
    #                                 name=name18)
    
    # # Case19 (neu_rec_60_all)
    # print(name19)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_neu_60, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_60_all,
    #                                 name=name19)
    
    # # Case20 (neu_rec_60_sep)
    # print(name20)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=rec_cluster_remove_neu_60, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_60_seperate,
    #                                 name=name20)


    # # # Case21 (neu_send_15_all)
    # print(name21)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu15all, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_15, 
    # #                                 output_dir=output_dir_neu,
    # #                                 order_txt=neu_genes_send_15_all,
    # #                                 name=name21)
    
    # # # Case22 (neu_send_15_sep)
    # print(name22)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu15sep, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_15, 
    # #                                 output_dir=output_dir_neu,
    # #                                 order_txt=neu_genes_send_15_seperate,
    # #                                 name=name22)
    
    # # # Case23 (neu_send_60_all)
    # print(name23)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu60all, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_60, 
    # #                                 output_dir=output_dir_meningeal,
    # #                                 order_txt=neu_genes_send_60_all,
    # #                                 name=name23)
    
    # # # Case24 (neu_send_60_sep)
    # print(name24)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu60spe, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_60, 
    # #                                 output_dir=output_dir_neu,
    # #                                 order_txt=neu_genes_send_60_seperate,
    # #                                 name=name24)

    # # Case24 (costum collagen)
    # print(name25)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesneucustom, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=send_custom_collagen_order, 
    #                                 output_dir=output_dir_costum,
    #                                 order_txt=mev_custome_genes_coll,
    #                                 name=name25)
    

    # ##################################################


    # # Case1 (imm_rec_15_all)
    # print(name1)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_15_all,
    #                                 name=name1,
    #                                 order_name=name999)
    
    # # Case2 (imm_rec_15_sep)
    # print(name2)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_15_seperate,
    #                                 name=name2,
    #                                 order_name=name999)
    
    # # Case3 (imm_rec_60_all)
    # print(name3)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_60_all,
    #                                 name=name3,
    #                                 order_name=name999)
    
    # # Case4 (imm_rec_60_sep)
    # print(name4)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=recgenesimm60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_rec_60_seperate,
    #                                 name=name4,
    #                                 order_name=name999)


    # # Case5 (imm_send_15_all)
    # print(name5)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_15_all,
    #                                 name=name5,
    #                                 order_name=name999)
    
    # # Case6 (imm_send_15_sep)
    # print(name6)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_15_seperate,
    #                                 name=name6,
    #                                 order_name=name999)
    
    # # Case7 (imm_send_60_all)
    # print(name7)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_60_all,
    #                                 name=name7,
    #                                 order_name=name999)
    
    # # Case8 (imm_send_60_sep)
    # print(name8)
    # create_dotplots_with_thresholds(adata=filtered_adataimm, 
    #                                 genes=sendgenesimm60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=imm_cluster_order_sigs, 
    #                                 output_dir=output_dir_immune,
    #                                 order_txt=immune_genes_send_60_seperate,
    #                                 name=name8,
    #                                 order_name=name999)

    # # Case9 (mev_rec_15_all)
    # print(name9)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_15_all,
    #                                 name=name9,
    #                                 order_name=name999)
    
    # # Case10 (mev_rec_15_sep)
    # print(name10)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_15_seperate,
    #                                 name=name10,
    #                                 order_name=name999)
    
    # # Case11 (mev_rec_60_all)
    # print(name11)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_60_all,
    #                                 name=name11,
    #                                 order_name=name999)
    
    # # Case12 (mev_rec_60_sep)
    # print(name12)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=recgenesmev60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_rec_60_seperate,
    #                                 name=name12,
    #                                 order_name=name999) 


    # # Case13 (mev_send_15_all)
    # print(name13)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_15_all,
    #                                 name=name13,
    #                                 order_name=name999)
    
    # # Case14 (mev_send_15_sep)
    # print(name14)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_15_seperate,
    #                                 name=name14,
    #                                 order_name=name999)
    
    # # Case15 (mev_send_60_all)
    # print(name15)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_60_all,
    #                                 name=name15,
    #                                 order_name=name999)
    
    # # Case16 (mev_send_60_sep)
    # print(name16)
    # create_dotplots_with_thresholds(adata=filtered_adatamev, 
    #                                 genes=sendgenesmev60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=mev_cluster_order_sigs, 
    #                                 output_dir=output_dir_meningeal,
    #                                 order_txt=meningeal_genes_send_60_seperate,
    #                                 name=name16,
    #                                 order_name=name999)
    
    # # Case17 (neu_rec_15_all)
    # print(name17)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu15all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=neu_cluster_order_sigs, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_15_all,
    #                                 name=name17,
    #                                 order_name=name999)
    
    # # Case18 (neu_rec_15_sep)
    # print(name18)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu15sep, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=neu_cluster_order_sigs, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_15_seperate,
    #                                 name=name18,
    #                                 order_name=name999)
    
    # # Case19 (neu_rec_60_all)
    # print(name19)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu60all, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=neu_cluster_order_sigs, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_60_all,
    #                                 name=name19,
    #                                 order_name=name999)
    
    # # Case20 (neu_rec_60_sep)
    # print(name20)
    # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    #                                 genes=recgenesneu60spe, 
    #                                 thresholds=pts_thresholds, 
    #                                 user_order=neu_cluster_order_sigs, 
    #                                 output_dir=output_dir_neu,
    #                                 order_txt=neu_genes_rec_60_seperate,
    #                                 name=name20,
    #                                 order_name=name999)


    # # # Case21 (neu_send_15_all)
    # print(name21)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu15all, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_15, 
    # #                                 output_dir=output_dir_neu,
    # #                                 order_txt=neu_genes_send_15_all,
    # #                                 name=name21)
    
    # # # Case22 (neu_send_15_sep)
    # print(name22)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu15sep, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_15, 
    # #                                 output_dir=output_dir_neu,
    # #                                 order_txt=neu_genes_send_15_seperate,
    # #                                 name=name22)
    
    # # # Case23 (neu_send_60_all)
    # print(name23)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu60all, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_60, 
    # #                                 output_dir=output_dir_meningeal,
    # #                                 order_txt=neu_genes_send_60_all,
    # #                                 name=name23)
    
    # # # Case24 (neu_send_60_sep)
    # print(name24)
    # # create_dotplots_with_thresholds(adata=filtered_adataneu, 
    # #                                 genes=sendgenesneu60spe, 
    # #                                 thresholds=pts_thresholds, 
    # #                                 user_order=send_cluster_remove_neu_60, 
    # #                                 output_dir=output_dir_neu,
    # #                                 order_txt=neu_genes_send_60_seperate,
    # #                                 name=name24)

    # Case25 (costum collagen)
    print(name25)
    create_dotplots_with_thresholds(adata=filtered_adatamev, 
                                    genes=sendgenesneucustom, 
                                    thresholds=pts_thresholds, 
                                    user_order=mev_cluster_order_sigs, 
                                    output_dir=output_dir_costum,
                                    order_txt=mev_custome_genes_coll,
                                    name=name25,
                                    order_name=name999)
    
    # Case26 (costum 15_collagen)
    print(name26)
    create_dotplots_with_thresholds(adata=adatamerged, 
                                    genes=custom15d_coll, 
                                    thresholds=pts_thresholds, 
                                    user_order=merged_custom_cluster_order, 
                                    output_dir=output_15d_coll,
                                    order_txt=d15_coll,
                                    name=name26)
    
    # Case27 (costum 60_collagen)
    print(name27)
    create_dotplots_with_thresholds(adata=adatamerged, 
                                    genes=custom60d_coll, 
                                    thresholds=pts_thresholds, 
                                    user_order=merged_custom_cluster_order, 
                                    output_dir=output_60d_coll,
                                    order_txt=d60_coll,
                                    name=name27)
    
    # Case28 (costum 15_endo)
    print(name28)
    create_dotplots_with_thresholds(adata=adatamerged, 
                                    genes=custom15d_endo, 
                                    thresholds=pts_thresholds, 
                                    user_order=merged_custom_cluster_order, 
                                    output_dir=output_15d_endo,
                                    order_txt=d15_endo,
                                    name=name28)
    
    # Case29 (costum 15_epend)
    print(name29)
    create_dotplots_with_thresholds(adata=adatamerged, 
                                    genes=custom15d_epend, 
                                    thresholds=pts_thresholds, 
                                    user_order=merged_custom_cluster_order, 
                                    output_dir=output_15d_epend,
                                    order_txt=d15_epend,
                                    name=name29)
    
    # Case30 (costum 60_epend)
    print(name30)
    create_dotplots_with_thresholds(adata=adatamerged, 
                                    genes=custom60d_epend, 
                                    thresholds=pts_thresholds, 
                                    user_order=merged_custom_cluster_order, 
                                    output_dir=output_60d_epend,
                                    order_txt=d60_epend,
                                    name=name30)