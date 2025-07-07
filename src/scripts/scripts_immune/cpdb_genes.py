# Import necessary packages
import os
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


def remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    
    mask_NA = adata.obs['leiden_fusion'] != 'Imm.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2

# Step 4: Remove any cluster
def remove_clusters_by_suffix(cluster_dfs, suffix):
    """
    Remove clusters whose names end with a specific suffix corresponding to the resolution number, or if the case is the NA cluster, all are integrated in the dictionary with every cluster.

    Parameters:
    cluster_dfs (dict): Dictionary of clusters and their data.
    suffix (str): The suffix that determines which clusters to remove (e.g., 'NA', '8.0.2').

    Returns:
    dict: Updated dictionary with specified clusters removed.
    """
    print(f"\nRemoving {suffix} clusters")
    clusters_to_delete = [cluster for cluster in cluster_dfs if cluster.endswith(suffix)]

    for cluster in clusters_to_delete: # Searches the specific cluster in the dic
        del cluster_dfs[cluster]

    return cluster_dfs


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

def create_dotplots_with_thresholds(adata, genes, thresholds, cluster_order, output_dir="canonical/canonical_immune/updated_pts2"):
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
        
    # Ensure leiden_fusion is categorical and reorder it
    #adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype('category')
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.reorder_categories(cluster_order, ordered=True)

    for threshold in thresholds:
        print(f"\nProcessing pts threshold: {threshold}")

        # Extract DGE data
        gene_names, pts = extract_dge_data(adata)

        # Create cluster DataFrames with the current threshold
        cluster_dfs = create_cluster_dfs(gene_names, pts, pts_threshold=threshold)

        # Remove NA clusters
        cluster_dfs = remove_clusters_by_suffix(cluster_dfs, "NA")

        # Compare canonical genes with cluster-specific genes
        filtered_genes = compare_canonical(genes, cluster_dfs)

        # Aggregate filtered genes by gene group
        top_genes_names = top_gene_names(filtered_genes, genes)

        # Example user-defined gene group order
        user_gene_group_order = ["Homeostatic", "MHCII", "Interferon", "DAM", "PVM", "Proliferative"]

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

        # Save dotplots with appropriate filenames
        output_scaled_no_dendro = os.path.join(output_dir, f"dotplot_scaled_no_dendro_{threshold}.png")
        

        dotplot_scaled_no_dendro.savefig(output_scaled_no_dendro, bbox_inches="tight")
        

        plt.close()
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

def check_cluster_order(adata, cluster_order):
    """
    Check for mismatches between existing clusters and the user-defined cluster order.

    Parameters:
    adata (AnnData): The AnnData object.
    cluster_order (list): The user-defined list of clusters to reorder.

    Returns:
    None
    """
    # Extract existing categories in 'leiden_fusion'
    existing_categories = list(adata.obs['leiden_fusion'].cat.categories)

    print("\n--- Existing Categories in 'leiden_fusion' ---")
    print(existing_categories)

    print("\n--- User-Defined Cluster Order ---")
    print(cluster_order)

    # Check for clusters in custom order that don't exist in the data
    missing_in_data = [cluster for cluster in cluster_order if cluster not in existing_categories]
    
    # Check for clusters in the data that aren't in the custom order
    missing_in_order = [cluster for cluster in existing_categories if cluster not in cluster_order]

    if missing_in_data:
        print("\n Clusters in custom order but NOT in 'leiden_fusion':")
        print(missing_in_data)
    
    if missing_in_order:
        print("\n Clusters in 'leiden_fusion' but NOT in custom order:")
        print(missing_in_order)

    if not missing_in_data and not missing_in_order:
        print("\n All categories match! Reordering should work.")





# Main execution block
if __name__ == "__main__":
    # Load data
    adataimm = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad")
    adatamev = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")
    adataneu = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Neu_CentralCanal_raw_norm_ranked_copy_copy.h5ad")

    filtered_adataimm = remove_NA_cat(adataimm)
    filtered_adatamev = remove_NA_cat(adatamev)
    filtered_adataneu = remove_NA_cat(adataneu)

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

    
    recgenesimm15all = load_canonical_from_dir(can_dir_imm_rec_15_all)
    recgenesimm15sep = load_canonical_from_dir(can_dir_imm_rec_15_seperate)
    recgenesimm60all = load_canonical_from_dir(can_dir_imm_rec_60_all)
    recgenesimm60spe = load_canonical_from_dir(can_dir_imm_rec_60_seperate)

    sengenesimm15all = load_canonical_from_dir(can_dir_imm_send_15_all)
    sengenesimm15sep = load_canonical_from_dir(can_dir_imm_send_15_seperate)
    sengenesimm60all = load_canonical_from_dir(can_dir_imm_send_60_all)
    sengenesimm60spe = load_canonical_from_dir(can_dir_imm_send_60_seperate)


    recgenesmev15all = load_canonical_from_dir(can_dir_mev_rec_15_all)
    recgenesmev15sep = load_canonical_from_dir(can_dir_mev_rec_15_seperate)
    recgenesmev60all = load_canonical_from_dir(can_dir_mev_rec_60_all)
    recgenesmev60spe = load_canonical_from_dir(can_dir_mev_rec_60_seperate)

    sengenesmev15all = load_canonical_from_dir(can_dir_mev_send_15_all)
    sengenesmev15sep = load_canonical_from_dir(can_dir_mev_send_15_seperate)
    sengenesmev60all = load_canonical_from_dir(can_dir_mev_send_60_all)
    sengenesmev60spe = load_canonical_from_dir(can_dir_mev_send_60_seperate)


    recgenesneu15all = load_canonical_from_dir(can_dir_neu_rec_15_all)
    recgenesneu15sep = load_canonical_from_dir(can_dir_neu_rec_15_seperate)
    recgenesneu60all = load_canonical_from_dir(can_dir_neu_rec_60_all)
    recgenesneu60spe = load_canonical_from_dir(can_dir_neu_rec_60_seperate)

    sengenesneu15all = load_canonical_from_dir(can_dir_neu_send_15_all)
    sengenesneu15sep = load_canonical_from_dir(can_dir_neu_send_15_seperate)
    sengenesneu60all = load_canonical_from_dir(can_dir_neu_send_60_all)
    sengenesneu60spe = load_canonical_from_dir(can_dir_neu_send_60_seperate)



    # # Create dendogram ot the top genes
    # dendogram_sc(filtered_adata)

    # Define thresholds
    pts_thresholds = [0, 0.2, 0.3]

####

    rec_custom_cluster_order_imm_15 = [
    "Imm.M0Like.1", 
    "Imm.Interferon.0",
    "Imm.PVM.0"]

    send_custom_cluster_order_imm_15 = [ 
    "Imm.Interferon.0", "Imm.DAM.0", 
    "Imm.PVM.0"]

    rec_custom_cluster_order_imm_60 = [
    "Imm.M0Like.1", 
    "Imm.Interferon.0", 
    "Imm.DAM.1"]

    send_custom_cluster_order_imm_60 = [
    "Imm.M0Like.1", 
    "Imm.Interferon.0", "Imm.DAM.0", 
    "Imm.DAM.1"]

####

    rec_custom_cluster_order_mev_15 = [
    "MeV.Endothelial.1", 
    "MeV.Endothelial.2"]

    send_custom_cluster_order_mev_15 = [
    "MeV.Endothelial.2", "MeV.FibCollagen.1", 
    "MeV.FibCollagen.2", "MeV.FibCollagen.3"]

    rec_custom_cluster_order_mev_60 = []

    send_custom_cluster_order_mev_60 = [
    "MeV.FibCollagen.1", 
    "MeV.FibCollagen.3"]

####

    rec_custom_cluster_order_neu_15 = [
    "Neu.Epend.0"]

    rec_custom_cluster_order_neu_60 = [
    "Neu.Epend.0"]

#### Customs orders

    send_collagen_order_15 = [
        "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"
    ]

    send_collagen_order_60 = [
        "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"
    ]


    #### Name of the list of genes txts

    immune_genes_rec_15_all = ["rec|Interferon", "rec|M0.1", "rec|PVM"]

    immune_genes_rec_15_seperate = ["FibColl.2|Interferon", "FibColl.2|M0.1", "FibColl.3|PVM"]
    
    immune_genes_rec_60_all = ["rec|DAM.1", "rec|Intreferon", "rec|M0.1"]

    immune_genes_rec_60_seperate = ["FibColl.3|DAM.1", "FibColl.3|Intreferon", "FibColl.3|M0.1"]
##
    immune_genes_send_15_all = ["DAM.0|send", "Intreferon|send", "PVM|send"]

    immune_genes_send_15_seperate = ["DAM.0|Endo.2", "DAM.0|Epend", "Intreferon|Endo.1", "PVM|Epend"]
    
    immune_genes_send_60_all = ["DAM|send", "Intreferon|send", "M0.1|send"]

    immune_genes_send_60_seperate = ["DAM.0|Epend", "DAM.1|Epend", "Intreferon|Epend", "M0.1|Epend"]
###
    meningeal_genes_rec_15_all = ["rec|Endo.1", "rec|Endo.2"]

    meningeal_genes_rec_15_seperate = ["DAM.0|Endo.2", "Interferon|Endo.1"]
    
    meningeal_genes_rec_60_all = []

    meningeal_genes_rec_60_seperate = []
##
    meningeal_genes_send_15_all = ["Endo.2|send", "FibColl.1|send", "FibColl.2|send", "FibColl.3|send"]

    meningeal_genes_send_15_seperate = ["FibColl.3|send", "FibColl.1|Endo.0", "FibColl.1|Pericytes", "FibColl.2|Interferon", "FibColl.2|M0.1", "FibColl.3|Endo.2", "FibColl.3|PVM"]
    
    meningeal_genes_send_60_all = ["FibColl.1|send", "FibColl.3|send"]

    meningeal_genes_send_60_seperate = ["FibColl.1|Pericytes", "FibColl.3|DAM.1", "FibColl.3|Endo.0", "FibColl.3|Epend", "FibColl.3|Interferon", "FibColl.3|M0.1"]
###
    neu_genes_rec_15_all = ["rec|Epend"]

    neu_genes_rec_15_seperate = ["DAM.0|Epend", "PVM|Epend"]
    
    neu_genes_rec_60_all = ["rec|Epend"]

    neu_genes_rec_60_seperate = ["DAM.0|Epend", "DAM.1|Epend", "FibColl.3|Epend", "Interferon|Epend", "M0.1|Epend"]
##
    neu_genes_send_15_all = []

    neu_genes_send_15_seperate = []
    
    neu_genes_send_60_all = []

    neu_genes_send_60_seperate = []

    # Check for mismatches before reordering
    check_cluster_order(filtered_adata, custom_cluster_order)
    

    # Generate dotplots for each threshold
    create_dotplots_with_thresholds(filtered_adata, genes, pts_thresholds, custom_cluster_order)