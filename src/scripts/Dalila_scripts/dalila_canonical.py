# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# This script will display the gene expression values of Dalila's canonical genes from Meningeal dataset


# Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.ss
    """
    print("Loading h5ad file")
    return sc.read_h5ad(file_path)


def split_adata_by_injury_day(adata: sc.AnnData):
    """
    Split an AnnData object into four sub-objects based on the 'injury_day' column in obs.

    Returns:
        tuple: (adata_uninjured_0, adata_sham_15, adata_injured_15, adata_injured_60)
    """
    print("Splitting AnnData by injury_day...")
    adata_control = adata[adata.obs['injury_day'].isin(["uninjured_0", "sham_15"])].copy()
    adata_injured_15  = adata[adata.obs['injury_day'] == "injured_15"].copy()
    adata_injured_60  = adata[adata.obs['injury_day'] == "injured_60"].copy()
    adata_injured = adata[adata.obs['injury_day'].isin(["injured_15", "injured_60"])].copy()

    return adata_control, adata_injured_15, adata_injured_60, adata_injured


def remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    
    mask_NA = adata.obs['leiden_fusion'] != 'MeV.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2


def remove_clusters(adata: sc.AnnData, clusters_to_remove: list, cluster_key: str = 'leiden_fusion'):
    """
    Remove specified clusters from the AnnData object.

    Parameters:
    adata (sc.AnnData): The annotated data matrix.
    clusters_to_remove (list): List of cluster names to remove.
    cluster_key (str): The key in adata.obs that contains cluster annotations. Default is 'leiden_fusion'.

    Returns:
    sc.AnnData: A new AnnData object with the specified clusters removed.
    """
    print(f"Removing clusters: {clusters_to_remove}")

    # Create a mask to filter out the specified clusters
    mask = ~adata.obs[cluster_key].isin(clusters_to_remove)
    
    # Apply the mask to create a new AnnData object
    adata_filtered = adata[mask].copy()
    
    return adata_filtered

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

def create_dotplots_with_thresholds(adata, genes, thresholds, cluster_order, name, prefix ,output_dir="canonical/canonical_dalila/new_sep"):
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


    adata = imm_keep_only_selected_clusters(adata=adata, clusters_to_keep=cluster_order)

    #print(adata['leiden_fusion'].cat.categories.to_list())
    
    # Convert to categorical and reorder only existing categories
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype('category')

    # # Find actual categories present in this subset
    # present_categories = [cat for cat in cluster_order if cat in adata.obs['leiden_fusion'].cat.categories]

    # # Inform if some clusters are missing
    # missing_clusters = [cat for cat in cluster_order if cat not in present_categories]
    # if missing_clusters:
    #     print(f"Skipping missing clusters in dotplot: {missing_clusters}")

    # Reorder with only present ones
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.reorder_categories(cluster_order, ordered=True)


    for threshold in thresholds:
        print(f"\nProcessing pts threshold: {threshold}")

        # Extract DGE data
        gene_names, logfoldchanges, pvals_adj, pts = extract_dge_data(adata)

        # Create cluster DataFrames with the current threshold
        cluster_dfs = create_cluster_dfs(gene_names, logfoldchanges, pvals_adj, pts, pts_threshold=threshold)

        # # Remove NA clusters
        # cluster_dfs = remove_clusters_by_suffix(cluster_dfs, "NA")

        # Compare canonical genes with cluster-specific genes
        filtered_genes = compare_canonical(genes, cluster_dfs)

        # # Export filtered genes to Excel
        # export_to_excel(filtered_genes, pts_threshold=threshold, prefix=prefix)

        # Aggregate filtered genes by gene group
        top_genes_names = top_gene_names(filtered_genes, genes)

        # Example user-defined gene group order
        user_gene_group_order = ["Mylip_Idol", "Mesenchymal_state"]

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
        output_scaled_no_dendro = os.path.join(output_dir, f"{prefix}_dotplot_scaled_{threshold}_{name}.pdf")
        #output_normal_no_dendro = os.path.join(output_dir, f"{prefix}_dotplot_normal_{threshold}_{name}.png")


        dotplot_scaled_no_dendro.savefig(output_scaled_no_dendro, bbox_inches="tight")
        #dotplot_normal_no_dendro.savefig(output_normal_no_dendro, bbox_inches="tight")

        plt.close()
        print(f"Saved dotplots for threshold {threshold}:")
        print(f"  - {output_scaled_no_dendro}")
        #print(f"  - {output_normal_no_dendro}")
        


def extract_dge_data(adata):
    dge_fusion = adata.uns.get('rank_genes_groups_leiden_fusion', None)
    
    gene_names = pd.DataFrame(dge_fusion.get('names', []))
    logfoldchanges = pd.DataFrame(dge_fusion.get('logfoldchanges', []))
    pvals_adj = pd.DataFrame(dge_fusion.get('pvals_adj', []))
    pts = pd.DataFrame(dge_fusion.get('pts', []))

    return gene_names, logfoldchanges, pvals_adj, pts



# Step 3: Create cluster dataframes with filtered data
def create_cluster_dfs(gene_names, logfoldchanges, pvals_adj, pts, pts_threshold=0):
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
            'logfoldchange': logfoldchanges.iloc[:, i].values,
            'pval_adj': pvals_adj.iloc[:, i].values,
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

def export_to_excel(filtered_genes, pts_threshold, prefix, output_dir="excels/dalila"):
    """
    Export filtered genes into separate Excel files, one for each gene group and condition.
    Each sheet in the Excel file represents a cluster, containing gene expression data for that gene group.

    Parameters:
    - filtered_genes (dict): Dictionary of filtered genes grouped by gene group and cluster.
    - pts_threshold (float): The threshold used for filtering.
    - prefix (str): Prefix to indicate injury condition.
    - output_dir (str): Directory where the Excel files will be saved.

    Returns:
    - None
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    for gene_group, clusters in filtered_genes.items():
        output_file = os.path.join(output_dir, f"{prefix}_{gene_group}_{pts_threshold}.xlsx")
        print(f"\n Exporting data for gene group '{gene_group}' ({prefix}) to: {output_file}")

        with pd.ExcelWriter(output_file) as writer:
            empty_sheets = True  # Flag to track if we have at least one sheet
            
            for cluster_name, df in clusters.items():
                if not df.empty:
                    df = df.sort_index()  # Sort genes alphabetically
                    df.to_excel(writer, sheet_name=cluster_name[:31])  # Sheet names max 31 chars
                    empty_sheets = False
                else:
                    print(f" WARNING: No data for {cluster_name}. Skipping sheet.")
            
            if empty_sheets:
                placeholder_df = pd.DataFrame({"Message": ["No data available"]})
                placeholder_df.to_excel(writer, sheet_name="Placeholder")

        print(f" Saved Excel file: {output_file}")



def filter_cells_by_gene_expression(adata: sc.AnnData, gene_name: str):
    """
    Filter the AnnData object to retain only cells that express a specific gene (expression > 0).

    Parameters:
    adata (AnnData): The AnnData object containing gene expression data.
    gene_name (str): The gene to filter on.

    Returns:
    AnnData: A new AnnData object containing only cells where the specified gene is expressed.
    """
    if gene_name not in adata.var_names:
        print(f" WARNING: Gene '{gene_name}' not found in the dataset. Returning original dataset.")
        return adata

    print(f"Filtering cells that express the gene '{gene_name}'...")

    # Create mask: Select cells where expression of the gene is greater than 0
    mask = adata[:, gene_name].X > 0

    # Apply mask to create new filtered dataset
    filtered_adata = adata[mask].copy()

    print(f"Filtered dataset contains {filtered_adata.n_obs} cells (out of {adata.n_obs}) that express '{gene_name}'.")

    return filtered_adata



if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")

    # Remove NA categories
    filtered_adata = remove_NA_cat(adata)

    # Split into injury_day groups
    adata_control, adata_injured_15, adata_injured_60, adata_injured = split_adata_by_injury_day(filtered_adata)
    print(adata_injured_15.obs["leiden_fusion"])

    # Group into a dictionary for looping
    adata_groups = {
        "control": adata_control,
        "injured_15": adata_injured_15,
        "injured_60": adata_injured_60,
        "injured": adata_injured
        # "fulldays": filtered_adata
    }

    # Load canonical gene lists from a directory
    canonical_genes_dir = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Dalila"
    genes = load_canonical_from_dir(canonical_genes_dir)

    # Thresholds and cluster order
    pts_thresholds = [0]
    custom_cluster_order_all = ["MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Epithelial.0",
                            "MeV.SMC.0", "MeV.Pericytes.0", "MeV.VLMC.0", "MeV.VLMC.1" , "MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3",
                            "MeV.FibLaminin.0", "MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2", "MeV.Fib.5", "MeV.Fib.3", "MeV.Fib.4", "MeV.FibProlif.0"]

    custom_cluster_order_mylip = ["MeV.Pericytes.0", "MeV.FibCollagen.1", "MeV.FibProlif.0"]


    cond = "all"
    cond1= "control"
    cond2= "injured15"
    cond3= "injured60"
    cond4= "injured"

    
    # Generate dotplots
    create_dotplots_with_thresholds(filtered_adata, genes, pts_thresholds, custom_cluster_order_mylip, None ,prefix=cond)

    create_dotplots_with_thresholds(adata_control, genes, pts_thresholds, custom_cluster_order_mylip, None,prefix=cond1)

    create_dotplots_with_thresholds(adata_injured_15, genes, pts_thresholds, custom_cluster_order_mylip, None,prefix=cond2)

    create_dotplots_with_thresholds(adata_injured_60, genes, pts_thresholds, custom_cluster_order_mylip, None,prefix=cond3)

    create_dotplots_with_thresholds(adata_injured, genes, pts_thresholds, custom_cluster_order_mylip, None,prefix=cond4)


    # # Process each subset
    # for label, ad in adata_groups.items():
    #     print(f"\n--- Processing condition: {label} ---")

    #     # Filter cells by gene expression
    #     gene_filtered_adata = filter_cells_by_gene_expression(ad, "Mylip")

    #     # Remove unwanted clusters
    #     clusters_to_remove = ['MeV.ImmuneDoublets.0', 'MeV.LowQuality.0', "MeV.FibUnknown.6", "MeV.EndoUnknow.4"]
    #     adata_filtered = remove_clusters(ad, clusters_to_remove)


    #     # Check cluster order
    #     check_cluster_order(adata_filtered, custom_cluster_order_mylip)

    #     # Generate dotplots
    #     create_dotplots_with_thresholds(adata_filtered, genes, pts_thresholds, custom_cluster_order_mylip, name1, prefix=label)

    #     create_dotplots_with_thresholds(adata_control, genes, pts_thresholds, custom_cluster_order_mylip, name1, prefix=cond1)

    #     create_dotplots_with_thresholds(adata_injured_15, genes, pts_thresholds, custom_cluster_order_mylip, name1, prefix=cond2)

    #     create_dotplots_with_thresholds(adata_injured_60, genes, pts_thresholds, custom_cluster_order_mylip, name1, prefix=cond3)

    #     create_dotplots_with_thresholds(adata_injured, genes, pts_thresholds, custom_cluster_order_mylip, name1, prefix=cond4)
