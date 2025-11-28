# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    return sc.read_h5ad(file_path)


def umap_reso_cluster(adata, resolution_name, output_dir="reso/reso_meninge/updates"):
    """
    Plot UMAP for a specific resolution with cluster numbers displayed at the centroid of each cluster.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    resolution_name (str): The resolution name to be used for coloring the UMAP plot (e.g., 'leiden_fusion).

    Returns:
    None
    """

    # Plot UMAP for the specified resolution
    ax = sc.pl.umap(adata, color=resolution_name, title=f"UMAP - {resolution_name}", return_fig=True)
    
    ax_ondata = sc.pl.umap(adata, 
                           color=resolution_name, 
                           title=f"UMAP - Immune Clusters",
                           legend_loc = 'on data',
                           legend_fontsize=10,  
                           legend_fontweight = 'bold',  
                           legend_fontoutline = 1,
                           size=50,
                           return_fig=True)

    
    # Save the UMAP plot as an image (optional)
    output_path = os.path.join(output_dir, f"umap_Immune_{resolution_name}_n.pdf") 
    output_path_leg = os.path.join(output_dir, f"umap_Immune_{resolution_name}_l.pdf")
    ax.figure.savefig(output_path, bbox_inches="tight")
    ax_ondata.figure.savefig(output_path_leg, bbox_inches="tight")
    # print(f"UMAP plot saved as {output_path}")
    plt.close()  # Close the plot to avoid overlap



# Step 2: Extract Differentially Expressed Genes (DEGs)
def extract_dge_data(adata):
    """
    Extract the different type of data from the AnnData object like the differential gene expression values from each cell.

    Parameters:
    adata (AnnData): The AnnData object containing the DGE data.

    Returns:
    tuple: DataFrames for gene names, logfoldchanges, adjusted p-values, scores, and pts.
    """
    dge_fusion = adata.uns['rank_genes_groups_leiden_fusion']
    
    # Convert the extracted data into DataFrames
    gene_names = pd.DataFrame(dge_fusion['names'])
    logfoldchanges = pd.DataFrame(dge_fusion['logfoldchanges'])
    pvals_adj = pd.DataFrame(dge_fusion['pvals_adj'])
    scores = pd.DataFrame(dge_fusion['scores'])
    pts = pd.DataFrame(dge_fusion['pts'])
    
    return gene_names, logfoldchanges, pvals_adj, scores, pts


# Step 3: Create cluster dataframes with filtered data
def create_cluster_dfs(gene_names, logfoldchanges, pvals_adj, scores, pts, sort_by_logfc=False, pts_threshold=0):
    """
    Create a dictionary of dataframes per cluster, with the respective gene expression data.
    By default this function will not order the genes by fold change and the default minimun pts is 0

    Parameters:
    gene_names (DataFrame): DataFrame with gene names.
    logfoldchanges (DataFrame): DataFrame with logfoldchange values for each gene from each cell.
    pvals_adj (DataFrame): DataFrame with adjusted p-values for each gene from each cell.
    scores (DataFrame): DataFrame with willcoxon scores for each gene from each cell.
    pts (DataFrame): DataFrame with p values for each gene from each cell.
    sort_by_logfc (bool): If True, sort the filtered DataFrame by log fold change in descending order.
    pts_threshold (float): The minimum percentage of expression value to filter genes.

    logFoldChange: Magnitude of the difference in gene expression between the compared groups.
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
            'score': scores.iloc[:, i].values,
            'logfoldchange': logfoldchanges.iloc[:, i].values,
            'pvals_adj': pvals_adj.iloc[:, i].values,
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

        # Sort by log fold change if sort_by_logfc is True
        if sort_by_logfc:
            filtered_df = filtered_df.sort_values(by='logfoldchange', ascending=False)

        # Assigns the cluster name as the name of the resolution atributed
        cluster_names = gene_names.columns[i]
        cluster_dfs[cluster_names] = filtered_df  # Store the respective clusters in the dictionary

    return cluster_dfs


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
    clusters_to_delete = [cluster for cluster in cluster_dfs if cluster.endswith(suffix)]

    for cluster in clusters_to_delete: # Searches the specific cluster in the dic
        del cluster_dfs[cluster]

    return cluster_dfs


# Step 5: Filter, sort, and select top 5 genes per cluster

def filter_pval(df):
    """
    First step of 6 in order to select the top genes:
    Filter genes by adj p-value and select the top 5 statisticaly significant genes.

    Parameters:
    df (DataFrame): DataFrame containing genes for a cluster.

    Returns:
    DataFrame: Filtered DataFrame with top 5 statisticaly significant genes .
    """
    mask_pval = df['pvals_adj'] < 0.05
    filtered_df = df[mask_pval]
    return filtered_df.head(5)


def has_neg_lfc(top_genes):
    """
    Step 2 of 6 to determinate if the cluster have any negative log fold change genes (under expressed).

    Parameters:
    top_genes (DataFrame): DataFrame of the top 5 genes for a cluster.

    Returns:
    bool: True if negative log fold change exists, False otherwise.
    """
    return (top_genes['logfoldchange'] < 0).any()


def count_neg_lfc(top_genes):
    """
    Step 3 of 6 to check how many of negative log fold change genes exist (under expressed).

    Parameters:
    top_genes (DataFrame): DataFrame of the top 5 genes for a cluster.

    Returns:
    int: Number of genes with negative log fold changes.
    """
    return (top_genes['logfoldchange'] < 0).sum()


def replace_neg_lfc(df, top_genes, num_nega_lfc):
    """
    Step 4 of 6 that replaces the negative log fold change genes with the positive ones from the unfiltered df

    Parameters:
    df (DataFrame): Original DataFrame of the cluster.
    top_genes (DataFrame): DataFrame of the top 5 genes for the cluster.
    num_nega_lfc (int): Number of genes with negative log fold changes.

    Returns:
    DataFrame: DataFrame with the top 5 genes.
    """
    # Replace the selected genes with the genes from the original df
    top_genes_positive = top_genes[top_genes['logfoldchange'] > 0]

    replacements = []
    for gene, values in df.iterrows():
        if gene not in top_genes_positive.index:
            replacements.append(values) # Prevent of appending duplicated genes
        if len(replacements) == num_nega_lfc:
            break  # Stop once enough replacement genes are found

    replacements_df = pd.DataFrame(replacements)
    return pd.concat([top_genes_positive, replacements_df]) # concact the 2 lists


def replace_all(df):
    """
    Step 5 of 6 where the program selects the top 5 genes as a fallback if no filtered genes are found.

    Parameters:
    df (DataFrame): Original DataFrame of the cluster.

    Returns:
    DataFrame: DataFrame with the top 5 fallback genes.
    """
    return df.head(5)


def select_top_genes(cluster_dfs):
    """
    Final step of the selection of apropriated genes 6 of 6:
    Main function to select the top 5 genes per cluster, replacing genes with negative log fold change.

    Parameters:
    cluster_dfs (dict): Dictionary containing dataframes of each cluster.

    Returns:
    dict: Dictionary containing the top 5 genes per cluster.
    """
    top_genes_cluster = {}
    
    for cluster, df in cluster_dfs.items():
        top_genes = filter_pval(df) # First remove the not significant genes 

        if has_neg_lfc(top_genes): # If there is any neg lfc genes
            num_nega_lfc = count_neg_lfc(top_genes)
            top_genes_cluster[cluster] = replace_neg_lfc(df, top_genes, num_nega_lfc) # replace the neg genes
        elif top_genes.empty:
            top_genes_cluster[cluster] = replace_all(df) # replace all if it is empty
        else:
            top_genes_cluster[cluster] = top_genes # if the cluster doesnt have neg genes just continues

    return top_genes_cluster


# Step 6: Collect top gene names for each cluster
def top_gene_names(top_genes_cluster):
    """
    Create a dictionary of top gene names for each cluster.

    Parameters:
    top_genes_cluster (dict): Dictionary containing top genes for each cluster.

    Returns:
    dict: Dictionary of top gene names per cluster.
    """
    top_genes_names = {}

    for cluster, df in top_genes_cluster.items():
        top_genes_names[cluster] = df.index.tolist()

    return top_genes_names


def create_dotplots_with_thresholds(adata, thresholds, clusters_to_remove, cluster_order, output_dir="dotplots/meningeal/leiden_fusion"):
    """
    Create and save dotplots for different pts thresholds, with and without dendrograms.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    thresholds (list of float): List of pts thresholds to generate dotplots for.
    output_dir (str): Directory to save the dotplots.

    Returns:
    None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype('category')
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.reorder_categories(cluster_order, ordered=True)


    for threshold in thresholds:
        print(f"\n🔹 Processing dotplots for pts threshold: {threshold}")

        # Extract DGE data dynamically for each threshold
        print("   - Extracting differential gene expression data...")
        gene_names, logfoldchanges, pvals_adj, scores, pts = extract_dge_data(adata)

        # Create cluster DataFrames with current threshold
        print(f"   - Applying pts threshold: {threshold}")
        cluster_dfs = create_cluster_dfs(
            gene_names, logfoldchanges, pvals_adj, scores, pts, 
            sort_by_logfc=True, pts_threshold=threshold
        )

        # Remove NA clusters
        cluster_dfs = remove_clusters_by_suffix(cluster_dfs, "NA")

        filtered_clusters = remove_clusters2(cluster_dfs, clusters_to_remove)

        ordered_clusters = order_clusters(filtered_clusters, cluster_order)

        export_to_excel_all(ordered_clusters, threshold,string = "all")

        # Create dendrogram (if not already present)
        print("   - Checking and creating dendrogram if necessary...")
        # Preform Dendrogram
        #dendogram_sc(adata)
        #plot_dendrogram(adata)

        # Select the top genes for each cluster
        print("   - Selecting top genes for each cluster...")
        top_genes_cluster = select_top_genes(cluster_dfs)

        # Add an asterisk to clusters with non-significant genes
        top_genes_cluster = addasterix(top_genes_cluster)

        filtered_clustersss = remove_clusters2(top_genes_cluster, clusters_to_remove)

        ordered_clusterss = order_clusters(filtered_clustersss, cluster_order)

        print("printaaa")
        print(filtered_clustersss)

        export_to_excel_all(ordered_clusterss, threshold, string = "top")

        # Collect top gene names for visualization
        top_genes_names = top_gene_names(top_genes_cluster)

        top_genes_names = filter_and_order_clusters(top_genes_names, clusters_to_remove, cluster_order)


        # Reorder clusters (ON and OFF for dendrogram)
        # print("   - Reordering clusters based on dendrogram...")
        # ordered_genes_dendro = reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram=True)
        # ordered_genes_no_dendro = reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram=False)

        print("   - Generating dotplots...")

        # (1) With Dendrogram
        dotplot_normal = sc.pl.rank_genes_groups_dotplot(
            adata,
            var_names=top_genes_names,
            groupby='leiden_fusion',
            key='rank_genes_groups_leiden_fusion',
            cmap='bwr',
            vmin=-4,
            vmax=4,
            values_to_plot='logfoldchanges',
            colorbar_title='Magnitude of expression',
            use_raw=False,
            dendrogram=False,
            return_fig=True
        )

        # # (2) scaled
        # dotplot_scaled = sc.pl.rank_genes_groups_dotplot(
        #     adata,
        #     var_names=top_genes_names,
        #     groupby='leiden_fusion',
        #     key='rank_genes_groups_leiden_fusion',
        #     cmap='Greys',
        #     use_raw=False,
        #     standard_scale='var',
        #     dendrogram=False,
        #     return_fig=True
        # )

        # Save plots
        #output_scale = os.path.join(output_dir, f"dotplot_scale_{threshold}.pdf")
        output_normal = os.path.join(output_dir, f"dotplot_normal_{threshold}.pdf")

        #dotplot_scaled.savefig(output_scale, bbox_inches="tight")
        dotplot_normal.savefig(output_normal, bbox_inches="tight")

        plt.close()
        #print(f" Saved: {output_scale}")
        print(f" Saved: {output_normal}")



# Print of gene names and clusters
def print_gene_names(top_genes_names):

    for cluster, df in top_genes_names.items():
        print(f"\nTop genes in {cluster}:")
        print(df)

def print_clusters(top_genes_cluster):

    for cluster, df in top_genes_cluster.items():
        print(f"\nTop genes in {cluster}:")
        print(df)

# Remove NA
def remove_NA_cat(adata: sc.AnnData):
    
    mask_NA = adata.obs['leiden_fusion'] != 'MeV.NA' #creates mask for remove NA cells
    #print(mask_NA)    
    adata2 = adata[mask_NA] #apply mask

    # print(adata2.obs['leiden_fusion'].cat.categories.to_list())
    # adata2.obs['leiden_fusion'] = adata2.obs['leiden_fusion'].cat.remove_unused_categories()
    # print(adata2.obs['leiden_fusion'].cat.categories.to_list())
    # print(adata.obs['leiden_fusion'])

    return adata2

# add an astrik to non sign genes 
def nonsigngene(top_genes):
    
    return (top_genes['pvals_adj'] > 0.05 ).any()

def addasterix(top_genes_cluster):
    
    updated_cluster = {}

    for cluster, df in top_genes_cluster.items():

        # Check if the cluster has any non-significant genes
        if nonsigngene(df):
            # Add an asterisk to the cluster name
            updated_cluster[cluster + '*'] = df
        else:
            updated_cluster[cluster]  = df

    return updated_cluster



def dendogram_sc(adata):
    """
    Compute and save the dendrogram for a specific group and store it in a specified key.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    dendrogram_key (str): The key to save the dendrogram in `adata.uns`.

    Returns:
    None
    """

    # Compute the dendrogram
    print(f"Computing dendrogram for leiden_fusion...")
    sc.tl.dendrogram(
        adata,
        groupby='leiden_fusion',
        use_rep='X_pca',
        cor_method='spearman',
        linkage_method='ward',
        use_raw=False
    )

def plot_dendrogram(adata, output_dir="dendrogram/dendrogram_meningeal"):
    """
    
    """

    # Create the directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Plot the dendrogram
    print(f"Plotting dendrogram for leiden_fusion...")
    sc.pl.dendrogram(
        adata,
        groupby='leiden_fusion',
        dendrogram_key=None,
        orientation='top',
        show=False
    )

    # Save the plot
    output_path = os.path.join(output_dir, f"leiden_fusion_dendrogram.png")
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap

def reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram, dendrogram_key = 'dendrogram_leiden_fusion'):
    """
    Reorder clusters based on dendrogram order if reorder is True.

    Parameters:
    adata (AnnData): The AnnData object containing the dendrogram information.
    top_genes_names (dict): Dictionary of top genes per cluster.
    reorder (bool): Flag indicating whether to reorder clusters based on dendrogram.
    dendrogram_key (str): The key in `adata.uns` containing the dendrogram information.

    Returns:
    dict: Reordered dictionary of top genes per cluster, or the original dictionary if reorder is False.
    """

    if not dendrogram:
        print("Reordering is disabled. Returning original cluster order.")
        return top_genes_names

    # Extract dendrogram order
    print("Reordering clusters based on dendrogram...")
    ivl = adata.uns[dendrogram_key]['dendrogram_info']['ivl']
    print(f"Dendrogram order: {ivl}")

    # Prepare a reordered dictionary
    reordered_dict = {}

    # Iterate through the dendrogram order and match clusters (ignoring '*')
    for cluster in ivl:
        # Handle cluster names with and without '*'
        matched_cluster = next(
            (key for key in top_genes_names if key.rstrip('*') == cluster), None
        )
        if matched_cluster:
            reordered_dict[matched_cluster] = top_genes_names[matched_cluster]

    return reordered_dict

def export_top_genes_to_txt(top_genes_cluster, threshold, output_dir="excels/meningeal/updates"):
    """
    Export the top genes from each cluster to a single text file.

    Parameters:
    top_genes_cluster (dict): Dictionary containing top genes for each cluster.
    threshold (float): The threshold value used for filtering, included in the filename.
    output_dir (str): Directory to save the text file.

    Returns:
    None
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Define the output file path
    output_file = os.path.join(output_dir, f"top_genes_clusters_{threshold}.txt")

    print(f"\nExporting top genes to text file: {output_file}")

    with open(output_file, 'w') as f:
        for cluster, df in top_genes_cluster.items():
            f.write(f"Cluster: {cluster}\n")  # Write cluster name
            f.write(f"Top Genes:\n")
            for gene in df.index.tolist():
                f.write(f"  - {gene}\n")  # Write each gene name
            f.write("\n" + "-"*40 + "\n")  # Separator between clusters

    print(f"Text file saved: {output_file}")


def export_to_excel_all(cluster_dfs, threshold,string, output_dir="excels/meningeal/new_tese/dge"):
    """
    Export all differential expression data per cluster to an Excel file,
    including logfoldchanges, pvals_adj, scores, and pts.

    Parameters:
    cluster_dfs (dict): Dictionary where keys = cluster names and values = DataFrames
                        containing ['logfoldchanges', 'pvals_adj', 'scores', 'pts'].
    threshold (float): The threshold value used for filtering.
    output_dir (str): Directory where the Excel file will be saved.

    Returns:
    None
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define file path dynamically using the threshold
    output_file = os.path.join(output_dir, f"{string}_genes_clusters_{threshold}.xlsx")
    print(f"\n🟩 Exporting full DGE results to Excel: {output_file}")

    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for cluster, df in cluster_dfs.items():
                # Ensure DataFrame has all expected columns
                expected_cols = ["logfoldchange", "pvals_adj", "score", "pts"]
                missing_cols = [col for col in expected_cols if col not in df.columns]
                if missing_cols:
                    print(f"⚠️  Warning: Missing columns {missing_cols} in cluster {cluster}")

                # Sort columns consistently
                cols = [col for col in expected_cols if col in df.columns]
                df_to_write = df[cols].copy()
                df_to_write.index.name = "Gene"

                # Write to Excel sheet named after the cluster
                sheet_name = str(cluster)[:31]  # Excel sheet name limit
                df_to_write.to_excel(writer, sheet_name=sheet_name)

                print(f"   ✅ Wrote {len(df_to_write)} genes for cluster '{cluster}'")

        print(f"✅ Excel file successfully saved at: {output_file}")

    except Exception as e:
        print(f"❌ Error exporting Excel file: {e}")
        import traceback
        traceback.print_exc()


def order_clusters(cluster_dict, custom_order):
    """
    Reorder clusters in a dictionary according to a custom order.

    Parameters
    ----------
    cluster_dict : dict
        Dictionary of cluster_name → DataFrame
    custom_order : list
        List of cluster names in desired order

    Returns
    -------
    dict :
        Ordered dictionary according to custom_order. Clusters not in custom_order
        are appended at the end in their original order.
    """
    ordered = {k: cluster_dict[k] for k in custom_order if k in cluster_dict}

    # Append clusters not in custom_order at the end
    remaining = [k for k in cluster_dict if k not in custom_order]
    for k in remaining:
        ordered[k] = cluster_dict[k]

    return ordered

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

def filter_and_order_clusters(top_genes_dict, clusters_to_remove, custom_order):
    """
    Remove specific clusters and reorder the remaining keys according to a custom order.

    Parameters:
    top_genes_dict (dict): Dictionary mapping cluster names to gene lists.
    clusters_to_remove (list): List of cluster names (keys) to remove.
    custom_order (list): List defining the desired cluster order.

    Returns:
    dict: A new dictionary filtered and ordered by the custom order.
    """
    print("\n🧹 Filtering and ordering clusters...")

    # 1️⃣ Remove unwanted clusters
    filtered_dict = {k: v for k, v in top_genes_dict.items() if k not in clusters_to_remove}
    
    removed = [c for c in clusters_to_remove if c in top_genes_dict]
    missing = [c for c in clusters_to_remove if c not in top_genes_dict]

    if removed:
        print(f"✅ Removed {len(removed)} clusters: {removed}")
    if missing:
        print(f"⚠️  Skipped {len(missing)} (not found): {missing}")
    
    # 2️⃣ Reorder according to custom_order
    ordered_dict = {k: filtered_dict[k] for k in custom_order if k in filtered_dict}

    # 3️⃣ Add any remaining clusters not in custom_order (at the end)
    remaining = [k for k in filtered_dict.keys() if k not in custom_order]
    if remaining:
        print(f"ℹ️ Added {len(remaining)} clusters not in custom order: {remaining}")
        for k in remaining:
            ordered_dict[k] = filtered_dict[k]

    print(f"📊 Final cluster count: {len(ordered_dict)}")
    return ordered_dict

def remove_clusters2(cluster_dict, clusters_to_remove):
    """
    Remove unwanted clusters (dictionary keys) from a dict of DataFrames.

    Parameters
    ----------
    cluster_dict : dict
        Mapping cluster → DataFrame
    clusters_to_remove : list
        List of keys to remove from the dictionary

    Returns
    -------
    dict
        Filtered dictionary with specified clusters removed
    """

    print("\n🧹 Removing unwanted clusters...")

    filtered = {k: v for k, v in cluster_dict.items() if k not in clusters_to_remove}

    removed = [c for c in clusters_to_remove if c in cluster_dict]
    missing = [c for c in clusters_to_remove if c not in cluster_dict]

    if removed:
        print(f"✅ Removed {len(removed)} clusters: {removed}")
    if missing:
        print(f"⚠️  Skipped {len(missing)} (not found): {missing}")

    print(f"📊 Final cluster count: {len(filtered)}")

    return filtered

def no_filter_dotplot(adata, cluster_order, output_dir="dotplots/meningeal/leiden_fusion"):


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Reorder adata categories
    adata.obs['leiden_fusion_old1'] = adata.obs['leiden_fusion_old1'].cat.reorder_categories(cluster_order, ordered=True)

    # 2. Extract top 5 genes + remove NA clusters
    top_genes_names = get_top_genes_for_dotplot(
        adata,
        key='rank_genes_groups_leiden_fusion_old1',
        n_genes=5,
        cluster_order=cluster_order,
        remove_suffix="NA"
    )

    # 3. Reorder dict (optional)
    top_genes_names = order_clusters(top_genes_names, cluster_order)

    # 4. Plot
    dotplot_normal = sc.pl.rank_genes_groups_dotplot(
        adata,
        var_names=top_genes_names,
        groupby='leiden_fusion_old1',
        key='rank_genes_groups_leiden_fusion_old1',
        cmap='bwr',
        vmin=-4, vmax=4,
        use_raw=False,
        values_to_plot='logfoldchanges',
        dendrogram=False,
        return_fig=True
    )


    # Save plots
    output_normal = os.path.join(output_dir, f"dotplot_dge_no_filter.png")
    dotplot_normal.savefig(output_normal, bbox_inches="tight")

    plt.close()
    print(f" Saved: {output_normal}")


def get_top_genes_for_dotplot(adata, 
                              key='rank_genes_groups_leiden_fusion_old1',
                              n_genes=5,
                              cluster_order=None,
                              remove_suffix="NA"):
    """
    Extract top n_genes used in Scanpy dotplot and apply ordering + NA removal.
    """

    rgg = adata.uns[key]
    all_clusters = list(rgg['names'].dtype.names)

    # If no custom order provided, keep original
    if cluster_order is None:
        cluster_order = all_clusters

    top_dict = {}

    # Extract top n genes in desired order
    for cluster in cluster_order:
        if cluster not in all_clusters:
            continue

        genes = list(rgg['names'][cluster][:n_genes])
        top_dict[cluster] = genes

    # Remove suffix clusters (like NA)
    if remove_suffix:
        top_dict = remove_varname_clusters_by_suffix(top_dict, remove_suffix)

    return top_dict


def order_clusters(cluster_dict, cluster_order):
    """
    Reorder a dict according to a list of cluster_order.
    Keys not in cluster_order are appended at the end.
    """
    ordered = {k: cluster_dict[k] for k in cluster_order if k in cluster_dict}

    # Add remaining clusters at the bottom
    for k in cluster_dict:
        if k not in ordered:
            ordered[k] = cluster_dict[k]

    return ordered

def remove_varname_clusters_by_suffix(varname_dict, suffix="NA"):
    """
    Remove clusters from a var_names dictionary (cluster → list of genes) 
    if their cluster name ends with a specified suffix.

    Parameters:
        varname_dict (dict): {cluster: [gene1, gene2, ...]}
        suffix (str): Suffix to remove (e.g. "NA")

    Returns:
        dict: Filtered dictionary
    """
    print(f"\n🧹 Removing clusters ending with suffix: '{suffix}'")

    filtered = {
        k: v for k, v in varname_dict.items()
        if not k.endswith(suffix)
    }

    removed = [k for k in varname_dict if k.endswith(suffix)]
    if removed:
        print(f"   Removed {len(removed)} clusters: {removed}")
    else:
        print("   No clusters removed.")

    return filtered


# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")


    filtered_adata = remove_NA_cat(adata)

    print(filtered_adata.obs["leiden_fusion_old1"].cat.categories.to_list())

    clusters_to_remove = ['MeV.ImmuneDoublets.0', 'MeV.LowQuality.0', 'MeV.EndoUnknow.4', 'MeV.FibUnknown.6']
    adatas_filtered = remove_clusters(filtered_adata, clusters_to_remove)

    # #Create cluster resolutions UMAP
    # umap_reso_cluster(adatas_filtered, 'leiden_fusion')

    pts_thresholds = [0.3, 0.4, 0.5]

    custom_cluster_order = ["MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Epithelial.0",
                            "MeV.SMC.0", "MeV.Pericytes.0", "MeV.VLMC.0", "MeV.VLMC.1" , "MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3",
                            "MeV.FibLaminin.0", "MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2", "MeV.Fib.5", "MeV.Fib.3", "MeV.Fib.4", "MeV.FibProlif.0"]
    

    custom_cluster_order_no_filter = ["MeV.1.4.1", "MeV.1.4.15", "MeV.1.4.5", "MeV.4.21", "MeV.1.4.20",
                                      "MeV.4.26", "MeV.4.31", "MeV.4.1", 
                                      "MeV.4.4", "MeV.1.4.0", "MeV.4.12" , 
                                      "MeV.2.8", "MeV.1.4.12", "MeV.4.30", "MeV.1.4.4", "MeV.3.17", 
                                      "MeV.2.1", "MeV.1.4.7", "MeV.1.4.6", "MeV.1.4.11", "MeV.4.34", "MeV.1.4.13", "MeV.1.4.2",
                                      "MeV.1.4.21", "MeV.1.4.8", "MeV.3.30"]
    


    
    # # Check for mismatches before reordering
    # check_cluster_order(filtered_adata, custom_cluster_order)

    
    # Create dotplot of the top genes
    #create_dotplots_with_thresholds(adatas_filtered, pts_thresholds, clusters_to_remove, custom_cluster_order)

    # print("----")
    # print(adata.obs['leiden_fusion'].cat.categories.to_list())
    # print("----")

    no_filter_dotplot(filtered_adata, custom_cluster_order_no_filter)

    

    print("\n********\n* DONE *\n********")  # Indicate completion

    # Prints
    #print_gene_names(top_genes_names)
    #print_clusters(top_genes_cluster)