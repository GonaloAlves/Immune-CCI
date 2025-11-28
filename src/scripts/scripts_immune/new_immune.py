# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Step 0: Extract all resolutions
def umap_reso_cluster(adata, resolution_name, output_dir="reso/reso_immune/updates"):
    """
    Plot UMAP for a specific resolution with cluster numbers displayed at the centroid of each cluster.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    resolution_name (str): The resolution name to be used for coloring the UMAP plot (e.g., 'leiden_fusion').

    Returns:
    None
    """

    # Plot UMAP for the specified resolution
    ax = sc.pl.umap(adata, 
                    color=resolution_name, 
                    title=f"UMAP - {resolution_name}",
                    size=20, 
                    return_fig=True)
    
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
    output_path= os.path.join(output_dir, f"umap_Immune_{resolution_name}_n.pdf") 
    output_dir_leg = os.path.join(output_dir, f"umap_Immune_{resolution_name}_l.pdf")
    ax.figure.savefig(output_path, bbox_inches="tight")
    ax_ondata.figure.savefig(output_dir_leg, bbox_inches="tight")
    # print(f"UMAP plot saved as {output_path}")
    plt.close()  # Close the plot to avoid overlap



def load_data(file_path):
    """
    Load the data from a .h5ad file.
    """
    print(f"\n🔹 Loading h5ad file: {file_path}")
    adata = sc.read_h5ad(file_path)
    print(f"Data loaded successfully! Shape: {adata.shape}")
    return adata


def extract_dge_data(adata):
    """
    Extract differential gene expression (DGE) data from the AnnData object.
    """
    print("\n🔹 Extracting DGE data...")
    
    dge_fusion = adata.uns['rank_genes_groups_leiden_fusion_old1']
    
    gene_names = pd.DataFrame(dge_fusion['names'])
    logfoldchanges = pd.DataFrame(dge_fusion['logfoldchanges'])
    pvals_adj = pd.DataFrame(dge_fusion['pvals_adj'])
    scores = pd.DataFrame(dge_fusion['scores'])
    pts = pd.DataFrame(dge_fusion['pts'])

    print(f"Extracted {gene_names.shape[0]} genes across {gene_names.shape[1]} clusters.")
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

    print(f"\nConverting into Dataframes")
    for i in range(len(gene_names.columns)):  # For to read each cluster
        cluster_names = gene_names.columns[i]
        print(f"   - Processing cluster: {cluster_names}")
        
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
        cluster_dfs[cluster_names] = filtered_df  # Store the respective clusters in the dictionary
        print(f"   Cluster {cluster_names}: {filtered_df.shape[0]} genes after filtering")
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
    print(f"\nRemoving {suffix} clusters")
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
    print(f"\nSelecting statisticaly significant genes from clusters")
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
    top_genes_negative = top_genes[top_genes['logfoldchange'] < 0]
    
    for gene, values in top_genes_negative.iterrows():
        print(f"\nGenes flagged: {gene}")

    replacements = []
    for gene, values in df.iterrows():
        if gene not in top_genes_positive.index:
            replacements.append(values) # Prevent of appending duplicated genes
            print(f"\nGenes added to DF: {gene} ")
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


def select_top_genes(cluster_dfs: dict):
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
        print(f"\nStarting the gene filtering in {cluster}")
        top_genes = filter_pval(df) # First remove the not significant genes 
        print(f"\nNon significative genes removed from {cluster} ")

        print(f"\nChecking if {cluster} has negative genes")
        if has_neg_lfc(top_genes): # If there is any neg lfc genes
            print(f"\n{cluster} has negative genes")
            num_nega_lfc = count_neg_lfc(top_genes)
            print(f"\nReplacing genes in {cluster}")
            top_genes_cluster[cluster] = replace_neg_lfc(df, top_genes, num_nega_lfc) # replace the neg genes
        elif top_genes.empty:
            print(f"\n{cluster} has no positive genes, replacting with new genes")
            top_genes_cluster[cluster] = replace_all(df) # replace all if it is empty
        else:
            print(f"\n{cluster} has no negative flagged genes")
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

def no_filter_dotplot(adata, cluster_order, output_dir="dotplots/immune/leiden_fusion"):


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
    output_normal = os.path.join(output_dir, f"dotplot_dge_no_filter.pdf")
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


def create_dotplots_with_thresholds(adata, thresholds, cluster_order, output_dir="dotplots/immune/leiden_fusion"):
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

    adata.obs['leiden_fusion_old1'] = adata.obs['leiden_fusion_old1'].astype('category')
    adata.obs['leiden_fusion_old1'] = adata.obs['leiden_fusion_old1'].cat.reorder_categories(cluster_order, ordered=True)

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

        print(cluster_dfs)
        
        ordered_clusters = order_clusters(cluster_dfs, cluster_order)

        #export_to_excel_all(ordered_clusters, threshold,string = "all")

        # Select the top genes for each cluster
        print("   - Selecting top genes for each cluster...")
        top_genes_cluster = select_top_genes(ordered_clusters)

        # Add an asterisk to clusters with non-significant genes
        #top_genes_cluster = addasterix(top_genes_cluster)

        ordered_clusterss = order_clusters(top_genes_cluster, cluster_order)

        #export_to_excel_all(ordered_clusterss, threshold, string = "top")

        # Collect top gene names for visualization
        top_genes_names = top_gene_names(ordered_clusterss)

        top_genes_names = order_clusters(top_genes_names, cluster_order)

        # Reorder clusters (ON and OFF for dendrogram)
        print("   - Reordering clusters based on dendrogram...")
        # ordered_genes_dendro = reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram=True)
        #ordered_genes_no_dendro = reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram=False)

        print("   - Generating dotplots...")

        # (1) With Dendrogram
        dotplot_normal = sc.pl.rank_genes_groups_dotplot(
            adata,
            var_names=top_genes_names,
            groupby='leiden_fusion_old1',
            key='rank_genes_groups_leiden_fusion_old1',
            cmap='bwr',
            vmin=-4,
            vmax=4,
            values_to_plot='logfoldchanges',
            colorbar_title='Magnitude of expression',
            use_raw=False,
            dendrogram=False,
            return_fig=True
        )

        # Save plots
        output_normal = os.path.join(output_dir, f"dotplot_dge_filter_{threshold}.pdf")

        dotplot_normal.savefig(output_normal, bbox_inches="tight")

        plt.close()
        print(f" Saved: {output_normal}")



def print_gene_names(top_genes_names):

    print("Printing top genes from each cluster")
    for cluster, df in top_genes_names.items():
        print(f"\nTop genes in {cluster}:")
        print(df)

def print_clusters(top_genes_cluster):

    print("Printing top genes from each cluster")
    for cluster, df in top_genes_cluster.items():
        print(f"\nTop genes in {cluster}:")
        print(df)

def remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    mask_NA = adata.obs['leiden_fusion_old1'] != 'Imm.NA' #creates mask for remove NA cells
    #print(mask_NA)    
    adata2 = adata[mask_NA] #apply mask

    # print(adata2.obs['leiden_fusion'].cat.categories.to_list())
    # adata2.obs['leiden_fusion'] = adata2.obs['leiden_fusion'].cat.remove_unused_categories()
    # print(adata2.obs['leiden_fusion'].cat.categories.to_list())
    # print(adata.obs['leiden_fusion'])

    return adata2

def nonsigngene(top_genes):
    
    return (top_genes['pvals_adj'] > 0.05 ).any()

def addasterix(top_genes_cluster):
    
    updated_cluster = {}
    print("Adding asterisk in non significative genes")

    for cluster, df in top_genes_cluster.items():
        # Check if the cluster has any non-significant genes
        if nonsigngene(df):
            # Add an asterisk to the cluster name
            print (f"\n{cluster} got targeted")
            updated_cluster[cluster + '*'] = df
        else:
            updated_cluster[cluster]  = df

    return updated_cluster

def export_to_excel(top_genes_cluster, threshold, output_dir="excels/immune/new_tese"):
    """
    Export all genes per cluster to an Excel file, using the threshold in the filename.

    Parameters:
    top_genes_cluster (dict): Dictionary containing DataFrames of genes for each cluster.
    threshold (float): The threshold value used for filtering.
    output_dir (str): Directory where the Excel file will be saved.

    Returns:
    None
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Define file path dynamically using the threshold
    output_file = os.path.join(output_dir, f"top_genes_cluster_{threshold}.xlsx")

    print(f"\n Exporting genes to Excel: {output_file}")

    # Create an Excel writer object
    with pd.ExcelWriter(output_file) as writer: 
        # Loop through each cluster and its corresponding DataFrame
        for cluster, df in top_genes_cluster.items():
            # Write each DataFrame to a different sheet, named after the cluster
            df.to_excel(writer, sheet_name=cluster)

    print(f"Excel file saved: {output_file}")

def export_to_excel_all(cluster_dfs, threshold,string, output_dir="excels/immune/new_tese/dge"):
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

def plot_dendogram(adata, output_dir="dendogram_immune"):
    """
    
    """

    # Create the directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Plot the dendrogram
    print(f"Plotting dendrogram for leiden_fusion...")
    sc.pl.dendrogram(
        adata,
        groupby='leiden_fusion_old1',
        dendrogram_key='dendrogram_leiden_fusion_old1',
        #dendrogram_key=None
        orientation='top',
        show=False
    )

    # Save the plot
    output_path = os.path.join(output_dir, f"leiden_fusion_dendro.pdf")
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap
    
def reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram, dendrogram_key = 'dendrogram_leiden_fusion_old1'):
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

def order_clusters(top_genes_dict, custom_order):
    """
    Reorder the dictionary keys according to a custom order.

    Parameters:
    top_genes_dict (dict): Dictionary mapping cluster names to gene lists.
    custom_order (list): List defining the desired cluster order.

    Returns:
    dict: A new dictionary ordered by the custom order.
    """
    print("\n🔃 Ordering clusters...")

    # 1️⃣ Add clusters in the custom order (only if they exist)
    ordered_dict = {k: top_genes_dict[k] for k in custom_order if k in top_genes_dict}

    # 2️⃣ Add any remaining clusters not present in the custom order
    remaining = [k for k in top_genes_dict.keys() if k not in custom_order]
    if remaining:
        print(f"ℹ️ Adding {len(remaining)} clusters not found in custom order: {remaining}")
        for k in remaining:
            ordered_dict[k] = top_genes_dict[k]

    print(f"📊 Final ordered cluster count: {len(ordered_dict)}")
    return ordered_dict

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
        groupby='leiden_fusion_old1',
        use_rep='X_pca',
        cor_method='spearman',
        linkage_method='ward',
        use_raw=False
    )

# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad")

    filtered_adata = remove_NA_cat(adata)

    print(adata)
    # print(adata.obs['leiden_fusion'].cat.categories.to_list())
    # print(adata.obs['leiden_fusion'])
    print(adata.obs['leiden_fusion_old1'])
    
    
    # #Create cluster resolutions UMAP
    # umap_reso_cluster(filtered_adata, 'leiden_fusion')

    pts_thresholds = [0.4]

    custom_cluster_order = [
    "Imm.M0Like.0", "Imm.M0Like.1", 
    "Imm.M0Like.2", "Imm.MHCII.0" ,
    "Imm.Interferon.0", "Imm.DAM.0", 
    "Imm.DAM.1", "Imm.PVM.0", "Imm.Proliferative.0"]

    custom_cluster_order_no_filter = [
    "Imm.0.8.0", "Imm.0.8.1", "Imm.0.8.2", 
    "Imm.0.8.3", 
    "Imm.1.2.13", 
    "Imm.1.2.14",
    "Imm.0.8.6", "Imm.1.2.12", "Imm.1.2.5", "Imm.1.2.4",
    "Imm.1.2.15"]   

    # dendogram_sc(filtered_adata)
    # plot_dendogram(filtered_adata)

    # Create dotplot of the top genes
    create_dotplots_with_thresholds(filtered_adata, pts_thresholds, custom_cluster_order_no_filter)
    #no_filter_dotplot(filtered_adata, custom_cluster_order_no_filter)

    print("\n********\n* DONE *\n********")  

    # Prints
    #print_gene_names(top_genes_names)
    #print_clusters(top_genes_cluster)
