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
    ax = sc.pl.umap(adata, color=resolution_name, title=f"UMAP - {resolution_name}", return_fig=True)
    ax_ondata = sc.pl.umap(adata, color=resolution_name, title=f"UMAP - {resolution_name}",legend_loc = 'on data',legend_fontweight = 'bold', legend_fontsize = 12, legend_fontoutline = 1,return_fig=True)

    
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
    
    dge_fusion = adata.uns['rank_genes_groups_leiden_fusion']
    
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


def create_dotplots_with_thresholds(adata, thresholds, output_dir="dotplots/immune/leiden_fusion"):
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

        # Create dendrogram (if not already present)
        print("   - Checking and creating dendrogram if necessary...")
        # plot_dendogram(adata)

        # Select the top genes for each cluster
        print("   - Selecting top genes for each cluster...")
        top_genes_cluster = select_top_genes(cluster_dfs)

        # Add an asterisk to clusters with non-significant genes
        top_genes_cluster = addasterix(top_genes_cluster)

        # Collect top gene names for visualization
        top_genes_names = top_gene_names(top_genes_cluster)

        # Reorder clusters (ON and OFF for dendrogram)
        print("   - Reordering clusters based on dendrogram...")
        # ordered_genes_dendro = reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram=True)
        #ordered_genes_no_dendro = reorder_clusters_to_dendrogram(adata, top_genes_names, dendrogram=False)

        print("   - Generating dotplots...")

        # # (1) With Dendrogram
        # dotplot_dendro = sc.pl.rank_genes_groups_dotplot(
        #     adata,
        #     var_names=ordered_genes_dendro,
        #     groupby='leiden_fusion_old1',
        #     key='rank_genes_groups_leiden_fusion_old1',
        #     cmap='bwr',
        #     vmin=-4,
        #     vmax=4,
        #     values_to_plot='logfoldchanges',
        #     colorbar_title='log fold change',
        #     use_raw=False,
        #     dendrogram='dendrogram_leiden_fusion_old1',
        #     return_fig=True
        # )

        # (2) Without Dendrogram
        dotplot_no_dendro = sc.pl.rank_genes_groups_dotplot(
            adata,
            var_names=top_genes_names,
            groupby='leiden_fusion',
            key='rank_genes_groups_leiden_fusion',
            cmap='Greys',
            vmin=0,
            vmax=1,
            values_to_plot='logfoldchanges',
            colorbar_title='log fold change',
            use_raw=False,
            dendrogram=False,
            return_fig=True
        )

        # dotplot_scaled_no_dendro = sc.pl.dotplot(
        #     adata,
        #     var_names=top_genes_names,
        #     groupby='leiden_fusion',
        #     cmap='Greys',
        #     use_raw=False,
        #     standard_scale='var',
        #     dendrogram=False,
        #     dot_min=0,
        #     dot_max=1,
        #     return_fig=True
        # )

        # Save plots
        # output_dendro = os.path.join(output_dir, f"dotplot_dendro_{threshold}.png")
        output_no_dendro = os.path.join(output_dir, f"dotplot_no_dendro_{threshold}.pdf")

        # dotplot_dendro.savefig(output_dendro, bbox_inches="tight")
        dotplot_no_dendro.savefig(output_no_dendro, bbox_inches="tight")

        plt.close()
        # print(f" Saved: {output_dendro}")
        print(f" Saved: {output_no_dendro}")


        # export_to_excel(top_genes_cluster, threshold)



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

def export_to_excel(top_genes_cluster, threshold, output_dir="excels/immune/updates"):
    """
    Export the top_genes_cluster dictionary to an Excel file, using the threshold in the filename.

    Parameters:
    top_genes_cluster (dict): Dictionary containing DataFrames of top genes for each cluster.
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

    print(f"\n Exporting top genes to Excel: {output_file}")

    # Create an Excel writer object
    with pd.ExcelWriter(output_file) as writer: 
        # Loop through each cluster and its corresponding DataFrame
        for cluster, df in top_genes_cluster.items():
            # Write each DataFrame to a different sheet, named after the cluster
            df.to_excel(writer, sheet_name=cluster)

    print(f"Excel file saved: {output_file}")



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
    output_path = os.path.join(output_dir, f"leiden_fusion_dendro.png")
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


# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad")

    filtered_adata = remove_NA_cat(adata)

    #print(adata)
    print(adata.obs['leiden_fusion'].cat.categories.to_list())
    print(adata.obs['leiden_fusion'])
    #print(adata.obs['leiden_fusion_old1'])
    
    # #Create cluster resolutions UMAP
    umap_reso_cluster(filtered_adata, 'leiden_fusion')

    pts_thresholds = [0.3, 0.4, 0.5]

    # Create dotplot of the top genes
    create_dotplots_with_thresholds(filtered_adata, pts_thresholds)

    print("\n********\n* DONE *\n********")  

    # Prints
    #print_gene_names(top_genes_names)
    #print_clusters(top_genes_cluster)