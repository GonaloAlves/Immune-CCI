# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Load the dataset
def load_data(file_path):
    """
    Load the single-cell data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    return sc.read_h5ad(file_path)


# Step 2: Extract Differentially Expressed Genes (DEGs)
def extract_dge_data(adata):
    """
    Extract the differential gene expression data from the AnnData object.

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
def create_cluster_dfs(gene_names, logfoldchanges, pvals_adj, scores, pts, sort_by_logfc=True, pts_threshold=0.5):
    """
    Create a dictionary of dataframes per cluster, with filtered and optionally sorted genes.

    Parameters:
    gene_names (DataFrame): DataFrame with gene names.
    logfoldchanges (DataFrame): DataFrame with log fold changes for each gene.
    pvals_adj (DataFrame): DataFrame with adjusted p-values for each gene.
    scores (DataFrame): DataFrame with scores for each gene.
    pts (DataFrame): DataFrame with "pts" values for each gene.
    sort_by_logfc (bool): If True, sort the filtered DataFrame by log fold change in descending order.
    pts_threshold (float): The minimum 'pts' value to filter genes.

    Returns:
    dict: Dictionary of filtered DataFrames for each cluster.
    """
    cluster_dfs = {}

    for i in range(len(gene_names.columns)):  # For each cluster
        gene_reindex = gene_names.iloc[:, i]  # Reindex genes for pts
        pts_reindexed = pts.iloc[:, i].reindex(gene_reindex.values)  # Reindex the pts

        # Create dataframe with relevant data
        cluster_df = pd.DataFrame({
            'score': scores.iloc[:, i].values,
            'logfoldchange': logfoldchanges.iloc[:, i].values,
            'pvals_adj': pvals_adj.iloc[:, i].values,
            'pts': pts_reindexed.values},
            index=gene_reindex.values
        )

        # Filter by 'pts' using the user-specified threshold
        mask_pts = (cluster_df['pts'] >= pts_threshold)
        filtered_df = cluster_df[mask_pts]

        # Optionally sort by log fold change if sort_by_logfc is True
        if sort_by_logfc:
            filtered_df = filtered_df.sort_values(by='logfoldchange', ascending=False)

        cluster_names = gene_names.columns[i]
        cluster_dfs[cluster_names] = filtered_df  # Store in dictionary

    return cluster_dfs


# Step 4: Remove any cluster
def remove_clusters_by_suffix(cluster_dfs, suffix):
    """
    Remove clusters whose names end with a specific suffix from the cluster_dfs dictionary.

    Parameters:
    cluster_dfs (dict): Dictionary of clusters and their data.
    suffix (str): The suffix that determines which clusters to remove (e.g., 'NA').

    Returns:
    dict: Updated dictionary with specified clusters removed.
    """
    clusters_to_delete = [cluster for cluster in cluster_dfs if cluster.endswith(suffix)]

    for cluster in clusters_to_delete:
        del cluster_dfs[cluster]

    return cluster_dfs


# Step 5: Filter, sort, and select top 5 genes per cluster

def filter_pval(df):
    """
    Filter genes by p-value and select the top 5 genes.

    Parameters:
    df (DataFrame): DataFrame containing genes for a cluster.

    Returns:
    DataFrame: Filtered DataFrame with top 5 genes.
    """
    mask_pval = df['pvals_adj'] < 0.05
    filtered_df = df[mask_pval]
    return filtered_df.head(5)


def has_neg_lfc(top_genes):
    """
    Check if the top genes contain any negative log fold changes.

    Parameters:
    top_genes (DataFrame): DataFrame of the top 5 genes for a cluster.

    Returns:
    bool: True if negative log fold change exists, False otherwise.
    """
    return (top_genes['logfoldchange'] < 0).any()


def count_neg_lfc(top_genes):
    """
    Count the number of genes with negative log fold changes.

    Parameters:
    top_genes (DataFrame): DataFrame of the top 5 genes for a cluster.

    Returns:
    int: Number of genes with negative log fold changes.
    """
    return (top_genes['logfoldchange'] < 0).sum()


def replace_neg_lfc(df, top_genes, num_nega_lfc):
    """
    Replace genes with negative log fold change with new genes.

    Parameters:
    df (DataFrame): Original DataFrame of the cluster.
    top_genes (DataFrame): DataFrame of the top 5 genes for the cluster.
    num_nega_lfc (int): Number of genes with negative log fold changes.

    Returns:
    DataFrame: DataFrame with negative log fold change genes replaced.
    """
    pos_lfc = df.head(num_nega_lfc)  # Replace with genes from the original df
    top_genes_positive = top_genes[top_genes['logfoldchange'] > 0]

    replacements = []
    for gene, values in df.iterrows():
        if gene not in top_genes_positive.index:
            replacements.append(values)
        if len(replacements) == num_nega_lfc:
            break  # Stop once enough replacement genes are found

    replacements_df = pd.DataFrame(replacements)
    return pd.concat([top_genes_positive, replacements_df])


def replace_all(df):
    """
    Select the top 5 genes as a fallback if no filtered genes are found.

    Parameters:
    df (DataFrame): Original DataFrame of the cluster.

    Returns:
    DataFrame: DataFrame with the top 5 fallback genes.
    """
    return df.head(5)


def select_top_genes(cluster_dfs):
    """
    Main function to select the top 5 genes per cluster, replacing genes with negative log fold change.

    Parameters:
    cluster_dfs (dict): Dictionary containing dataframes of each cluster.

    Returns:
    dict: Dictionary containing the top 5 genes per cluster.
    """
    top_genes_cluster = {}
    
    for cluster, df in cluster_dfs.items():
        top_genes = filter_pval(df)

        if has_neg_lfc(top_genes):
            num_nega_lfc = count_neg_lfc(top_genes)
            top_genes_cluster[cluster] = replace_neg_lfc(df, top_genes, num_nega_lfc)
        elif top_genes.empty:
            top_genes_cluster[cluster] = replace_all(df)
        else:
            top_genes_cluster[cluster] = top_genes

    return top_genes_cluster

# def select_top_genes(cluster_dfs):
#     """
#     Select the top 5 genes for each cluster, replacing genes with negative log fold change.

#     Parameters:
#     cluster_dfs (dict): Dictionary containing dataframes of each cluster.

#     Returns:
#     dict: Dictionary containing the top 5 genes per cluster.
#     """
#     top_genes_cluster = {}

#     for cluster, df in cluster_dfs.items():
#         mask_pval = df['pvals_adj'] < 0.05  # Filter by p-value
#         filtered_df = df[mask_pval]
#         top_genes = filtered_df.head(5)

#         if (top_genes['logfoldchange'] < 0).any():  # Check for negative log fold changes
#             num_nega_lfc = (top_genes['logfoldchange'] < 0).sum()
#             pos_lfc = df.head(num_nega_lfc)

#             mask1 = top_genes['logfoldchange'] > 0
#             top_genes_positive = top_genes[mask1]

#             replacements = []
#             for gene, values in df.iterrows():
#                 if gene not in top_genes_positive.index:
#                     replacements.append(values)

#                 if len(replacements) == num_nega_lfc:
#                     break  # Stop once we have enough replacement genes

#             replacements_df = pd.DataFrame(replacements)

#             combined_top_genes = pd.concat([top_genes_positive, replacements_df])
#             top_genes_cluster[cluster] = combined_top_genes

#         elif top_genes.empty:
#             top_genes_new = df.head(5)
#             top_genes_cluster[cluster] = top_genes_new

#         else:
#             top_genes_cluster[cluster] = top_genes

#     return top_genes_cluster

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


# Step 7: Visualize the DotPlots of the DGE's
def create_dotplot(adata, top_genes_names, output_dir="dotplots"):
    """
    Create and save a dotplot of the top genes per cluster.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    top_genes_names (dict): Dictionary of top gene names per cluster.
    output_dir (str): Directory to save the dotplot.

    Returns:
    None
    """
    # Create the directory if it doesn't exist
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Generate the dotplot
    dotplot = sc.pl.rank_genes_groups_dotplot(
        adata,
        var_names=top_genes_names,
        groupby='leiden_fusion',
        key='rank_genes_groups_leiden_fusion',
        cmap='bwr',
        vmin=-4,
        vmax=4,
        values_to_plot='logfoldchanges',
        colorbar_title='log fold change',
        use_raw=False,
        dendrogram=False,
        return_fig=True
    )

    output_path = os.path.join(output_dir, "dotplot_1.png")
    dotplot.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap


# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad")
    
    # Extract DGE data
    gene_names, logfoldchanges, pvals_adj, scores, pts = extract_dge_data(adata)
    
    # Create cluster DataFrames
    cluster_dfs = create_cluster_dfs(gene_names, logfoldchanges, pvals_adj, scores, pts, sort_by_logfc=True, pts_threshold=0.3)    
    
    # Remove NA clusters
    cluster_dfs = remove_clusters_by_suffix(cluster_dfs, "NA")
    
    # Select the top genes for each cluster
    top_genes_cluster = select_top_genes(cluster_dfs)
    
    # Collect top gene names for visualization
    top_genes_names = top_gene_names(top_genes_cluster)
    
    # Create dotplot of the top genes
    create_dotplot(adata, top_genes_names)


#test