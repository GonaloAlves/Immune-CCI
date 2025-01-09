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


def umap_reso_cluster(adata, resolution_name):
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
    ax_ondata = sc.pl.umap(adata, color=resolution_name, title=f"UMAP - {resolution_name}",legend_loc = 'on data',legend_fontweight = 'bold', legend_fontsize = 8, legend_fontoutline = 1,return_fig=True)

    
    # Save the UMAP plot as an image (optional)
    output_path = f"umap_Meningeal_{resolution_name}_n.png"
    output_path_leg = f"umap_Meningeal_{resolution_name}_l.png"
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


# Step 7: Visualize the DotPlots of the DGE's
def create_dotplot(adata, top_genes_names, output_dir="dotplots_meningeal"):
    """
    Create and save a dotplot of the top genes per cluster.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    top_genes_names (dict): Dictionary of top gene names per cluster.
    output_dir (str): Directory to save the dotplot.
    filtered_adata = remove_NA_cat(adata)

    Returns:
    None
    """
    # #Create the directory if it doesn't exist
    # if os.path.exists(output_dir):
    #     shutil.rmtree(output_dir)
    # os.makedirs(output_dir)



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
        dendrogram='dendrogram_leiden_fusion',
        #dendrogram=False,
        return_fig=True
    )

    output_path = os.path.join(output_dir, "dotplot_0.3_dendro.png")
    dotplot.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap

    

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

# export to excel
def export_to_excel(top_genes_cluster, output_file="Meningeal_clusters.xlsx"):
    """
    Export the top_genes_cluster dictionary to an Excel file.

    Parameters:
    top_genes_cluster (dict): Dictionary containing DataFrames of top genes for each cluster.
    output_file (str): Path to the output Excel file.

    Returns:
    None
    """
    # Create an Excel writer object
    with pd.ExcelWriter(output_file) as writer: 
        # Loop through each cluster and its corresponding DataFrame
        for cluster, df in top_genes_cluster.items():
            # Write each DataFrame to a different sheet, named after the cluster
            df.to_excel(writer, sheet_name=cluster)

    # print(f"Data successfully exported to {output_file}")

def dendogram_sc(adata, output_dir="dendogram_meningeal"):
    """
    
    """

    # Create the directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Compute the dendrogram
    print(f"Computing dendrogram for leiden_fusion...")
    sc.tl.dendrogram(
        adata,
        groupby='leiden_fusion',
        use_rep= 'X_pca',
        cor_method= 'spearman',
        linkage_method='ward',
        use_raw=False
    )

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
    output_path = os.path.join(output_dir, f"leiden_fusion_dendro_0.3.png")
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap

# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy.h5ad")

    #print(adata)
    
    # Create cluster resolutions UMAP
    #umap_reso_cluster(adata, 'leiden_fusion')

    # Do a inicial filter in the a data
    filtered_adata = remove_NA_cat(adata)

    # Extract DGE data
    gene_names, logfoldchanges, pvals_adj, scores, pts = extract_dge_data(filtered_adata)
    
    # Create cluster DataFrames
    cluster_dfs = create_cluster_dfs(gene_names, logfoldchanges, pvals_adj, scores, pts, sort_by_logfc=True, pts_threshold=0.3)
    
    # Remove NA clusters
    cluster_dfs = remove_clusters_by_suffix(cluster_dfs, "NA")

    # Create dendogram ot the top genes
    dendogram_sc(filtered_adata)
    
    # Select the top genes for each cluster
    top_genes_cluster = select_top_genes(cluster_dfs)

    # Add the asterisk to cluster names with non-significant genes
    top_genes_cluster = addasterix(top_genes_cluster)
    
    # Collect top gene names for visualization
    top_genes_names = top_gene_names(top_genes_cluster)

    # Create dotplot of the top genes
    create_dotplot(filtered_adata, top_genes_names)

    print("Done")

    # export_to_excel(top_genes_cluster, output_file="top_genes_cluster_0.5.xlsx")

    # Prints
    #print_gene_names(top_genes_names)
    #print_clusters(top_genes_cluster)