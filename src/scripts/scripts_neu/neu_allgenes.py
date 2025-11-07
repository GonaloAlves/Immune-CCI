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
    return filtered_df

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

        # if has_neg_lfc(top_genes): # If there is any neg lfc genes
        #     num_nega_lfc = count_neg_lfc(top_genes)
        #     top_genes_cluster[cluster] = replace_neg_lfc(df, top_genes, num_nega_lfc) # replace the neg genes
        # elif top_genes.empty:
        #     top_genes_cluster[cluster] = replace_all(df) # replace all if it is empty
        # else:
        top_genes_cluster[cluster] = top_genes # if the cluster doesnt have neg genes just continues

    return top_genes_cluster


# export to excel
def export_to_excel(top_genes_cluster, threshold, output_dir="excels/neuron/updates"):
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


def start(adata, thresholds, output_dir="excels/neuron/updates"):
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

        # Select the top genes for each cluster
        print("   - Selecting top genes for each cluster...")
        top_genes_cluster = select_top_genes(cluster_dfs)


        export_to_excel(top_genes_cluster, threshold)

    # Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Neu_CentralCanal_raw_norm_ranked_copy_copy.h5ad")

    print(adata)

    pts_thresholds = [0, 0.2, 0.3, 0.4, 0.5]

    # Create dotplot of the top genes
    start(adata, pts_thresholds)
    

    print("\n********\n* DONE *\n********")  # Indicate completion

    # Prints
    #print_gene_names(top_genes_names)
    #print_clusters(top_genes_cluster)