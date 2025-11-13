# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# This script will display the gene expression values of the canonical genes from Meningeal dataset


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


def create_dotplots_with_thresholds(adata, genes, thresholds, cluster_order, gene, output_dir="canonical/canonical_meningeal/updated_pts4"):
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

    #print(adata['leiden_fusion'].cat.categories.to_list())
    print(adata)
    # Ensure leiden_fusion is categorical and reorder it
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].astype('category')
    adata.obs['leiden_fusion'] = adata.obs['leiden_fusion'].cat.reorder_categories(cluster_order, ordered=True)



    for threshold in thresholds:
        print(f"\nProcessing pts threshold: {threshold}")

        # Extract DGE data
        gene_names, logfoldchanges, pts = extract_dge_data(adata)

        # Create cluster DataFrames with the current threshold
        cluster_dfs = create_cluster_dfs(gene_names, logfoldchanges, pts, pts_threshold=threshold)
        
        # Remove NA clusters
        cluster_dfs = remove_clusters_by_suffix(cluster_dfs, "NA")

        # Compare canonical genes with cluster-specific genes
        filtered_genes = compare_canonical(genes, cluster_dfs)

        # Aggregate filtered genes by gene group
        top_genes_names = top_gene_names(filtered_genes, genes)

        # Example user-defined gene group order
        user_gene_group_order = ["Endothelial","Epithelial","SMC","Pericytes" ,"VLMC", 
                                  "Fibroblasts","ECM_Laminin","ECM_Collagen","ECM_Secreted_Signaling",
                                  "Proliferative"]

        # Reorder the dictionary based on user order
        top_genes_names = {key: top_genes_names[key] for key in user_gene_group_order}

        excel_order = { 
            "Endothelial": ["MeV.Endothelial.0","MeV.Endothelial.1","MeV.Endothelial.2","MeV.Endothelial.3"],
            "Epithelial": ["MeV.Epithelial.0"],
            "SMC": ["MeV.SMC.0"],
            "Pericytes": ["MeV.Pericytes.0"],
            "VLMC": ["MeV.VLMC.0","MeV.VLMC.1"],
            "Fibroblasts": ["MeV.FibCollagen.0","MeV.FibCollagen.1", "MeV.FibCollagen.2","MeV.FibCollagen.3","MeV.FibLaminin.0","MeV.Fib.0","MeV.Fib.1",
                            "MeV.Fib.2","MeV.Fib.5","MeV.Fib.3","MeV.Fib.4"],
            "ECM_Laminin": ["MeV.FibLaminin.0"],
            "ECM_Collagen": ["MeV.FibCollagen.0","MeV.FibCollagen.1", "MeV.FibCollagen.2","MeV.FibCollagen.3"],
            #"ECM_Secreted_Signaling": ["MeV.FibProlif.0"],
            "Proliferative": ["MeV.FibProlif.0"]
        }

        # Export only the genes that match these clusters per group
        export_canonical_to_excel(filtered_genes, excel_order, threshold)

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

        # # (2) Scaled expression (Greys) with dendrogram
        # dotplot_scaled_dendro = sc.pl.dotplot(
        #     adata,
        #     var_names=top_genes_names,
        #     groupby='leiden_fusion',
        #     cmap='Greys',
        #     colorbar_title='Scaled expression',
        #     use_raw=False,
        #     standard_scale='var',
        #     dendrogram='dendrogram_leiden_fusion',
        #     return_fig=True
        # )

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

        # # (4) Raw expression (Reds) with dendrogram
        # dotplot_normal_dendro = sc.pl.dotplot(
        #     adata,
        #     var_names=top_genes_names,
        #     groupby='leiden_fusion',
        #     cmap='Reds',
        #     use_raw=False,
        #     dendrogram='dendrogram_leiden_fusion',
        #     return_fig=True
        # )

        # Save dotplots with appropriate filenames
        output_scaled_no_dendro = os.path.join(output_dir, f"{gene}dotplot_scaled_no_dendro_{threshold}.pdf")
        #output_scaled_dendro = os.path.join(output_dir, f"{gene}dotplot_scaled_dendro_{threshold}.png")
        output_normal_no_dendro = os.path.join(output_dir, f"{gene}dotplot_normal_no_dendro_{threshold}.pdf")
        #output_normal_dendro = os.path.join(output_dir, f"{gene}dotplot_normal_dendro_{threshold}.png")

        dotplot_scaled_no_dendro.savefig(output_scaled_no_dendro, bbox_inches="tight")
        #dotplot_scaled_dendro.savefig(output_scaled_dendro, bbox_inches="tight")
        dotplot_normal_no_dendro.savefig(output_normal_no_dendro, bbox_inches="tight")
        #dotplot_normal_dendro.savefig(output_normal_dendro, bbox_inches="tight")

        plt.close()
        print(f"Saved dotplots for threshold {threshold}:")
        print(f"  - {output_scaled_no_dendro}")
        #print(f"  - {output_scaled_dendro}")
        print(f"  - {output_normal_no_dendro}")
        #print(f"  - {output_normal_dendro}")


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
    logfoldchanges = pd.DataFrame(dge_fusion['logfoldchanges'])
    pts = pd.DataFrame(dge_fusion['pts'])
    
    return gene_names, logfoldchanges, pts


# Step 3: Create cluster dataframes with filtered data
def create_cluster_dfs(gene_names, logfoldchanges, pts, pts_threshold=0):
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

         # Create dataframe for this cluster
        cluster_df = pd.DataFrame({
            'logfoldchange': logfoldchanges.iloc[:, i].values,
            'pts': pts_reindexed.values
        }, index=gene_reindex.values)

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

def export_cluster_cell_counts(adata, output_dir="excels/meningeal/updates"):
    """
    Export the number of cells per cluster to an Excel file.
    
    Parameters:
    - adata (AnnData): Annotated data matrix.
    - output_dir (str): Directory where the Excel file will be saved.
    
    Returns:
    - None
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Count the number of cells per cluster
    cell_counts = adata.obs['leiden_fusion'].value_counts().reset_index()
    cell_counts.columns = ['Cluster', 'Cell Count']
    
    # Define output file path
    output_file = os.path.join(output_dir, "cluster_cell_counts_all.xlsx")
    
    # Save to Excel
    print(f"\n Exporting cell counts to: {output_file}")
    cell_counts.to_excel(output_file, index=False)
    print(f"Saved Excel file: {output_file}")


def export_canonical_to_excel(filtered_genes, excel_order, threshold, output_dir="excels/meningeal/new_tese/canonical"):
    """
    Export canonical gene expression results to an Excel file,
    showing only the clusters defined in 'excel_order' for each gene group.
    """

    print("\n🟦 [CHECKPOINT 1] Starting export_canonical_to_excel")
    print(f"Output dir: {output_dir}")
    print(f"Threshold: {threshold}")
    print(f"Filtered gene groups: {list(filtered_genes.keys())}")
    print(f"Excel order groups: {list(excel_order.keys())}")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("🟩 [CHECKPOINT 2] Output directory created.")
    else:
        print("🟩 [CHECKPOINT 2] Output directory already exists.")

    output_file = os.path.join(output_dir, f"canonical_genes_filtered_{threshold}.xlsx")
    print(f"🟦 [CHECKPOINT 3] Preparing Excel file: {output_file}")

    try:
        with pd.ExcelWriter(output_file) as writer:
            for group_name, allowed_clusters in excel_order.items():
                print(f"\n🔹 [CHECKPOINT 4] Processing group: {group_name}")

                if group_name not in filtered_genes:
                    print(f"⚠️ Skipping {group_name}: not found in filtered_genes")
                    continue

                group_clusters = filtered_genes[group_name]
                print(f"   ➜ Available clusters for {group_name}: {list(group_clusters.keys())}")

                combined_data = []

                for cluster_name in allowed_clusters:
                    print(f"   ➜ Checking cluster: {cluster_name}")

                    if cluster_name in group_clusters:
                        df = group_clusters[cluster_name].copy()
                        df["Cluster"] = cluster_name
                        combined_data.append(df)
                        print(f"      ✅ Added {len(df)} genes from {cluster_name}")
                    else:
                        print(f"      ⚠️ No data for {cluster_name} in {group_name}")

                if combined_data:
                    print(f"   🔸 [CHECKPOINT 5] Combining {len(combined_data)} DataFrames for {group_name}")
                    try:
                        group_df = pd.concat(combined_data)
                        cols = ["Cluster"] + [col for col in group_df.columns if col != "Cluster"]
                        group_df = group_df[cols]
                        group_df.to_excel(writer, sheet_name=group_name)
                        print(f"   ✅ Wrote sheet '{group_name}' with {group_df.shape[0]} rows.")
                    except Exception as e:
                        print(f"   ❌ ERROR writing sheet {group_name}: {e}")
                else:
                    print(f"   ⚠️ No combined data for {group_name} — writing info sheet instead.")
                    pd.DataFrame(
                        {"Info": [f"No genes found for these clusters at threshold {threshold}"]}
                    ).to_excel(writer, sheet_name=group_name)

        print(f"\n✅ [CHECKPOINT 6] Excel file saved successfully: {output_file}")

    except Exception as e:
        print(f"❌ [ERROR] Failed while writing Excel file: {e}")
        import traceback
        traceback.print_exc()



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
                           title=f"UMAP - Meningeal Clusters",
                           legend_loc = 'on data',
                           legend_fontsize=7,  
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


# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")

    # print(adata.obs["leiden_fusion"])
    # print(adata.uns["rank_genes_groups_leiden_fusion"])

    filtered_adata = remove_NA_cat(adata)

    # gene_filtered_adata = filter_cells_by_gene_expression(filtered_adata, "Mylip")

    clusters_to_remove = ['MeV.ImmuneDoublets.0', 'MeV.LowQuality.0', 'MeV.EndoUnknow.4', 'MeV.FibUnknown.6']
    adatas_filtered = remove_clusters(filtered_adata, clusters_to_remove)

    #Create cluster resolutions UMAP
    umap_reso_cluster(adatas_filtered, 'leiden_fusion')

    # #preform dendrogram
    # dendogram_sc(adatas_filtered)

    # # Load canonical gene lists from a directory
    # canonical_genes_dir = "/home/makowlg/Documents/Immune-CCI/src/canonical/canonical_txt/Meningeal"
    # genes = load_canonical_from_dir(canonical_genes_dir)


    # # Define thresholds
    # pts_thresholds = [0, 0.2, 0.3]

    # custom_cluster_order = ["MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Epithelial.0",
    #                         "MeV.SMC.0", "MeV.Pericytes.0", "MeV.VLMC.0", "MeV.VLMC.1" , "MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3",
    #                         "MeV.FibLaminin.0", "MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2", "MeV.Fib.5", "MeV.Fib.3", "MeV.Fib.4", "MeV.FibProlif.0"]
    

    # # Check for mismatches before reordering
    # check_cluster_order(adatas_filtered, custom_cluster_order)
    
    # # Generate dotplots for each threshold
    # create_dotplots_with_thresholds(adatas_filtered, genes, pts_thresholds, custom_cluster_order, "")

    # #export_cluster_cell_counts(adatas_filtered)


    