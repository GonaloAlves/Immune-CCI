# def run_gsea_for_cluster(tstat_file, gene_sets, output_dir):
#     """
#     Run GSEA analysis for a given t-statistic file.

#     Parameters:
#     - tstat_file (str): Path to the t-stat Excel file.
#     - gene_sets (str): Gene set database to use in the GSEA analysis.
#     - output_dir (str): Directory to save the GSEA results.
#     """
#     # Load the t-stat file into a DataFrame
#     df = pd.read_excel(tstat_file)

#     # Prepare the data for GSEA analysis
#     # Extract Gene Names and t-statistics for ranking
#     gsea_data = df[['Gene.Name', 't']]  # Replace 't' with the correct column containing t-statistics
    
#     # Sort by t-statistic for ranking
#     gsea_data = gsea_data.sort_values(by='t', ascending=False)

#     # Define output directory for the cluster
#     cluster_name = os.path.basename(tstat_file).split("_tstat")[0]
#     cluster_output_dir = os.path.join(output_dir, cluster_name)
#     if not os.path.exists(cluster_output_dir):
#         os.makedirs(cluster_output_dir)

#     # Run GSEA analysis using the t-statistics
#     print(f"Running GSEA for cluster {cluster_name}...")
#     results = gp.prerank(
#         rnk=gsea_data.values,  # Use Gene Name and t-statistics as ranking
#         gene_sets=gene_sets,
#         outdir=cluster_output_dir,  # Directory to save GSEA results
#         min_size=15,        # Minimum gene set size to consider
#         max_size=500,       # Maximum gene set size to consider
#         permutation_num=100,  # Number of permutations
#         seed=42  # For reproducibility
#     )
    
#     print(f"GSEA results saved for {tstat_file} in {cluster_output_dir}")
#     # Clean up memory
#     gc.collect()


# # Main function to process all t-stat files
# def start_gsea_tstat():
#     """
#     Start GSEA analysis for all clusters using available t-statistic files.
#     """
#     # Iterate through each t-statistic file and run GSEA analysis
#     for tstat_file in tstat_files:
#         # Run GSEA for each t-stat file
#         run_gsea_for_cluster(tstat_file, gene_sets, gsea_dir)