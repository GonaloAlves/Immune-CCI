import pandas as pd

# File paths for source files
file_0_3 = '/home/makowlg/Documents/Immune-CCI/src/top_genes_cluster_0.3.xlsx'
file_0_4 = '/home/makowlg/Documents/Immune-CCI/src/top_genes_cluster_0.4.xlsx'
file_0_5 = '/home/makowlg/Documents/Immune-CCI/src/top_genes_cluster_0.5.xlsx'
output_file = '/home/makowlg/Documents/Immune-CCI/src/output_excel.xlsx'  
merged_output_file = '/home/makowlg/Documents/Immune-CCI/src/merged_top_genes.xlsx'

# Load the top genes Excel files
df_0_3 = pd.read_excel(file_0_3, sheet_name=None)
df_0_4 = pd.read_excel(file_0_4, sheet_name=None)
df_0_5 = pd.read_excel(file_0_5, sheet_name=None)

# Function to collect genes, corresponding pts values, and analysis sources

def collect_genes_from_source(sources):

    """
    Collect the required variables (genes, pts and the analysis of origin (0.3, 0.4 or 0.5).

    Parameters:
    sources (list): List of tuples where each element contains a dataframe and the corresponding analysis name.

    Returns:
    dict: Dictionary with genes as keys and a dictionary containing 'pts' and 'analysis' list as values.
    """

    all_genes = {}
    # Iterate through each source file (0.3, 0.4, 0.5)
    for source, name in sources:
        for cluster, df in source.items(): # Iterate through each cluster in the source file
            if cluster not in all_genes: #??
                all_genes[cluster] = {}

            # Iterate through rows in the cluster's DataFrame to collect gene data
            for _, row in df.iterrows():
                gene = row['Unnamed: 0']  # gene name, assumed to be the first column # why unnamed: 0
                pts = row['pts']  # the corresponding pts value

                # If gene is not already collected for this cluster, add it
                if gene not in all_genes[cluster]:
                    all_genes[cluster][gene] = {
                        'pts': round(pts, 2),
                        'analysis': [name]
                    }  # Round pts and store the analysis source
                else:
                    # Update pts if a higher value is found for the same gene
                    all_genes[cluster][gene]['pts'] = round(max(all_genes[cluster][gene]['pts'], pts), 2)
                    # Append the current analysis to the analysis list if not already included
                    if name not in all_genes[cluster][gene]['analysis']:
                        all_genes[cluster][gene]['analysis'].append(name)
    return all_genes

# Collect genes from all three sources and merge them into one dictionary
all_genes = collect_genes_from_source([(df_0_3, "0.3"), (df_0_4, "0.4"), (df_0_5, "0.5")])

# Export the merged data into a single Excel file (optional step, but kept for merging)
with pd.ExcelWriter(merged_output_file, engine='openpyxl') as writer:
    for cluster, genes in all_genes.items():
        cluster_df = pd.DataFrame({
            'Gene': list(genes.keys()),
            'pts': [gene_info['pts'] for gene_info in genes.values()],
            'pts of origin': [', '.join(gene_info['analysis']) for gene_info in genes.values()]
        })
        cluster_df.to_excel(writer, sheet_name=cluster, index=False)

print(f"Merged data successfully exported to {merged_output_file}")

# Now let's gather all the data into one sheet

# Load the merged Excel file (each sheet represents a cluster)
merged_data = pd.read_excel(merged_output_file, sheet_name=None)

# Initialize an empty DataFrame to store the combined data
combined_df = pd.DataFrame()

# Iterate through each cluster (sheet) in the merged Excel
for cluster, df in merged_data.items():
    # Add three columns for each cluster: one for Gene, one for pts, and one for 'pts of origin'
    combined_df[cluster] = df['Gene']
    combined_df[cluster + "_pts"] = df['pts'].round(2)  # Round pts values to 2 decimal places
    combined_df[cluster + "_pts of origin"] = df['pts of origin']  # Add the analysis source column

# Now rename the 'pts' and 'pts of origin' columns to remove the cluster name prefix
for col in combined_df.columns:
    if "_pts" in col or "_pts of origin" in col:
        new_col_name = col.split("_", 1)[-1]  # Remove everything before the first underscore
        combined_df.rename(columns={col: new_col_name}, inplace=True)

# Save the combined DataFrame into a new Excel file
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    combined_df.to_excel(writer, sheet_name="Combined_Genes", index=False)

print(f"Data successfully combined and saved to {output_file}")


