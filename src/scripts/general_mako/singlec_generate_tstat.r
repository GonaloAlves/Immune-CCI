# Output moderate t-statistic tables for GSEA - LeonorSaude10x
#
# Use limma to generate moderate t-statistics values for each subpopulation and cluster within,
# idellay *_raw_norm_rank.h5ad stage of analysis shoud be used
# Requires to have cluster information generated in the .h5ad file
# The input for limma is log-normalized X matrix 
#
# Daniel Ribeiro, 2023, ft Gonçalo Alves 2024

library("limma")
library("Matrix")
library("anndata")
library("R.utils")

# Set WD
argv = commandArgs(trailingOnly = FALSE)
base_dir = dirname(substring(argv[grep("--file=", argv)], 8))
setwd(base_dir)
print(getwd())
options(Ncpus = parallel::detectCores())

# Dataset suffix
output_dir = "/home/makowlg/Documents/Immune-CCI/src"

# Define dirs
checkpoint_dir = "/home/makowlg/Documents/Immune-CCI/h5ad_files/"
tstat_dir = sprintf("%s/tstat_files_mako", output_dir)

# Create the directory for t-stat files if it doesn't exist
if (!dir.exists(tstat_dir)) {
    dir.create(tstat_dir, recursive = TRUE)
    print(sprintf("Directory %s created.", tstat_dir))
} else {
    print(sprintf("Directory %s already exists.", tstat_dir))
}

# Dataset
d = "Meningeal_Vascular"
resolution = "leiden_mako"

print("Export moderated t-statistic values...")

# Load data
print("Load expression data...")
start_time = Sys.time()     # time execution
dest = sprintf("%s/adata_final_%s_raw_norm_ranked_copy_copy.h5ad", checkpoint_dir, d)
print(dest)

ann = read_h5ad(dest)
print("H5AD file contents:")
print(ann)

### Cluster analysis ###
clusters = sprintf("leiden_mako")
printf("Analysing cluster %s...\n", clusters)
printf("Clusters: %s \n", paste(levels(ann$obs[[clusters]]), collapse = " "))
for (c in levels(ann$obs[[clusters]]))
{
  printf("\nt-stat for cluster %s\n", c)
  printf("Cells in test: %i, Cells in rest: %i\n",  sum(ann$obs[[clusters]] == c), sum(ann$obs[[clusters]] != c))
      
  # Comparison groups. group0 = test group, group1 = rest of cells of other clusters
  group = rep(0, dim(ann)[1])
  group[which(ann$obs[[clusters]] == c)] = 0
  group[which(ann$obs[[clusters]] != c)] = 1
  group = factor(group)
  
  # Model matrix
  mm = model.matrix(~0 + group)
  
  # Fitting linear models in limma #
  # lmFit fits a linear model using weighted least squares
  y = t(as.data.frame(ann))   # Convert to df and transpose to be R-style: rows=genes, cols=cells
  fit = lmFit(y, mm)
  
  # Comparisons between groups (Log2FC) are defined by contrasts
  # Test contrasts Rest
  contrast = makeContrasts(group0 - group1, levels = colnames(coef(fit)))
  # Estimate contrasts for each gene
  tmp = contrasts.fit(fit, contrast)
  # Perform empirical Bayes smoothing of standard errors
  tmp = eBayes(tmp)
  # See most significant genes
  top_table = topTable(tmp, sort.by = "P", n = Inf)
  
  # Add Gene names and Ensembl IDs
  top_table["Gene.Name"] = rownames(top_table)
  top_table["Ensembl.ID"] = ann$var[["ensembl_gene_id"]][match(rownames(top_table), colnames(ann))]
  top_table = top_table[, c("Gene.Name", "Ensembl.ID", colnames(top_table)[1:(length(colnames(top_table)) - 2)])] 
  # Calculate how many genes are DE by adj.P.Val (<0.05)
  deg = length(which(top_table['adj.P.Val'] < 0.05))
  printf("DEGs: %i\n", deg)
  
  # Write to file #
  printf("Write table...\n")
  dest = sprintf("%s/%s_%s_c%s_tstat.xlsx", tstat_dir, d, clusters, c)
  write.table(top_table, file = dest, row.names = F, sep = "\t", quote = F)
  printf("%s\n", dest)

  rm(top_table, tmp, y, fit, contrast)
  gc()

}
  
  # Free memory
  rm(ann)
  gc()
#}

printf("\n********\n* DONE *\n********\n")
