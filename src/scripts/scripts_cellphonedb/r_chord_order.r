library("circlize")
library("readr")

# Load the CSV files
edge_list_15 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_15.csv")
edge_list_60 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_60.csv")

# Ensure 'value' is numeric
edge_list_15$value <- as.numeric(edge_list_15$value)
edge_list_60$value <- as.numeric(edge_list_60$value)

cat("Unique labels in edge list:\n")
print(edge_list_15)
print(edge_list_60)

# Define biological cluster groups
biological_groups <- list(
  Imm_Resting = c("Imm.M0Like.0", "Imm.M0Like.1", "Imm.M0Like.2"), 
  Imm_Other = c("Imm.MHCII.0", "Imm.PVM.0"), 
  Imm_Injury = c("Imm.Interferon.0", "Imm.DAM.0", "Imm.DAM.1"),
  Neu = c("Neu.CSFcN.0", "Neu.Epend.0"),
  MeV_Vascular = c("MeV.Endothelial.0", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.Endothelial.3", "MeV.Pericytes.0", "MeV.SMC.0"),
  MeV_Epithelial = c("MeV.Epithelial.0"),
  MeV_Fibroblast = c("MeV.Fib.0", "MeV.Fib.1", "MeV.Fib.2", "MeV.Fib.3", "MeV.Fib.4", "MeV.Fib.5"),
  MeV_Fib_Col = c("MeV.FibCollagen.0", "MeV.FibCollagen.1", "MeV.FibCollagen.2", "MeV.FibCollagen.3"),
  MeV_Fib_Lam = c("MeV.FibLaminin.0"),
  VLMC = c("MeV.VLMC.0", "MeV.VLMC.1"),
  Prolifs = c("Imm.Proliferative.0", "MeV.FibProlif.0")
)

cat("Biological groups:\n")
print(biological_groups)

# Flatten biological_groups to a lookup table
cluster_to_group <- unlist(biological_groups)
names(cluster_to_group) <- rep(names(biological_groups), lengths(biological_groups))
cat("briuh:\n")
print(cluster_to_group)

# Assign colors to biological groups
group_palette <- rainbow(length(biological_groups))
names(group_palette) <- names(biological_groups)

# Function to get biological group
get_bio_group <- function(label) {
  return(cluster_to_group[label])
}

# =============================
# Function to Plot Chord Diagram
# =============================
plot_chord <- function(edge_list, file_path) {


  circos.clear()

  cat("\n======== Starting plot for:", file_path, "========\n")

  # Get all unique labels
  all_labels <- unique(c(edge_list$from, edge_list$to))
  cat("Total unique labels from 'from' and 'to':", length(all_labels), "\n")
  cat("Unique labels in edge list:\n")
  print(all_labels)

  # Assign biological groups
  group_assignments <- sapply(all_labels, get_bio_group)
  cat("Group assignments (non-NA):", sum(!is.na(group_assignments)), "\n")
  cat("Group assignments (NA):", sum(is.na(group_assignments)), "\n")
  cat("Mapped biological groups:\n")
  print(group_assignments)

  # Filter only valid labels
  valid_labels <- all_labels[!is.na(group_assignments)]
  group_assignments <- group_assignments[!is.na(group_assignments)]

  cat("Valid labels:\n")
  print(valid_labels)

  cat("Number of valid labels after filtering:", length(valid_labels), "\n")

  if (length(valid_labels) == 0) {
    warning("No valid labels found, skipping plot.")
    return(NULL)
  }

  # Filter edge list
  edge_list <- edge_list[edge_list$from %in% valid_labels & edge_list$to %in% valid_labels, ]
  cat("Number of edges after filtering:", nrow(edge_list), "\n")

  if (nrow(edge_list) == 0) {
    warning(paste("No valid edges for:", file_path))
    return(NULL)
  }

  # Assign colors
  grid_colors <- setNames(group_palette[group_assignments], valid_labels)

  # Check color mapping
  cat("Assigned colors to sectors:\n")
  print(head(grid_colors))

  # Start PDF output
  pdf(file_path, width = 10, height = 10)

  chordDiagram(
    edge_list,
    grid.col = grid_colors,
    transparency = 0.4,
    directional = 1,
    direction.type = c("diffHeight", "arrows"),
    annotationTrack = "grid",
    link.arr.type = "big.arrow",
    link.sort = TRUE,
    link.decreasing = FALSE,
    preAllocateTracks = 1
  )

  # Add sector labels
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(5),
                  sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
    }, bg.border = NA
  )

  # Add biological group legend
  legend("topleft", legend = names(group_palette), fill = group_palette, cex = 0.7, bty = "n")

  dev.off()

  cat("======== Plot complete for:", file_path, "========\n")
}

# =============================
# Run for Injured_15 and Injured_60
# =============================
plot_chord(edge_list_15, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_15_order.png")
plot_chord(edge_list_60, "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_60_order.png")
