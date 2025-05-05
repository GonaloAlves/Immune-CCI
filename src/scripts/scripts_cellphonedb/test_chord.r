library("circlize")
library("readr")

# Load the CSV files
edge_list_15 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_15.csv")
edge_list_60 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_60.csv")

# Ensure 'value' is numeric
edge_list_15$value <- as.numeric(edge_list_15$value)
edge_list_60$value <- as.numeric(edge_list_60$value)

# Define biological cluster groups
color_groups <- list(
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

# Define color for each subgroup
color_map <- c(
  "Imm_Resting" = "#1f77b4", 
  "Imm_Other" = "#ff7f0e", 
  "Imm_Injury" = "#2ca02c", 
  "Neu" = "#d62728", 
  "MeV_Vascular" = "#9467bd", 
  "MeV_Epithelial" = "#8c564b", 
  "MeV_Fibroblast" = "#e377c2", 
  "MeV_Fib_Col" = "#7f7f7f", 
  "MeV_Fib_Lam" = "#bcbd22", 
  "VLMC" = "#17becf", 
  "Prolifs" = "#9e9e9e"
)

# Function to assign group based on the label
assign_group <- function(label) {
  for (group in names(color_groups)) {
    if (label %in% color_groups[[group]]) {
      return(group)
    }
  }
  return(NA)
}

# Assign groups based on the custom groupings
edge_list_15$group_from <- sapply(edge_list_15$from, assign_group)
edge_list_15$group_to <- sapply(edge_list_15$to, assign_group)

edge_list_60$group_from <- sapply(edge_list_60$from, assign_group)
edge_list_60$group_to <- sapply(edge_list_60$to, assign_group)

# Build color map for edge_list_15
all_labels_15 <- unique(c(edge_list_15$from, edge_list_15$to))

# Map each label to the appropriate color
grid.col_15 <- sapply(all_labels_15, function(label) {
  group <- assign_group(label)
  color_map[group]
})

# Build color map for edge_list_60
all_labels_60 <- unique(c(edge_list_60$from, edge_list_60$to))

# Map each label to the appropriate color
grid.col_60 <- sapply(all_labels_60, function(label) {
  group <- assign_group(label)
  color_map[group]
})

# =============================
# Plot for Injured_15
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_injured_15.pdf", width = 18, height = 18)

# Create the chord diagram for Injured_15 with correct group colors
chordDiagram(
  edge_list_15,
  grid.col = grid.col_15,
  transparency = 0.4,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = c("grid", "name", "axis"),
  preAllocateTracks = list(track.height = 0.1),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Highlight sectors (groups) and label them
highlight.sector(
  sector.index = c("Imm_Resting", "Imm_Other", "Imm_Injury", "Neu", "MeV_Vascular", "MeV_Epithelial", "MeV_Fibroblast", "MeV_Fib_Col", "MeV_Fib_Lam", "VLMC", "Prolifs"),
  track.index = 1,
  col = c(
    "#1f77b4", 
    "#ff7f0e", 
    "#2ca02c", 
    "#d62728", 
    "#9467bd", 
    "#8c564b", 
    "#e377c2", 
    "#7f7f7f", 
    "#bcbd22", 
    "#17becf", 
    "#9e9e9e"
  ),
  text = c(
    "Imm Resting", 
    "Imm Other", 
    "Imm Injury", 
    "Neu", 
    "MeV Vascular", 
    "MeV Epithelial", 
    "MeV Fibroblast", 
    "MeV Fib Collagen", 
    "MeV Fib Laminin", 
    "VLMC", 
    "Proliferative"
  ),
  cex = 0.8,
  text.col = "white", 
  niceFacing = TRUE
)

# Close the PDF output
dev.off()

# =============================
# Plot for Injured_60
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_injured_60.pdf", width = 18, height = 18)

# Create the chord diagram for Injured_60 with correct group colors
chordDiagram(
  edge_list_60,
  grid.col = grid.col_60,
  transparency = 0.4,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = c("grid", "name", "axis"),
  preAllocateTracks = list(track.height = 0.1),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Highlight sectors (groups) and label them
highlight.sector(
  sector.index = c("Imm_Resting", "Imm_Other", "Imm_Injury", "Neu", "MeV_Vascular", "MeV_Epithelial", "MeV_Fibroblast", "MeV_Fib_Col", "MeV_Fib_Lam", "VLMC", "Prolifs"),
  track.index = 1,
  col = c(
    "#1f77b4", 
    "#ff7f0e", 
    "#2ca02c", 
    "#d62728", 
    "#9467bd", 
    "#8c564b", 
    "#e377c2", 
    "#7f7f7f", 
    "#bcbd22", 
    "#17becf", 
    "#9e9e9e"
  ),
  text = c(
    "Imm Resting", 
    "Imm Other", 
    "Imm Injury", 
    "Neu", 
    "MeV Vascular", 
    "MeV Epithelial", 
    "MeV Fibroblast", 
    "MeV Fib Collagen", 
    "MeV Fib Laminin", 
    "VLMC", 
    "Proliferative"
  ),
  cex = 0.8,
  text.col = "white", 
  niceFacing = TRUE
)

# Close the PDF output
dev.off()
