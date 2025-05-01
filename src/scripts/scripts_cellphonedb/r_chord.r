library("circlize")
library("readr")

# Load the CSV files
edge_list_15 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list_injured_15.csv")
edge_list_60 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list_injured_60_copy.csv")

# Ensure 'value' is numeric
edge_list_15$value <- as.numeric(edge_list_15$value)
edge_list_60$value <- as.numeric(edge_list_60$value)

# Define group colors
group_colors <- setNames(
  c("red", "blue", "green"),
  c("Imm", "MeV", "Neu")
)

# Function to extract group prefix
get_group <- function(label) sub("\\..*$", "", label)

# =============================
# Plot for Injured_15
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_15.pdf", width = 10, height = 10)
circos.clear()

# Build color map
all_labels_15 <- unique(c(edge_list_15$from, edge_list_15$to))
grid.col_15 <- setNames(
  group_colors[get_group(all_labels_15)],
  all_labels_15
)

# Draw chord diagram
chordDiagram(
  edge_list_15,
  grid.col = grid.col_15,
  transparency = 0.4,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = "grid",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)
dev.off()

# =============================
# Plot for Injured_60
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_60_10lig.pdf", width = 10, height = 10)
circos.clear()

# Build color map
all_labels_60 <- unique(c(edge_list_60$from, edge_list_60$to))
grid.col_60 <- setNames(
  group_colors[get_group(all_labels_60)],
  all_labels_60
)

# Draw chord diagram
chordDiagram(
  edge_list_60,
  grid.col = grid.col_60,
  transparency = 0.4,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = "grid",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)
dev.off()
