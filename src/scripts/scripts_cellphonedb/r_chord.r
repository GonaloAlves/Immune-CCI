# install.packages("migest")
# install.packages("circlize")

library("circlize")
library("readr")

# Load the CSV file
edge_list_15 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list_injured_15_copy.csv")
edge_list_60 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list_injured_60.csv")


# Make sure 'value' is numeric
edge_list_15$value <- as.numeric(edge_list_15$value)


output_file <- "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/chord_diagram_15.pdf"

# Open PDF device
pdf(output_file, width = 10, height = 10)


# Clear previous plots
circos.clear()

# Example: define color groups by regex
group_colors <- setNames(
  c("red", "blue", "green"),
  c("Imm", "MeV", "Neu")
)

# Function to match prefix
get_group <- function(label) sub("\\..*$", "", label)

# Build color map
all_labels <- unique(c(edge_list_15$from, edge_list_15$to))
grid.col <- setNames(
  group_colors[get_group(all_labels)],
  all_labels
)


# Chord diagram
chordDiagram(
  edge_list_15,
  grid.col = grid.col,  
  transparency = 0.4,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = "grid",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Close the device and write the file
dev.off()
