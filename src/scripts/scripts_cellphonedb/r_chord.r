library("circlize")
library("readr")
library("magrittr")

# Load the CSV files
edge_list_15 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_15.csv")
edge_list_60 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/grouped_edge_list_60.csv")

# Ensure 'value' is numeric
edge_list_15$value <- as.numeric(edge_list_15$value)
edge_list_60$value <- as.numeric(edge_list_60$value)


# Function to extract group prefix
#get_group <- function(label) sub("\\..*$", "", label)

get_group <- function(label) {
  # Match everything up to the second dot (if exists), else return entire string
  m <- regexpr("^([^\\.]+\\.[^\\.]+)", label)
  ifelse(m > 0, regmatches(label, m), label)
}

# Define group colors
# group_colors <- setNames(
#   c("red", "blue", "green"),
#   c("Imm", "MeV", "Neu")
# )

get_main_group <- function(label) sub("\\..*$", "", label)
group_colors <- setNames(
  c("blue", "red", "green")[match(get_main_group(unique(sapply(c(edge_list_15$from, edge_list_15$to), get_group))), c("MeV", "Imm", "Neu"))],
  unique(sapply(c(edge_list_15$from, edge_list_15$to), get_group))
)




# Function to highlight sectors by group
highlight_groups <- function(labels, track_index = 1) {
  unique_groups <- unique(sapply(labels, get_group))
  for (group in unique_groups) {
    sector_labels <- labels[get_group(labels) == group]
    highlight.sector(
      sector_labels,
      track.index = track_index,
      col = group_colors[group],
      text = group,
      cex = 0.7,
      text.col = "black",
      niceFacing = TRUE
    )
  }
}

# # Add external group labels
# add_group_labels <- function(labels, group_colors) {
#   unique_groups <- unique(sapply(labels, get_group))
#   for (group in unique_groups) {
#     # Get all sectors for this group
#     sectors <- labels[get_group(labels) == group]
#     if (length(sectors) == 0) next

#     # Choose the first sector as representative (simplified approach)
#     sector <- sectors[1]
#     x_pos <- get.cell.meta.data("xcenter", sector.index = sector, track.index = 1)

#     # Add the text outside the sector
#     circos.text(
#       x = x_pos,
#       y = 1.5,
#       labels = group,
#       sector.index = sector,
#       track.index = 1,
#       facing = "bending.outside",
#       niceFacing = TRUE,
#       col = group_colors[group],
#       cex = 1.2
#     )
#   }
# }


# =============================
# Plot for Injured_15
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/order_chord_diagram_15__.pdf", width = 18, height = 18)
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
  annotationTrack = c("grid", "name", "axis"),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Add group labels
highlight_groups(all_labels_15)


dev.off()

# =============================
# Plot for Injured_60
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/order_chord_diagram_60__.pdf", width = 18, height = 18)
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
  annotationTrack = c("grid", "name", "axis"),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Add group labels
highlight_groups(all_labels_60)
#add_group_labels(all_labels_60, group_colors)

dev.off()
