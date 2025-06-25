library("circlize")
library("readr")
library("magrittr")

# Load the CSV files
edge_list_15 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/filtered_pvalues/grouped_edge_list_15.csv")
edge_list_60 <- read_csv("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/filtered_pvalues/grouped_edge_list_60.csv")

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

imm_red_colors <- c(  #immune verdes
  "#91cb3e33", # 33
  "#00000000", 
  "#5bba6f33", #
  "#73a94233", #
  "#11982233", #
  "#00000000", 
  "#00000000", 
  "#00000000",
  "#00000000", 
  "#00000000"  
)
mev_blue_colors <- c( # mevs azuis
  "#7b70d233", #7b70d2ff
  "#9c70c233", #9c70c2
  "#00000000",
  "#00000000",
  "#729fcf33", # 
  "#3f7cac33", #
  "#00000000",
  "#00000000",
  "#00000000"  
)
neu_green_colors <- c("#ff666bff", #neus red/purple ff666b
 "#00000000")


imma_red_colors <- c(  #immune verdes
  "#9c70c2", #
  "#729fcf", 
  "#7b70d2", #
  "#729fcf", #
  "#486ec7", #
  
  
  "#f57c73", #f57c73
  "#f68c70", #aa6554
  "#f6ac69",
  "#f6bc66", 
  "#f6bc66"  
)
meva_blue_colors <- c( # mevs azuis
  "#f6ac69", #
  "#ff666b", #

  "#7ec4cf",
  "#9cadce",

  "#fe9e9e", # 
  "#c24b72", #

  "#9898ab",
  "#8494c0",
  "#b2b2a4"  
)
neua_green_colors <- c("#00bf63", #neus red/purple
 "#bbdb44")


mev_groups <- c("MeV.Endothelial", "MeV.Pericytes", "MeV.SMC", "MeV.Epithelial", "MeV.Fib", "MeV.FibCollagen", "MeV.FibLaminin", "MeV.VLMC", "MeV.FibProlif")
imm_groups <- c("Imm.Homeostatic", "Imm.MHCII", "Imm.PVM", "Imm.Interferon", "Imm.DAM", "Imm.Proliferative")
neu_groups <- c("Neu.Epend", "Neu.CSFcN")



get_main_group <- function(label) sub("\\..*$", "", label)
# group_colors <- setNames(
#   c("blue", "red", "green")[match(get_main_group(unique(sapply(c(edge_list_15$from, edge_list_15$to), get_group))), c("MeV", "Imm", "Neu"))],
#   unique(sapply(c(edge_list_15$from, edge_list_15$to), get_group))
# )

group_colors <- setNames(
  c(mev_blue_colors[1:length(mev_groups)],
    imm_red_colors[1:length(imm_groups)],
    neu_green_colors[1:length(neu_groups)]),
  c(mev_groups, imm_groups, neu_groups)
)


# Function to highlight sectors by group, without text labels
highlight_groups <- function(labels, track_index = 1) {
  unique_groups <- unique(sapply(labels, get_group))
  for (group in unique_groups) {
    sector_labels <- labels[get_group(labels) == group]
    highlight.sector(
      sector_labels,
      track.index = track_index,
      col = group_colors[group],
      text = NULL,           # ← disables the label
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
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/sigs/order_chord_diagram_15.pdf", width = 18, height = 18)
circos.clear()

# Build color map
all_labels_15 <- unique(c(edge_list_15$from, edge_list_15$to))
grid.col_15 <- setNames(
  group_colors[sapply(all_labels_15, get_group)],
  all_labels_15
)


# Draw chord diagram
chordDiagram(
  edge_list_15,
  grid.col = grid.col_15,
  transparency = 0,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = "grid",
  annotationTrackHeight = c(0.02, 0.03),  # ← control size of grid + name tracks
  preAllocateTracks = list(track.height = 0.03),  # Optional: fine-tune sector track
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Add names manually with full control
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  circos.text(
    x = CELL_META$xcenter,
    y = CELL_META$ycenter + 1.5,  # ← Move name outward
    labels = sector_name,
    facing = "bending",
    niceFacing = TRUE,
    cex = 1.7,                    # ← Font size
    font = 2                     # ← Bold
  )
}, bg.border = NA)

# Add group labels
highlight_groups(all_labels_15)


dev.off()

# =============================
# Plot for Injured_60
# =============================
pdf("/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/chord/sigs/order_chord_diagram_60.pdf", width = 18, height = 18)
circos.clear()

# Build color map
all_labels_60 <- unique(c(edge_list_60$from, edge_list_60$to))
grid.col_60 <- setNames(
  group_colors[sapply(all_labels_60, get_group)],
  all_labels_60
)


# Draw chord diagram
chordDiagram(
  edge_list_60,
  grid.col = grid.col_60,
  transparency = 0,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  annotationTrack = "grid",
  annotationTrackHeight = c(0.02, 0.03),  # ← control size of grid + name tracks
  preAllocateTracks = list(track.height = 0.03),  # Optional: fine-tune sector track
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = FALSE
)

# Add names manually with full control
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector_name <- CELL_META$sector.index
  circos.text(
    x = CELL_META$xcenter,
    y = CELL_META$ycenter + 1.5,  # ← Move name outward
    labels = sector_name,
    facing = "bending",
    niceFacing = TRUE,
    cex = 1.25,  # ← Font size
    font = 2   # ← Bold
  )
}, bg.border = NA)

# Add group labels
highlight_groups(all_labels_60)
#add_group_labels(all_labels_60, group_colors)

dev.off()
