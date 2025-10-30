# Load libraries
```
library(devtools)
library(tidyverse)
library(ggtree)
library(ComplexHeatmap)
library(circlize)
library(ggtreeExtra)
library(treeio)
library(ggrepel)
```

# Creat function to clean names
```
clean_names <- function(x) sub("\\.[0-9]$", "", x)
```
# Lyngbya ANI/AAI
## Load the phylogenomic tree
```
tree = ggtree::read.tree("Lyngbya_astral.txt")
```
## Add fake branch lengths
This makes plotting easier 
```
tree$edge.length <- rep(1, length(tree$edge.length))
tree$edge.length <- tree$edge.length * 2
tree$edge.length <- seq_along(tree$edge.length)
```

## Plot tree
```
p <- ggtree(tree, branch.length = "branch.length") +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, color = "blue")
```

## Set tip order to organize heatmaps
```
tip_order <- p$data %>% filter(isTip) %>% pull(label)  # y increases from bottom to top
tip_order = tip_order %>% clean_names()
```
## Create dendrograms
```
dendrogram <- ape::chronos(tree)
row_cluster <- as.hclust(dendrogram)
col_cluster <- as.hclust(dendrogram)
```

## Load names for heatmap
```
names_map <- read.table("names2.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(names_map) <- c("ID", "Name")
```

## Load ANI values and create matrix
```
ani <- read.delim("ANI2.txt") %>%
  select(c(Ref_file, Query_file, ANI))%>%
  filter(Query_file %in% tip_order) %>%
  filter(Ref_file %in% tip_order) %>%
  pivot_wider(names_from = Query_file, values_from = ANI, values_fn = mean) %>%
  column_to_rownames(var = "Ref_file")
```
## Order matix
```
ordered_data_ANI <- ani[tip_order, tip_order]
```
## Load AAI values and create matrix
```
aai <- read.delim("~/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/MANUSCRIPTS/Lyngbya_confervoides/Analyses/Oscillatoriales/AAI2.txt") %>%
  select(Label.1, Label.2, AAI) %>%
  filter(Label.1 %in% tip_order) %>%
  filter(Label.2 %in% tip_order) %>%
  pivot_wider(names_from = Label.1, values_from = AAI, values_fn = mean) %>%
  column_to_rownames(var = "Label.2")
```
## Order matix
```
ordered_data_AAI <- aai[tip_order, tip_order]
```

## Combine the ANI and AAI matrices
```
combined_matrix <- as.matrix(ordered_data_ANI)
combined_matrix[upper.tri(combined_matrix)] <- as.matrix(ordered_data_AAI)[upper.tri(ordered_data_AAI)]
```
## Replace names
```
rownames(combined_matrix) <- names_map$Name[match(rownames(combined_matrix), names_map$ID)]
colnames(combined_matrix) <- names_map$Name[match(colnames(combined_matrix), names_map$ID)]
```
## Create breaks and define colors
```
breaks <- c(77, 80, 85, 95, 100)
group_colors <- c("green", "yellow", "orange", "red")
```

## Plot ANI/AAI heatmap
```
ComplexHeatmap::pheatmap(
  combined_matrix,
  legend_breaks = (c(80, 85, 95, 100)),
  #color = group_colors,
  #breaks = breaks,
  display_numbers = T,
  number_color = "black",
  cluster_rows = row_cluster,
  cluster_cols = col_cluster,
  fontsize_col = 8,
  fontsize_number = 6,
  name = "ANI/AAI",
  angle_col = "45",
  number_format = "%.1f",
  column_title = "ANI (lower left) / AAI (upper right)"
)
```

# OKEANIA ANI
## Load the phylogenomic tree
```
tree = ggtree::read.tree("Okeania_ASTRAL.txt")
```
## Add fake branch lengths
This makes plotting easier 
```
tree$edge.length <- rep(1, length(tree$edge.length))
tree$edge.length <- tree$edge.length * 2
tree$edge.length <- seq_along(tree$edge.length)
```
## Plot tree
```
p <- ggtree(tree, branch.length = "branch.length") +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, color = "blue")
```
## Set tip order to organize heatmaps
```
tip_order <- p$data %>% filter(isTip) %>% pull(label)  # y increases from bottom to top
tip_order = tip_order %>% clean_names()
```
## Create dendrograms
```
dendrogram <- ape::chronos(tree)
row_cluster <- as.hclust(dendrogram)
col_cluster <- as.hclust(dendrogram)
```
## Load names for heatmap
```
names_map <- read.table("names.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(names_map) <- c("ID", "Name")
```
## Load ANI values and create matrix
```
ani <- read.delim("ANI.txt") %>%
  select(c(Ref_file, Query_file, ANI))%>%
  filter(Query_file %in% tip_order) %>%
  filter(Ref_file %in% tip_order) %>%
  pivot_wider(names_from = Query_file, values_from = ANI, values_fn = mean) %>%
  column_to_rownames(var = "Ref_file")
```
## Order the matrix
```
ordered_data_ANI <- ani[tip_order, tip_order]
```
## Replace names
```
rownames(ordered_data_ANI) <- names_map$Name[match(rownames(ordered_data_ANI), names_map$ID)]
colnames(ordered_data_ANI) <- names_map$Name[match(colnames(ordered_data_ANI), names_map$ID)]
```

## Plot the matrix
```
ComplexHeatmap::pheatmap(
  as.matrix(ordered_data_ANI),
  legend_breaks = (c(80, 85, 95, 100)),
  display_numbers = T,
  number_color = "black",
  cluster_rows = row_cluster,
  cluster_cols = col_cluster,
  fontsize_col = 8,
  fontsize_number = 6,
  name = "ANI",
  angle_col = "45",
  number_format = "%.1f",
  column_title = "ANI"
)
```
