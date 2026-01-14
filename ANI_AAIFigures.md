# Load Libraries
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
## Function to clean names
```
clean_names <- function(x) sub("\\.[0-9]$", "", x)
```

# Lyngbya 
```
tree = ggtree::read.tree("Lyngbya_astral.txt")

#add fake branch lengths
tree$edge.length <- rep(1, length(tree$edge.length))
tree$edge.length <- tree$edge.length * 2
# Example: set first branch to 1, second to 2, etc.
tree$edge.length <- seq_along(tree$edge.length)


p <- ggtree(tree, branch.length = "branch.length") +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, color = "blue")

tip_order <- p$data %>% filter(isTip) %>% pull(label)  # y increases from bottom to top

tip_order = tip_order %>% clean_names()


dendrogram <- ape::chronos(tree)
row_cluster <- as.hclust(dendrogram)
col_cluster <- as.hclust(dendrogram)

names_map <- read.table("Lyngbya_names.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(names_map) <- c("ID", "Name")


ani <- read.delim("Lyngbya_ANI.txt") %>%
  select(c(Ref_file, Query_file, ANI))%>%
  filter(Query_file %in% tip_order) %>%
  filter(Ref_file %in% tip_order) %>%
  pivot_wider(names_from = Query_file, values_from = ANI, values_fn = mean) %>%
  column_to_rownames(var = "Ref_file")

ordered_data_ANI <- ani[tip_order, tip_order]

aai <- read.delim("Lyngbya_AAI.txt") %>%
  select(Label.1, Label.2, AAI) %>%
  filter(Label.1 %in% tip_order) %>%
  filter(Label.2 %in% tip_order) %>%
  pivot_wider(names_from = Label.1, values_from = AAI, values_fn = mean) %>%
  column_to_rownames(var = "Label.2")

ordered_data_AAI <- aai[tip_order, tip_order]


combined_matrix <- as.matrix(ordered_data_ANI)
combined_matrix[upper.tri(combined_matrix)] <- as.matrix(ordered_data_AAI)[upper.tri(ordered_data_AAI)]

rownames(combined_matrix) <- names_map$Name[match(rownames(combined_matrix), names_map$ID)]
colnames(combined_matrix) <- names_map$Name[match(colnames(combined_matrix), names_map$ID)]

breaks <- c(77, 80, 85, 95, 100)

# Define group colors (must be length = length(breaks) - 1)
group_colors <- c("green", "yellow", "orange", "red")


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

## Inflection points
```
ani_raw <- read.table("Lyngbya_ANI.txt", header=TRUE, sep="\t", check.names=FALSE)
ani_raw$Ref_file <- names_map$Name[match(ani_raw$Ref_file, names_map$ID)]
ani_raw$Query_file <- names_map$Name[match(ani_raw$Query_file, names_map$ID)]

ani_raw %>%
  filter(Ref_file == "Capilliphycus GCA_000169095.1") %>%
  ggplot(aes(x = ANI, y = Align_fraction_ref, label = Query_file)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(x = "ANI", y = "Alignment Fraction (Ref)") +
  theme_minimal()+
  scale_x_continuous(limits = c(70, 100), breaks = seq(70, 100, by = 10)) +
  scale_y_continuous(limits = c(30, 100), breaks = seq(30, 100, by = 10))

ani_raw %>%
  filter(Ref_file == "Lyngbya confervoides BLCC-M349") %>%
  ggplot(aes(x = ANI, y = Align_fraction_ref, label = Query_file)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(x = "ANI", y = "Alignment Fraction (Ref)") +
  theme_minimal()+
  scale_x_continuous(limits = c(70, 100), breaks = seq(70, 100, by = 10)) +
  scale_y_continuous(limits = c(30, 100), breaks = seq(30, 100, by = 10))
```

# Okeania ANI
```
tree = ggtree::read.tree("Okeania_ASTRAL.txt")

#add fake branch lengths
tree$edge.length <- rep(1, length(tree$edge.length))
tree$edge.length <- tree$edge.length * 2
# Example: set first branch to 1, second to 2, etc.
tree$edge.length <- seq_along(tree$edge.length)


p <- ggtree(tree, branch.length = "branch.length") +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, color = "blue")

tip_order <- p$data %>% filter(isTip) %>% pull(label)  # y increases from bottom to top

tip_order = tip_order %>% clean_names()

dendrogram <- ape::chronos(tree)
row_cluster <- as.hclust(dendrogram)
col_cluster <- as.hclust(dendrogram)

names_map <- read.table("Okeania_ANI_names.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(names_map) <- c("ID", "Name")

ani <- read.delim("Okeania_ANI.txt") %>%
  select(c(Ref_file, Query_file, ANI))%>%
  filter(Query_file %in% tip_order) %>%
  filter(Ref_file %in% tip_order) %>%
  pivot_wider(names_from = Query_file, values_from = ANI, values_fn = mean) %>%
  column_to_rownames(var = "Ref_file")

ordered_data_ANI <- ani[tip_order, tip_order]

rownames(ordered_data_ANI) <- names_map$Name[match(rownames(ordered_data_ANI), names_map$ID)]
colnames(ordered_data_ANI) <- names_map$Name[match(colnames(ordered_data_ANI), names_map$ID)]

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
