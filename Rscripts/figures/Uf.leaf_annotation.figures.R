library(tidyverse)
library("here")
library(devtools)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)



uf_annotation_final <- read.delim("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/uf_annot/uf.leaf_annotation.V1.meta.txt") %>% 
  mutate(across(c(LouvainClusters),factor))


colnames(uf_annotation_final)

# Generate Louvain Clustering Results  ------------------------------------

umap_cluster_colors <- c("#DEA940","#A2A763","#D6EEBE","#5894D4","#6AAD51","#6C8686","#FFD92F", "#7f2201","#b53310",
                         "#E7673B","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6","#386CB0","#FFFFB3", "#5dead5","#c11577",
                         "#CCCCCC", "#7570B3","#B15928", "#6A3D9A", "#d87c6a", "#0060e8", "#84c5ff", "#ffd460", "#d81397")


final_annotation_set <- unlist(unique(uf_annotation_final$LouvainClusters))
names(umap_cluster_colors)  <- levels(final_annotation_set)
col_scale <- scale_colour_manual(name = "grp", values = umap_cluster_colors)


uf.louvain_clusters <- ggplot(uf_annotation_final, aes(x=umap1, y = umap2, color = LouvainClusters)) +
  theme_half_open() +
  col_scale + geom_point(size = .25, alpha = .8) + ggtitle("sb Cell Annotation Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) + ggtitle("Zea mays Louvain Clusters")

uf.louvain_clusters


ggsave("uf.louvain_clusters.pdf", 
       plot = uf.louvain_clusters, 
       path = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 10, height = 10)

ggsave("uf.louvain_clusters.png", 
       plot = uf.louvain_clusters, 
       path = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 10, height = 10)


# Generate final UMAP for Figure 1 ----------------------------------------
umap_cluster_colors <- c("#DEA940","#A2A763","#D6EEBE","#5894D4","#6AAD51","#6C8686","#FFD92F", "#7f2201","#b53310",
                         "#E7673B","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6","#386CB0","#FFFFB3", "#5dead5","#c11577",
                         "#CCCCCC", "#7570B3","#B15928", "#6A3D9A", "#d87c6a", "#0060e8", "#84c5ff", "#ffd460", "#d81397")


hard_code_annotation_colors <- c(
  "epidermis" = "#A2A763",
  "protoderm" = "#D6EEBE",
  "unknown" = "#CAB2D6",
  "companion cells/sieve elements" = "#CCCCCC",
  "mesophyll" = "#E31A1C",
  "mesophyll developing" = "#FFD92F",
  "phloem_sieve_elements" = "#ffd460",
  "bundle sheath" = "#5894D4",
  "procambial meristem" = "#B15928")

colnames(uf_annotation_final)
uf_annotation_final <- uf_annotation_final %>% 
  mutate(final_annotation = annotation_v1) %>% 
  mutate(cell_type = case_when(final_annotation == "phloem_sieve_elements" ~ "companion_cells/sieve_elements",
                               final_annotation == "companion_cells" ~ "companion_cells/sieve_elements",
                               final_annotation == "companion_cells_sieve_elements" ~ "companion_cells/sieve_elements",
                               final_annotation == "unknown;1" ~ "unknown",
                               final_annotation == "procambium" ~ "procambial meristem",
                               final_annotation == "mesophyll;developing" ~ "mesophyll developing",
                               TRUE ~ final_annotation)) %>% 
  dplyr::filter(is.na(cell_type) != TRUE)  %>%  
  mutate(cell_type = str_replace_all(cell_type, "_", " "))  %>% 
  mutate_at(vars(cell_type), 
            list(factor))


unique(uf_annotation_final$cell_type)

final_annotation_set <- unlist(unique(uf_annotation_final$cell_type))
names(umap_cluster_colors)  <- levels(final_annotation_set)
col_scale <- scale_colour_manual(name = "grp", values = umap_cluster_colors)


fig_1_umap_pdf <- ggplot(uf_annotation_final, aes(x = umap1, y = umap2, color = cell_type)) +
  theme_half_open() +
  geom_point(size = 0.25, alpha = 0.8) +
  ggtitle("Urochloa fusca cell type annotation") +
  scale_color_manual(values = hard_code_annotation_colors) +
  guides(colour = guide_legend(override.aes = list(size = 5), title = "Cell type")) +
  theme(legend.title.align = 0.5)
fig_1_umap_pdf

ggsave("uf.cell_annot.pdf", 
       plot = fig_1_umap_pdf, 
       path = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 8, height = 5)

ggsave("uf.cell_annot.png", 
       plot = fig_1_umap_pdf, 
       path = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 10, height = 10)
