library(tidyverse)
library("here")
library(devtools)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(scales) 
library(patchwork)


# Generate Louvain Clustering Results  ------------------------------------

umap_cluster_colors <- c("#DEA940", "#A2A763", "#D6EEBE", "#5894D4", "#6AAD51", "#6C8686", "#FFD92F", "#7f2201", "#b53310",
                         "#E7673B", "#FB9A99", "#E31A1C", "#FDBF6F", "#CAB2D6", "#386CB0", "#FFFFB3", "#5dead5", "#c11577",
                         "#CCCCCC", "#7570B3", "#B15928", "#6A3D9A", "#d87c6a", "#0060e8", "#84c5ff", "#ffd460", "#d81397",
                         "#4d908e", "#b088f9", "#800020", "#30D5C8")


uf_annotation_final <- read.delim("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/uf_annot/uf.leaf_annotation.V1.meta.txt") %>% 
  mutate(across(c(LouvainClusters),factor))


colnames(uf_annotation_final)
colnames(uf_annotation_final)

library(ggplot2)
# Assuming uf_annotation_final and umap_cluster_colors are already defined

# Plot for 'tss'
lc_tss <- ggplot(uf_annotation_final, aes(x = LouvainClusters, y = log(tss), fill = LouvainClusters)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'total'
lc_total <- ggplot(uf_annotation_final, aes(x = LouvainClusters, y = log(total), fill = LouvainClusters)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'FRiP'
lc_frip <- ggplot(uf_annotation_final, aes(x = LouvainClusters, y = FRiP, fill = LouvainClusters)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'ptmt_ratio'
lc_ptmt_ratio <- ggplot(uf_annotation_final, aes(x = LouvainClusters, y = (ptmt/total), fill = LouvainClusters)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, .5)) # Extend y-axis to include zero

# Plot for 'log10nSites'
lc_lognsites <- ggplot(uf_annotation_final, aes(x = LouvainClusters, y = log10nSites, fill = LouvainClusters)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = umap_cluster_colors)

# Plot for 'doublet_score'
lc_doublet_score <- ggplot(uf_annotation_final, aes(x = LouvainClusters, y = doubletscore, fill = LouvainClusters)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1.5)) # Extend y-axis to include zero


colnames(uf_annotation_final_count)
# Stacked + percent
uf_annotation_final_count <- uf_annotation_final %>% 
  mutate(count = 1)

lc_cell_prop <- ggplot(uf_annotation_final_count, aes(fill=LouvainClusters, y=count, x=sampleID)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = umap_cluster_colors)  +
  theme_cowplot() +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero


lc_metrics <- (lc_tss + lc_total)/(lc_frip + lc_ptmt_ratio)/(lc_lognsites + lc_doublet_score)
lc_metrics

ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Uf/Uf.louvain_metrics.pdf",
       lc_metrics,
       width = 22, 
       height = 15,
       units = "in")


ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Uf/Uf.loucain_cell_prop.pdf",
       lc_cell_prop,
       width = 5, 
       height = 8,
       units = "in")


########### Same Graphs but for Annotations
#######
library(ggplot2)
library(dplyr)

# Assuming uf_annotation_final and umap_cluster_colors are already defined
colnames(uf_annotation_final)
uf_annotation_final <- uf_annotation_final %>%
  mutate(annotation_v1 = case_when(
    annotation_v1 == "companion_cells_sieve_elements" ~ "CC/SE",
    TRUE ~ annotation_v1 # Keeps all other values as they are
  ))
unique(uf_annotation_final$annotation_v1)


# Plot for 'tss'
anno_tss <- ggplot(uf_annotation_final, aes(x = annotation_v1, y = log(tss), fill = annotation_v1)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'total'
anno_total <- ggplot(uf_annotation_final, aes(x = annotation_v1, y = log(total), fill = annotation_v1)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'FRiP'
anno_frip <- ggplot(uf_annotation_final, aes(x = annotation_v1, y = FRiP, fill = annotation_v1)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'ptmt_ratio'
anno_ptmt_ratio <- ggplot(uf_annotation_final, aes(x = annotation_v1, y = (ptmt/total), fill = annotation_v1)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, .5)) # Extend y-axis to include zero

# Plot for 'log10nSites'
anno_lognsites <- ggplot(uf_annotation_final, aes(x = annotation_v1, y = log10nSites, fill = annotation_v1)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA))

# Plot for 'doubletscore'
anno_doublet_score <- ggplot(uf_annotation_final, aes(x = annotation_v1, y = doubletscore, fill = annotation_v1)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1.5)) # Extend y-axis to include zero


colnames(uf_annotation_final_count)
# Stacked + percent
uf_annotation_final_count <- uf_annotation_final %>% 
  mutate(count = 1)

anno_cell_prop <- ggplot(uf_annotation_final_count, aes(fill=annotation_v1, y=count, x=sampleID)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = umap_cluster_colors)  +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Extend y-axis to include zero


unique(uf_annotation_final_count$annotation_v1)
anno_metrics <- (anno_tss + anno_total)/(anno_frip + anno_ptmt_ratio)/(  anno_lognsites + anno_doublet_score)
anno_metrics


ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Uf/Uf.anno_metrics.pdf",
       anno_metrics,
       width = 10, 
       height = 15,
       units = "in")

ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Uf/Uf.anno_cell_prop.pdf",
       anno_cell_prop,
       width = 5, 
       height = 8,
       units = "in")


