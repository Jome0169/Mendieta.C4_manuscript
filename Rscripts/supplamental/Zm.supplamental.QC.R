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

umap_cluster_colors <- c("#DEA940","#A2A763","#D6EEBE","#5894D4","#6AAD51","#6C8686","#FFD92F", "#7f2201","#b53310",
                         "#E7673B","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6","#386CB0","#FFFFB3", "#5dead5","#c11577",
                         "#CCCCCC", "#7570B3","#B15928", "#6A3D9A", "#d87c6a", "#0060e8", "#84c5ff", "#ffd460", "#d81397")

zm_annotation_final <- read.delim("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_annot_final/Zm.leaf_annot.V5.meta.frozen.txt") %>% 
  mutate(across(c(LouvainClusters),factor))

colnames(zm_annotation_final)


lc_tss <- ggplot(zm_annotation_final, aes(x = LouvainClusters, y = log(tss), fill = LouvainClusters)) +
  geom_violin() + facet_grid(.~sampleID) + 
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values = umap_cluster_colors) + 
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

lc_total <- ggplot(zm_annotation_final, aes(x = LouvainClusters, y = log(total), fill = LouvainClusters)) +
  geom_violin() + facet_grid(.~sampleID) + 
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values = umap_cluster_colors) + 
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

lc_frip <- ggplot(zm_annotation_final, aes(x = LouvainClusters, y = (FRiP), fill = LouvainClusters)) +
  geom_violin() + facet_grid(.~sampleID) + 
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values = umap_cluster_colors) + 
  theme_cowplot() +
  theme(legend.position = "none") + 
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero



lc_ptmt_ratio <- ggplot(zm_annotation_final, aes(x = LouvainClusters, y = (ptmt_ratio), fill = LouvainClusters)) +
  geom_violin() + facet_grid(.~sampleID) + 
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, .25)) # Extend y-axis to include zero


lc_lognsites <- ggplot(zm_annotation_final, aes(x = LouvainClusters, y = log10nSites, fill = LouvainClusters)) +
  geom_violin() + facet_grid(.~sampleID) + 
  geom_boxplot(width=0.1) +
  theme_cowplot( ) + theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = umap_cluster_colors)


lc_doublet_score <- ggplot(zm_annotation_final, aes(x = LouvainClusters, y = doubletscore, fill = LouvainClusters)) +
  geom_violin() + facet_grid(.~sampleID) + 
  geom_boxplot(width=0.1) + scale_fill_manual(values = umap_cluster_colors) + 
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 1)) # Extend y-axis to include zero


colnames(zm_annotation_final_count)
# Stacked + percent
zm_annotation_final_count <- zm_annotation_final %>% 
    mutate(count = 1)

lc_cell_prop <- ggplot(zm_annotation_final_count, aes(fill=LouvainClusters, y=count, x=sampleID)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = umap_cluster_colors)  +
  theme_cowplot() +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero


lc_metrics <- (lc_tss + lc_total)/(lc_frip + lc_ptmt_ratio)/(lc_lognsites + lc_doublet_score)
lc_metrics

ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Zm/Zm.louvain_metrics.pdf",
       lc_metrics,
       width = 18, 
       height = 15,
       units = "in")
       

ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Zm/Zm.loucain_cell_prop.pdf",
       lc_cell_prop,
       width = 5, 
       height = 8,
       units = "in")


########### Same Graphs but for Annotations
#######
library(ggplot2)
library(dplyr)

# Assuming zm_annotation_final and umap_cluster_colors are already defined

  zm_annotation_final <- zm_annotation_final %>%
  mutate(reduce_resolution_annotation = case_when(
    reduce_resolution_annotation == "companion_cells_sieve_elements" ~ "CC/SE",
    TRUE ~ reduce_resolution_annotation # Keeps all other values as they are
  ))

  

# Plot for 'tss'
anno_tss <- ggplot(zm_annotation_final, aes(x = reduce_resolution_annotation, y = log(tss), fill = reduce_resolution_annotation)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'total'
anno_total <- ggplot(zm_annotation_final, aes(x = reduce_resolution_annotation, y = log(total), fill = reduce_resolution_annotation)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'FRiP'
anno_frip <- ggplot(zm_annotation_final, aes(x = reduce_resolution_annotation, y = FRiP, fill = reduce_resolution_annotation)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA)) # Extend y-axis to include zero

# Plot for 'ptmt_ratio'
anno_ptmt_ratio <- ggplot(zm_annotation_final, aes(x = reduce_resolution_annotation, y = ptmt_ratio, fill = reduce_resolution_annotation)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, .25)) # Extend y-axis to include zero

# Plot for 'log10nSites'
anno_lognsites <- ggplot(zm_annotation_final, aes(x = reduce_resolution_annotation, y = log10nSites, fill = reduce_resolution_annotation)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA))

# Plot for 'doubletscore'
anno_doublet_score <- ggplot(zm_annotation_final, aes(x = reduce_resolution_annotation, y = doubletscore, fill = reduce_resolution_annotation)) +
  geom_violin() +
  facet_grid(. ~ sampleID) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = umap_cluster_colors) +
  theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 1)) # Extend y-axis to include zero


colnames(zm_annotation_final_count)
# Stacked + percent
zm_annotation_final_count <- zm_annotation_final %>% 
  mutate(count = 1)

anno_cell_prop <- ggplot(zm_annotation_final_count, aes(fill=reduce_resolution_annotation, y=count, x=sampleID)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = umap_cluster_colors)  +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # Extend y-axis to include zero


unique(zm_annotation_final_count$reduce_resolution_annotation)
anno_metrics <- (anno_tss + anno_total)/(anno_frip + anno_ptmt_ratio)/(  anno_lognsites + anno_doublet_score)
anno_metrics


ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Zm/Zm.anno_metrics.pdf",
       anno_metrics,
       width = 10, 
       height = 15,
       units = "in")

ggsave(filename = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/figures/supplamental/01.QC_Libs/Zm/Zm.anno_cell_prop.pdf",
       anno_cell_prop,
       width = 5, 
       height = 8,
       units = "in")


