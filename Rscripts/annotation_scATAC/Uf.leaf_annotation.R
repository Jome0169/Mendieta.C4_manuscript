library(tidyverse)
library("here")
library(devtools)
library(Seurat)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(RColorBrewer)
library(patchwork)

`%ni%` <- Negate(`%in%`)

n <- 30
nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)



uf_subcluster_meta_file <- here("Mendieta_et_al_comparative_single_cell/metrics/annotations/uf_annot/Uf_leaf.merged_replicates.Subclustering.SVD.full.metadata.txt")
uf_subcluster_leaf <- read.table(uf_subcluster_meta_file)

uf_subcluster_leaf$LouvainClusters

uf_subcluster_leaf_annotated <- uf_subcluster_leaf %>% 
  mutate(subcluster_lc_safe = str_c("Lc", LouvainClusters, sep = "_")) %>% 
  dplyr::filter(subcluster_lc_safe != "Lc_6") %>% 
  mutate(annotation_v1 = case_when(subcluster_lc_safe == "Lc_1" ~ "epidermis", 
                                   subcluster_lc_safe == "Lc_2" ~ "mesophyll", 
                                   subcluster_lc_safe == "Lc_3" ~ "bundle_sheath", 
                                   subcluster_lc_safe == "Lc_4" ~ "mesophyll", 
                                   subcluster_lc_safe == "Lc_5" ~ "bundle_sheath")) %>% 
  dplyr::rename(umap1.subcluster = umap1, umap2.subcluster = umap2) %>% 
  dplyr::select(cellID, subcluster_lc_safe, annotation_v1,umap1.subcluster,umap2.subcluster)





ggplot(uf_subcluster_leaf_annotated, aes(x=umap1.subcluster, y = umap2.subcluster, color = as.factor(annotation_v1))) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation SubClsuter")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



uf_force_annotation_file <- here("Mendieta_et_al_comparative_single_cell/metrics/annotations/uf_annot/Uf_leaf.merged_replicates.SVD.full.metadata.txt")
uf_force_annotation <- read.table(uf_force_annotation_file)

uf_force_annotation$Louvain_cluster_safe

uf_force_annotation_removed_subcluster_cells <- uf_force_annotation %>% 
  mutate(Louvain_cluster_safe = str_c("LC", LouvainClusters, sep = "_")) %>% 
  dplyr::filter(cellID %ni% uf_subcluster_leaf_annotated$cellID) %>% 
  dplyr::filter(Louvain_cluster_safe != "LC_15") %>% 
  dplyr::filter(Louvain_cluster_safe != "LC_16") %>% 
  mutate(annotation_v1 = case_when(Louvain_cluster_safe == "LC_1" ~ "epidermis", 
                                   Louvain_cluster_safe == "LC_2" ~ "mesophyll", 
                                   Louvain_cluster_safe == "LC_3" ~ "bundle_sheath", 
                                   Louvain_cluster_safe == "LC_4" ~ "bundle_sheath", 
                                   Louvain_cluster_safe == "LC_5" ~ "bundle_sheath", 
                                   Louvain_cluster_safe == "LC_6" ~ "epidermis", 
                                   Louvain_cluster_safe == "LC_7" ~ "mesophyll", 
                                   Louvain_cluster_safe == "LC_8" ~ "companion_cells_sieve_elements",
                                   Louvain_cluster_safe == "LC_9" ~ "bundle_sheath"))


uf_Lcs <- ggplot(uf_force_annotation_removed_subcluster_cells, aes(x=umap1, y = umap2, color = Louvain_cluster_safe)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



uf_force_annotation_subcluster_cells <- uf_force_annotation %>% 
  dplyr::filter(cellID %in% uf_subcluster_leaf_annotated$cellID) %>% 
  left_join(.,uf_subcluster_leaf_annotated, by = "cellID")


uf_combined_annotation <- bind_rows(uf_force_annotation_removed_subcluster_cells, uf_force_annotation_subcluster_cells)

uf_annotation_ggplot <- ggplot(uf_combined_annotation, aes(x=umap1, y = umap2, color = annotation_v1)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

library(patchwork)
uf_Lcs + uf_annotation_ggplot


### Find Leftover Cells based on PCA Proximity

## Select down the cells to pull 
uf.PCAs <- read.table(here("Mendieta_et_al_comparative_single_cell/metrics/annotations/uf_annot/Uf_leaf.merged_replicates.SVD.full.reduced_dimensions.txt"))

uf.na_annots <- uf_combined_annotation %>% 
  filter(is.na(annotation_v1) == TRUE)
uf.annot <- uf_combined_annotation %>% 
  filter(is.na(annotation_v1) != TRUE)

annotate_PCAs <- uf.PCAs[uf.annot$cellID,]
non_annotated_PCAs <- uf.PCAs[uf.na_annots$cellID,]


knn_pred = class::knn(annotate_PCAs, non_annotated_PCAs, uf.annot$annotation_v1, k = 5, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
uf.na_annots$annotation_v1 <- knn_pred

View(uf.na_annots)
uf.total_annotation <- bind_rows(uf.na_annots, uf.annot)

uf_annotation_v1 <- ggplot(uf.total_annotation, aes(x=umap1, y = umap2, color = annotation_v1)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

uf.total_annotation_class_counts <- uf.total_annotation %>% 
  group_by(annotation_v1) %>% 
  summarise(counts = n())

combined_annotation <- left_join(uf.total_annotation, uf.total_annotation_class_counts, by = ("annotation_v1")) %>% 
  mutate(annotation_v1_ncell = str_c(annotation_v1, counts, sep = "_ncell_"))

unique(combined_annotation$annotation_v1)

unique(combined_annotation$annotation_v1_ncell)
write_delim(combined_annotation, "/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/uf_annot/uf.leaf_annotation.V1.meta.txt", 
            col_names = TRUE, quote = "none", delim = "\t")



