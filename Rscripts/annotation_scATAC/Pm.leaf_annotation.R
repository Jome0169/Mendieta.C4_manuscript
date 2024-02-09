library(tidyverse)
library("here")
library(devtools)
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


pm_subcluster_meta_file <- here("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/pm_annot/Pm_leaf.merged_replicates.Subclustering_vasculature.SVD.full.metadata.txt")
pm_subcluster_leaf <- read.table(pm_subcluster_meta_file)


pm_subcluster_leaf_annotated <- pm_subcluster_leaf %>% 
    mutate(subcluster_lc_safe = str_c("Lc", LouvainClusters, sep = "_")) %>% 
    dplyr::filter(subcluster_lc_safe != "Lc_6") %>% 
    mutate(annotation_v1 = case_when(subcluster_lc_safe == "Lc_1" ~ "bundle_sheath", 
                                        subcluster_lc_safe == "Lc_2" ~ "mesophyll;developing", 
                                        subcluster_lc_safe == "Lc_3" ~ "protoderm", 
                                        subcluster_lc_safe == "Lc_4" ~ "procambium", 
                                        subcluster_lc_safe == "Lc_5" ~ "companion_cells_sieve_elements")) %>% 
    dplyr::rename(umap1.subcluster = umap1, umap2.subcluster = umap2) %>% 
    dplyr::select(cellID, subcluster_lc_safe, annotation_v1,umap1.subcluster,umap2.subcluster)





ggplot(pm_subcluster_leaf_annotated, aes(x=umap1.subcluster, y = umap2.subcluster, color = as.factor(annotation_v1))) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



pm_force_annotation_file <- here("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/pm_annot/Pm_leaf_annotation.meta.txt")
pm_force_annotation <- read.table(pm_force_annotation_file)




pm_force_annotation_removed_subcluster_cells <- pm_force_annotation %>% 
    dplyr::filter(cellID %ni% pm_subcluster_leaf$cellID) %>% 
    dplyr::filter(Louvain_cluster_safe != "LC_2") %>% 
    mutate(annotation_v1 = case_when(Louvain_cluster_safe == "LC_1" ~ "mesophyll", 
                                        Louvain_cluster_safe == "LC_4" ~ "epidermis", 
                                        Louvain_cluster_safe == "LC_5" ~ "mesophyll;developing", 
                                      Louvain_cluster_safe == "LC_6" ~ "protoderm"))
  
Pm_Lcs <- ggplot(pm_force_annotation_removed_subcluster_cells, aes(x=umap1, y = umap2, color = Louvain_cluster_safe)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))
Pm_Lcs


pm_force_annotation_subcluster_cells <- pm_force_annotation %>% 
  dplyr::filter(cellID %in% pm_subcluster_leaf_annotated$cellID) %>% 
  left_join(.,pm_subcluster_leaf_annotated, by = "cellID")


pm_combined_annotation <- bind_rows(pm_force_annotation_removed_subcluster_cells, pm_force_annotation_subcluster_cells)

pm_annotation_ggplot <- ggplot(pm_combined_annotation, aes(x=umap1, y = umap2, color = annotation_v1)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

pm_annotation_ggplot


### Find Leftover Cells based on PCA Proximity
## Select down the cells to pull 
pm.PCAs <- read.table(here("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/pm_annot/Pm_leaf.merged_replicates.SVD.full.reduced_dimensions.txt"))

pm.na_annots <- pm_combined_annotation %>% 
  filter(is.na(annotation_v1) == TRUE)
pm.annot <- pm_combined_annotation %>% 
  filter(is.na(annotation_v1) != TRUE)

annotate_PCAs <- pm.PCAs[pm.annot$cellID,]
non_annotated_PCAs <- pm.PCAs[pm.na_annots$cellID,]


knn_pred = class::knn(annotate_PCAs, non_annotated_PCAs, pm.annot$annotation_v1, k = 5, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
pm.na_annots$annotation_v1 <- knn_pred

View(pm.na_annots)
pm.total_annotation <- bind_rows(pm.na_annots, pm.annot) %>% 
  mutate(annotation_v1_safe = str_replace_all(annotation_v1, ";", "_"))


pm_annotation_v1 <- ggplot(pm.total_annotation, aes(x=umap1, y = umap2, color = annotation_v1)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

pm_annotation_v1

pm.total_annotation_class_counts <- pm.total_annotation %>% 
    group_by(annotation_v1) %>% 
    summarise(counts = n())

combined_annotation <- left_join(pm.total_annotation, pm.total_annotation_class_counts, by = ("annotation_v1")) %>% 
    mutate(annotation_v1_ncell = str_c(annotation_v1, counts, sep = "_ncell_")) %>% 
    dplyr::mutate(reduce_resolution_annotation = case_when(annotation_v1 == "mesophyll;developing" ~ "mesophyll",
                                                         TRUE ~ annotation_v1))

pm_annotation_v1 <- ggplot(combined_annotation, aes(x=umap1, y = umap2, color = reduce_resolution_annotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))
pm_annotation_v1


write_delim(combined_annotation, "/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/pm_annot/Pm.leaf_annotation.V1.meta.txt", 
            col_names = TRUE, quote = "none", delim = "\t")



