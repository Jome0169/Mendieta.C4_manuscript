library("here")
library(devtools)
library(Seurat)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(RColorBrewer)

load_all('/Users/feilab/Programming/Socrates')

n <- 30
nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)
library(randomcoloR)


# Vasculature Annotation --------------------------------------------------


zm.vascualer.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/zm_v3","Zm_leaf_svd.V3.subclustering.vasculature_bins.compressed_markers.meta.txt"))
collapsed_meta_data <- as_tibble(zm.vascualer.collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))

colnames(collapsed_meta_data)

one_set <- unlist(unique(as.vector(collapsed_meta_data$cell_annotation_knn)))
second_set <- unlist(unique(as.vector(collapsed_meta_data$cell_annotation_smooth)))
third <- unlist(unique(as.vector(collapsed_meta_data$cell_annotation_enrich)))
final <- as.factor(unique(c(one_set,second_set,third)))

n <- length(final)
qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
c25 = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

names(c25)  <- levels(final)
col_scale <- scale_colour_manual(name = "grp", values = c25)

collapsed_markers_original_clusters <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cell_annotation <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

collapse_cell_annotation

collapse_cell_annotation_knn <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cluster_annotation <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

collapsed_markers_original_clusters


(collapsed_markers_original_clusters  + collapse_cell_annotation_knn + collapse_cell_annotation + collapse_cluster_annotation)




zm.vascualer.final_annotation <-  zm.vascualer.collapsed_meta_data_file  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("5"), "xylem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8"), "metaphloem_sieve_element", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1"), "companion_cells", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("2","3"), "vascular_parenchyma", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("7"), "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("6"), "epidermis", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("9"), "unknown", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("4"), "procambium", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("6"), "guard_cell", V3_final_annnotation))
  


ggplot(zm.vascualer.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Vasculature Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


zm.vascualer.final_annotation_writeable <- zm.vascualer.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(zm.vascualer.final_annotation_writeable, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf_svd.V3.subclustering.vasculature.annotations.txt", 
            col_names = TRUE, quote = "none", delim = "\t")





# Epidermal Subclustering Annotation --------------------------------------

zm.epidermal.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/zm_v3","Zm_leaf_svd.V3.subclustering.epidermal_bins.compressed_markers.meta.txt"))
collapsed_meta_data.epidermal <- as_tibble(zm.epidermal.collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))




epidermal.one_set <- unlist(unique(as.vector(collapsed_meta_data.epidermal$cell_annotation_knn)))
epidermal.second_set <- unlist(unique(as.vector(collapsed_meta_data.epidermal$cell_annotation_smooth)))
epidermal.third <- unlist(unique(as.vector(collapsed_meta_data.epidermal$cell_annotation_enrich)))
epidermal.final <- as.factor(unique(c(epidermal.one_set,epidermal.second_set,epidermal.third)))

epidermal.n <- length(epidermal.final)
epidermal.qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
epidermal.c25 = sample(unlist(mapply(brewer.pal, epidermal.qual_col_pals$maxcolors, rownames(epidermal.qual_col_pals))))

names(epidermal.c25)  <- levels(epidermal.final)
epidermal.col_scale <- scale_colour_manual(name = "grp", values = epidermal.c25)

epidermal.collapsed_markers_original_clusters <- ggplot(collapsed_meta_data.epidermal, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

epidermal.collapse_cell_annotation <- ggplot(collapsed_meta_data.epidermal, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  epidermal.col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

epidermal.collapse_cell_annotation_knn <- ggplot(collapsed_meta_data.epidermal, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  epidermal.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


epidermal.collapse_cluster_annotation <- ggplot(collapsed_meta_data.epidermal, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  epidermal.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


epidermal.collapsed_markers_original_clusters
epidermal.collapsed_markers_original_clusters + epidermal.collapse_cell_annotation + epidermal.collapse_cell_annotation_knn + epidermal.collapse_cluster_annotation


epidermal.collapsed_markers_original_clusters

zm.epidermal.final_annotation <-  collapsed_meta_data.epidermal  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1", "4", "5", "9", "6"), "epidermal", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8", "7", "14", "10", "11"), "protoderm", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("3", "2", "13"), "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("15"), "subsidary_cell", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("12", "17", '16'), "unknown", V3_final_annnotation))
  

final.epidermal_annotation_graph <- ggplot(zm.epidermal.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Vasculature Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


epidermal.collapsed_markers_original_clusters + final.epidermal_annotation_graph

zm.epidermal.final_annotation.writable <- zm.epidermal.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(zm.epidermal.final_annotation.writable, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf_svd.V3.subclustering.epidermal.annotations.txt", 
            col_names = TRUE, quote = "none", delim = "\t")

