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
# Epidermal Subclustering Annotation --------------------------------------

sb.epidermal.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.epidermal.meta.txt"))
collapsed_meta_data.epidermal <- as_tibble(sb.epidermal.collapsed_meta_data_file) %>% 
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



epidermal.collapsed_markers_original_clusters + epidermal.collapse_cell_annotation + epidermal.collapse_cell_annotation_knn + epidermal.collapse_cluster_annotation


epidermal.collapsed_markers_original_clusters

sb.epidermal.final_annotation <-  collapsed_meta_data.epidermal  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1", "2", "3", "4", "5", "11", "12", "13"), "epidermal", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("6"), "protoderm", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("7"), "ground_meristem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8", "10", "9"), "protoderm", V3_final_annnotation)) %>%
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1"), "unknown", V3_final_annnotation))


final.epidermal_annotation_graph <- ggplot(sb.epidermal.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Vasculature Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


epidermal.collapsed_markers_original_clusters + final.epidermal_annotation_graph

sb.epidermal.final_annotation.writable <- sb.epidermal.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(sb.epidermal.final_annotation.writable, 
            "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.epidermal.final_annotation.txt", 
            col_names = TRUE, quote = "none", delim = "\t")




# Vascualture Subclsutering  ----------------------------------------------
n <- 30
nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)
library(randomcoloR)


sb.vasculature.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.vasculature.meta.txt"))
collapsed_meta_data.vasculature <- as_tibble(sb.vasculature.collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))


vasculature.one_set <- unlist(unique(as.vector(collapsed_meta_data.vasculature$cell_annotation_knn)))
vasculature.second_set <- unlist(unique(as.vector(collapsed_meta_data.vasculature$cell_annotation_smooth)))
vasculature.third <- unlist(unique(as.vector(collapsed_meta_data.vasculature$cell_annotation_enrich)))
vasculature.final <- as.factor(unique(c(vasculature.one_set,vasculature.second_set,vasculature.third)))

vasculature.n <- length(vasculature.final)
vasculature.qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
vasculature.c25 = sample(unlist(mapply(brewer.pal, vasculature.qual_col_pals$maxcolors, rownames(vasculature.qual_col_pals))))

names(vasculature.c25)  <- levels(vasculature.final)
vasculature.col_scale <- scale_colour_manual(name = "grp", values = vasculature.c25)

vasculature.collapsed_markers_original_clusters <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

vasculature.collapse_cell_annotation <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  vasculature.col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

vasculature.collapse_cell_annotation_knn <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  vasculature.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


vasculature.collapse_cluster_annotation <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  vasculature.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cluster Annotation smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



vasculature.collapsed_markers_original_clusters + vasculature.collapse_cell_annotation + vasculature.collapse_cell_annotation_knn + vasculature.collapse_cluster_annotation


vasculature.collapsed_markers_original_clusters

sb.vasculature.final_annotation <-  collapsed_meta_data.vasculature  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1","13"), "unknown", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("2","3","14"), "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("6"), "vascular_parenchyma", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("5", "12"), "bundle_sheath", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8"), "sieve_element", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("11"), "procambial_meristem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("7", "9"), "companion_cell", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("4", "10"), "xylem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("13"), "unknown", V3_final_annnotation))


final.vasculature_annotation_graph <- ggplot(sb.vasculature.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Vasculature Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


vasculature.collapsed_markers_original_clusters + final.vasculature_annotation_graph

sb.vasculature.final_annotation.writable <- sb.vasculature.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(sb.vasculature.final_annotation.writable, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.vasculature.final_annotation.txt", 
            col_names = TRUE, quote = "none", delim = "\t")




# Mesophyll Subclsutering  ----------------------------------------------
n <- 30
nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)
library(randomcoloR)


sb.mesophyll.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.mesophyll.meta.txt"))
collapsed_meta_data.mesophyll <- as_tibble(sb.mesophyll.collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))


mesophyll.one_set <- unlist(unique(as.vector(collapsed_meta_data.mesophyll$cell_annotation_knn)))
mesophyll.second_set <- unlist(unique(as.vector(collapsed_meta_data.mesophyll$cell_annotation_smooth)))
mesophyll.third <- unlist(unique(as.vector(collapsed_meta_data.mesophyll$cell_annotation_enrich)))
mesophyll.final <- as.factor(unique(c(mesophyll.one_set,mesophyll.second_set,mesophyll.third)))

mesophyll.n <- length(mesophyll.final)
mesophyll.qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
mesophyll.c25 = sample(unlist(mapply(brewer.pal, mesophyll.qual_col_pals$maxcolors, rownames(mesophyll.qual_col_pals))))

names(mesophyll.c25)  <- levels(mesophyll.final)
mesophyll.col_scale <- scale_colour_manual(name = "grp", values = mesophyll.c25)

mesophyll.collapsed_markers_original_clusters <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

mesophyll.collapse_cell_annotation <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  mesophyll.col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

mesophyll.collapse_cell_annotation_knn <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  mesophyll.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


mesophyll.collapse_cluster_annotation <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  mesophyll.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



mesophyll.collapsed_markers_original_clusters + mesophyll.collapse_cell_annotation + mesophyll.collapse_cell_annotation_knn + mesophyll.collapse_cluster_annotation


sb.mesophyll.final_annotation <-  collapsed_meta_data.mesophyll  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("2", "3", "8", "11", "12", "14"), "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("6"), "developing_mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1", "13"), "unknown", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("5"), "bundle_sheath", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("4"), "epidermal", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("7", "9"), "ground_meristem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("10"), "vascular_parenchyma", V3_final_annnotation))


final.mesophyll_annotation_graph <- ggplot(sb.mesophyll.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation mesophyll Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


mesophyll.collapsed_markers_original_clusters + final.mesophyll_annotation_graph

sb.mesophyll.final_annotation.writable <- sb.mesophyll.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(sb.mesophyll.final_annotation.writable, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.mesophyll.final_annotation.txt", 
            col_names = TRUE, quote = "none", delim = "\t")



# Re-Annotating Subclustering Results: 2022-08-11 -------------------------


### Epidermal Annotation 
sb.epidermal.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.epidermal.V2.meta.txt"))
collapsed_meta_data.epidermal <- as_tibble(sb.epidermal.collapsed_meta_data_file) %>% 
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



epidermal.collapsed_markers_original_clusters + epidermal.collapse_cell_annotation + epidermal.collapse_cell_annotation_knn + epidermal.collapse_cluster_annotation


epidermal.collapsed_markers_original_clusters

sb.epidermal.final_annotation <-  collapsed_meta_data.epidermal  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("2", "3", "5", "14", "15", "10"), "epidermal", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("6", "7", "11", "9", "4"), "protoderm", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8", "12"), "ground_meristem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("13"), "subsidary_cell", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1"), "unknown", V3_final_annnotation))


final.epidermal_annotation_graph <- ggplot(sb.epidermal.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Vasculature Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


epidermal.collapsed_markers_original_clusters + final.epidermal_annotation_graph

sb.epidermal.final_annotation.writable <- sb.epidermal.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(sb.epidermal.final_annotation.writable, 
            "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.epidermal.final_annotation.txt", 
            col_names = TRUE, quote = "none", delim = "\t")



### Vasculature Annotations
n <- 30
nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)
library(randomcoloR)


sb.vasculature.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.vasculature.V2.meta.txt"))
collapsed_meta_data.vasculature <- as_tibble(sb.vasculature.collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))


vasculature.one_set <- unlist(unique(as.vector(collapsed_meta_data.vasculature$cell_annotation_knn)))
vasculature.second_set <- unlist(unique(as.vector(collapsed_meta_data.vasculature$cell_annotation_smooth)))
vasculature.third <- unlist(unique(as.vector(collapsed_meta_data.vasculature$cell_annotation_enrich)))
vasculature.final <- as.factor(unique(c(vasculature.one_set,vasculature.second_set,vasculature.third)))

vasculature.n <- length(vasculature.final)
vasculature.qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
vasculature.c25 = sample(unlist(mapply(brewer.pal, vasculature.qual_col_pals$maxcolors, rownames(vasculature.qual_col_pals))))

names(vasculature.c25)  <- levels(vasculature.final)
vasculature.col_scale <- scale_colour_manual(name = "grp", values = vasculature.c25)

vasculature.collapsed_markers_original_clusters <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

vasculature.collapse_cell_annotation <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  vasculature.col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

vasculature.collapse_cell_annotation_knn <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  vasculature.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


vasculature.collapse_cluster_annotation <- ggplot(collapsed_meta_data.vasculature, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  vasculature.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cluster Annotation smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



vasculature.collapsed_markers_original_clusters + vasculature.collapse_cell_annotation + vasculature.collapse_cell_annotation_knn + vasculature.collapse_cluster_annotation


vasculature.collapsed_markers_original_clusters

sb.vasculature.final_annotation <-  collapsed_meta_data.vasculature  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1","9"), "unknown", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("2","3", "6","17", "18"), "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("12", "13"), "vascular_parenchyma", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("4", "7", "11"), "bundle_sheath", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("5", "10","15", "16", "19"), "sieve_element", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8"), "companion_cell", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("14"), "epidermal", V3_final_annnotation))


final.vasculature_annotation_graph <- ggplot(sb.vasculature.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Vasculature Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


vasculature.collapsed_markers_original_clusters + final.vasculature_annotation_graph

sb.vasculature.final_annotation.writable <- sb.vasculature.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(sb.vasculature.final_annotation.writable, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.vasculature.final_annotation.txt", 
            col_names = TRUE, quote = "none", delim = "\t")

#### Mesophyll Subclsutering
n <- 30
nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)
library(randomcoloR)


sb.mesophyll.collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.mesophyll.V2.meta.txt"))
collapsed_meta_data.mesophyll <- as_tibble(sb.mesophyll.collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))


mesophyll.one_set <- unlist(unique(as.vector(collapsed_meta_data.mesophyll$cell_annotation_knn)))
mesophyll.second_set <- unlist(unique(as.vector(collapsed_meta_data.mesophyll$cell_annotation_smooth)))
mesophyll.third <- unlist(unique(as.vector(collapsed_meta_data.mesophyll$cell_annotation_enrich)))
mesophyll.final <- as.factor(unique(c(mesophyll.one_set,mesophyll.second_set,mesophyll.third)))

mesophyll.n <- length(mesophyll.final)
mesophyll.qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
mesophyll.c25 = sample(unlist(mapply(brewer.pal, mesophyll.qual_col_pals$maxcolors, rownames(mesophyll.qual_col_pals))))

names(mesophyll.c25)  <- levels(mesophyll.final)
mesophyll.col_scale <- scale_colour_manual(name = "grp", values = mesophyll.c25)

mesophyll.collapsed_markers_original_clusters <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

mesophyll.collapse_cell_annotation <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  mesophyll.col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

mesophyll.collapse_cell_annotation_knn <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  mesophyll.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


mesophyll.collapse_cluster_annotation <- ggplot(collapsed_meta_data.mesophyll, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  mesophyll.col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



mesophyll.collapsed_markers_original_clusters + mesophyll.collapse_cell_annotation + mesophyll.collapse_cell_annotation_knn + mesophyll.collapse_cluster_annotation


sb.mesophyll.final_annotation <-  collapsed_meta_data.mesophyll  %>% 
  mutate(V3_final_annnotation = "NA") %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("2", "6", "7"), "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("1"), "unknown", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("3", "9", "10"), "bundle_sheath", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("5","4"), "epidermal", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters %in% c("8"), "ground_meristem", V3_final_annnotation))
  


final.mesophyll_annotation_graph <- ggplot(sb.mesophyll.final_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation mesophyll Final")  +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_colour_manual(values=other)


mesophyll.collapsed_markers_original_clusters + final.mesophyll_annotation_graph

sb.mesophyll.final_annotation.writable <- sb.mesophyll.final_annotation %>% 
  dplyr::select(cellID,V3_final_annnotation)

write_delim(sb.mesophyll.final_annotation.writable, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.mesophyll.final_annotation.txt", 
            col_names = TRUE, quote = "none", delim = "\t")
