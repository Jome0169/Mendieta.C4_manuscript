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



collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4","sb_leaf_svd_strict_compressed_markers.meta.txt"))
colnames(collapsed_meta_data_file)
View(collapsed_meta_data_file)

collapsed_meta_data <- as_tibble(collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))


collapsed_meta_data %>% 
    summarise(count = n())


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
  scale_colour_manual(values=other)  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cell_annotation <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  col_scale  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cell_annotation_knn <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  col_scale  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  col_scale  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

(collapsed_markers_original_clusters + collapse_cell_annotation  + collapse_cell_annotation_knn)
collapsed_markers_original_clusters

colnames(collapsed_meta_data)

sb_final_annotation_assignment <- collapsed_meta_data %>% 
  mutate(V1_annotation = "NA") %>% 
  mutate(V1_annotation = if_else(cell_annotation_knn == cell_annotation_smooth,  cell_annotation_smooth, V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("11", "19"),  "bundle_sheath", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("5", "17", "8", "23"),  "epidermis", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("18", "25", "16", "14"),  "vasculature", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("15", "10", "3", "22", "7"),  "parenchyma", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("9"),  "protoderm", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("1"),  "xylem", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("7", "20"),  "companion_cells", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("21", "4"),  "mesophyll", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("13", "27"),  "vascular_sclerenchyma", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("22"),  "subsidiary_mother_cell", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("26"),  "subsidiary_mother_cell", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(V1_annotation == "NA",  LouvainClusters_t, V1_annotation))
  
  
  

# Generate Cluster Numbers To Append ON -----------------------------------

cluster_numbers <- sb_final_annotation_assignment %>% 
  group_by(V1_annotation) %>% 
  summarise(total_cell_count = n())

sum(cluster_numbers$total_cell_count)


sb_final_annotation_assignment.v1 <- left_join(sb_final_annotation_assignment, cluster_numbers, by = c("V1_annotation")) %>% 
  mutate(V1_annotation_n = str_c(V1_annotation, total_cell_count, sep = "_ncell_"))


updated_annotation <- ggplot(sb_final_annotation_assignment.v1, aes(x=umap1, y = umap2, color = V1_annotation)) +
  scale_colour_manual(values=other)  + guides(colour = guide_legend(override.aes = list(size=5))) + 
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Updated")

collapsed_markers_original_clusters + updated_annotation

updated_annotation



(collapsed_markers_original_clusters + updated_annotation) / (collapse_cell_annotation_knn + collapse_cell_annotation)



write_delim(sb_final_annotation_assignment.v1, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Sb_leaf_svd.knn_100_strict.meta.annotation_V1.txt", 
            col_names = TRUE, quote = "none", delim = "\t")






# Load In V3 Annotation ---------------------------------------------------

generate_color_scheme <- function(meta_file) {
  
  
  
  one_set <- unlist(unique(as.vector(meta_file$cell_annotation_knn)))
  second_set <- unlist(unique(as.vector(meta_file$cell_annotation_smooth)))
  third <- unlist(unique(as.vector(meta_file$cell_annotation_enrich)))
  fourth <- unlist(unique(as.vector(meta_file$V2_annotation)))
  final <- as.factor(unique(c(one_set,second_set,third)))
  
  n <- length(final)
  qual_col_pals =sample( brewer.pal.info[brewer.pal.info$category == 'qual',])
  c25 = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  
  names(c25)  <- levels(final)
  col_scale <- scale_colour_manual(name = "grp", values = c25)
  
  return(col_scale)
  
}

v3_collapsed_meta_data<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/sb_v3","Sb_leaf_svd.V3.subclustering.compressed_markers.meta.txt")) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))

# Use Subclustering Results from V3 ---------------------------------------

### Read in subclustering results
sb_epidermal_subcluster_annotations <- (read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.epidermal.final_annotation.txt", header = TRUE))
sb_vasculature_subcluster_annotations <- (read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.mesophyll.final_annotation.txt", header = TRUE))
sb_mesophyll_subcluster_annotations <- (read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb.vasculature.final_annotation.txt", header = TRUE))%>% 
    filter(cellID %ni% sb_vasculature_subcluster_annotations$cellID )

### Bind both subclustering annnotations and then left joing them
sb.subcluster_annot <- bind_rows(sb_epidermal_subcluster_annotations, sb_vasculature_subcluster_annotations, sb_mesophyll_subcluster_annotations)


sb.subcluster_annot.counts <- sb.subcluster_annot %>% 
  group_by(cellID) %>% 
  summarise(n_times = n()) %>% 
  filter(n_times > 1)

`%ni%` <- Negate(`%in%`)
different_annotations <- sb.subcluster_annot %>% 
    filter(cellID %in% sb.subcluster_annot.counts$cellID) %>% 
    ungroup() %>% 
    distinct() %>% 
    group_by(cellID) %>% 
    filter(n()>1) %>% 
    arrange(cellID) %>% 
    ungroup() %>% 
    filter(V3_final_annnotation != "unknown") %>% 
    group_by(cellID) %>% 
    filter(n()>1) %>% 
    View()
  

same_annotaation_multiple_times <- sb.subcluster_annot %>% 
      filter(cellID %in% sb.subcluster_annot.counts$cellID) %>% 
      ungroup() %>% 
      distinct() %>% 
      group_by(cellID) %>% 
      filter(n() == 1) %>% 
      ungroup() 

remaining_annotations <- sb.subcluster_annot %>% 
  filter(cellID %ni% same_annotaation_multiple_times$cellID) %>% 
  filter(cellID %ni% different_annotations$cellID) 
  



sb.subcluster_annot <- bind_rows(remaining_annotations, same_annotaation_multiple_times)
sb.final.v3_annotation <- left_join(v3_collapsed_meta_data, sb.subcluster_annot,
                                 by = c("cellID")) 

sb_final_colors <- generate_color_scheme(sb.final.v3_annotation)



collapsed_markers_original_clusters <-  ggplot(sb.final.v3_annotation, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other) + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


final_v3_cluster <- ggplot(sb.final.v3_annotation, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8, show.legend = TRUE) + theme_minimal() + ggtitle("Cell Annotation V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

collapse_cell_annotation_knn_v3 <- ggplot(sb.final.v3_annotation, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  sb_final_colors  + geom_point(size = .25, alpha = .8, show.legend = TRUE) + theme_minimal() + ggtitle("Cell Annotation Knn V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cell_annotation_smooth_v3 <- ggplot(sb.final.v3_annotation, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  sb_final_colors  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation smooth V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapsed_markers_original_clusters
collapsed_markers_original_clusters  + collapse_cell_annotation_knn_v3 + collapse_cell_annotation_smooth_v3 + final_v3_cluster
collapsed_markers_original_clusters + collapse_cell_annotation_knn_v3
collapsed_markers_original_clusters + final_v3_cluster


sb.final.v3_annotation %>% 
  filter(is.na(V3_final_annnotation) == TRUE) %>% 
  ggplot(., aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8, show.legend = TRUE) + theme_minimal() + ggtitle("Cell Annotation V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


sb.final.v3_annotation.final <- sb.final.v3_annotation %>% 
  #mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("28", "33"), "xylem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("31"), "bundle_sheath", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("28"), "xylem", V3_final_annnotation)) %>% 
  #mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("27"), "phloem", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(is.na(V3_final_annnotation) == TRUE, "unknown", V3_final_annnotation))

  
  
collapsed_markers_original_clusters

ggplot(sb.final.v3_annotation.final, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8, show.legend = TRUE) + theme_minimal() + ggtitle("Cell Annotation V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

sb.final.v3_annotation.final %>% 
  filter(V3_final_annnotation == "unknown") %>% 
  ggplot(., aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8, show.legend = TRUE) + theme_minimal() + ggtitle("Cell Annotation V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

colnames(sb.final.v3_annotation.final)


final_cell_count <- sb.final.v3_annotation.final %>% 
  group_by(V3_final_annnotation) %>% 
  summarise(total_cell_count = n())

View(final_cell_count)

sb_final_annotation_assignment.v3.final <- left_join(sb.final.v3_annotation.final, final_cell_count, by = c("V3_final_annnotation")) %>% 
  mutate(V3_annotation_n = str_c(V3_final_annnotation, total_cell_count, sep = "_ncell_")) 

unique(sb_final_annotation_assignment.v3.final$V3_annotation_n)

write_delim(sb_final_annotation_assignment.v3.final, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Sb_leaf.V3_final.txt", 
            col_names = TRUE, quote = "none", delim = "\t")



# Generate Final UMAPs ----------------------------------------------------

### note there is an annoying thing here wher the colors set in col_scale are generated by the maize 
### annotation script. So you NEED to run the maize script before you can get the correct colors here.
unique(sb_final_annotation_assignment.v3.final$V3_final_annnotation)

sb.fig_1_umap_pdf <- ggplot(sb_final_annotation_assignment.v3.final, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  theme_half_open() +
  col_scale + geom_point(size = .25, alpha = .8) + ggtitle("Sb Cell Annotation Final")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



sb.fig_1_umap_pdf
ggsave("sb.cell_annot.pdf", 
       plot = sb.fig_1_umap_pdf, 
       path = "/Users/feilab/Projects/05.comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 10, height = 10)



sb.fig_1_umap_png <- ggplot(sb_final_annotation_assignment.v3.final, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  theme_void() + 
  col_scale + geom_point(size = .5, alpha = .8, show.legend = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5)))


ggsave("sb.cell_annot.png", 
       plot = sb.fig_1_umap_png, 
       path = "/Users/feilab/Projects/05.comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 5, height = 5)


collapsed_markers_original_clusters + sb.fig_1_umap_pdf
sb_final_annotation_assignment.v3.final

sb_final_annotation_assignment.v3.final %>% 
  group_by(LouvainClusters_t, V3_final_annnotation) %>% 
  summarise(counts = n()) %>% 
  arrange(LouvainClusters_t, counts, desc = TRUE) %>% 
  View()


# Evo Chromo Figure -------------------------------------------------------


# Evo Chromo Figure
###,Colors,for,UMAP

library(scales)
library(randomcoloR)

?randomColor
show_col(randomColor(count = 100, hue = "pink"))


collapsed_meta_data_file<- read_delim(here("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4", "Sb.leaf_annot.V4.meta.final.2022-10-10.txt"),
                                      col_names = TRUE)


umap_cluster_colors <- c("#DEA940","#A2A763","#D6EEBE","#5894D4","#6AAD51","#6C8686","#FFD92F", "#7f2201","#b53310",
                         "#E7673B","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6","#386CB0","#FFFFB3","#589e1c", "#5dead5","#c11577",
                         "#CCCCCC", "#7570B3","#B15928", "#6A3D9A", "#d87c6a", "#0060e8", "#84c5ff", "#ffd460", "#d81307", "#589e1e")

my_pals <- show_col(umap_cluster_colors)

colnames(collapsed_meta_data_file)
library(extrafont)
sb_collapsed_markers_original_clusters <- 
  collapsed_meta_data_file %>%  
  mutate(across(c(LouvainClusters),
                factor)) %>% 
    ggplot(., aes(x=umap1, y = umap2, color = LouvainCluster_t)) +
  scale_colour_manual(values=umap_cluster_colors)  + geom_point(size = .5, alpha = 1) + ggtitle("Sorghum bicolor Intial Clustering") +
  guides(colour = guide_legend(override.aes = list(size=5))) + theme_cowplot()

sb_collapsed_markers_original_clusters




ggsave("sb.initial_clustering.2.pdf", 
       plot = sb_collapsed_markers_original_clusters, 
       path = "/Users/feilab/Documents/11.Presentations/EvoChromo_poster/figures",
       units = "in",
       width = 10, height = 10)


colnames(collapsed_meta_data_file)
sb.annotations <-  collapsed_meta_data_file %>%  
  ggplot(., aes(x=umap1, y = umap2, color = sb_v4_annot)) +
  scale_colour_manual(values=umap_cluster_colors)  + geom_point(size = .5, alpha = 1) + ggtitle("Sorghum bicolor Annotation") +
  guides(colour = guide_legend(override.aes = list(size=5))) + theme_cowplot()

sb.annotations
