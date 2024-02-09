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



collapsed_meta_data_file<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/zm_v1","Zm_leaf_svd_split_markers.knn_100_strict.meta.txt"))
colnames(collapsed_meta_data_file)


collapsed_meta_data <- as_tibble(collapsed_meta_data_file) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))



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

collapsed_markers_original_clusters


collapse_cell_annotation <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Smooth")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

collapse_cell_annotation

collapse_cell_annotation_knn <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  col_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cluster_annotation <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cluster_annotation_smooth)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Knn")  +
  guides(colour = guide_legend(override.aes = list(size=5)))




colnames((collapsed_meta_data))
enriched_cell_annotation <- ggplot(collapsed_meta_data, aes(x=umap1, y = umap2, color = cell_annotation_enrich)) +
  col_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Enriched Cell Annotation")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

enriched_cell_annotation


(collapsed_markers_original_clusters  + collapse_cell_annotation_knn + collapse_cell_annotation)


colnames(collapsed_meta_data)
unique(collapsed_meta_data$cell_annotation_knn)

zm_final_annotation_assignment <- collapsed_meta_data %>% 
  mutate(V1_annotation = "NA") %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("10", "2", "3"),  "mesophyll", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("8", "6"),  "bundle_sheath", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("9", "12", "17", "19", "18", "16"),  "vasculature", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("4"),  "sieve_element_precursors_companion_cells", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("12"),  "xylem_parenchya", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("13", "11","15"),  "leaf_primordia", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("16"),  "companion_cell", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("1"),  "epidermis", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("9"),  "phloem", V1_annotation))  %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("14"),  "protoderm", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("7"),  "subsidary_cells", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("18"),  "unknown", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("17"),  "subsidary_cells", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("16"),  "phloem_parenchyma", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("11"),  "mesophyll", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(LouvainClusters_t %in% c("5"),  "epidermis", V1_annotation)) %>% 
  mutate(V1_annotation = if_else(V1_annotation == "NA",  LouvainClusters_t, V1_annotation))
  
  #mutate(V1_annotation = if_else(cell_annotation_knn == cell_annotation_smooth,  cell_annotation_smooth, V1_annotation))
  
updated_annotation <- ggplot(zm_final_annotation_assignment, aes(x=umap1, y = umap2, color = V1_annotation)) +
  scale_colour_manual(values=other)  + guides(colour = guide_legend(override.aes = list(size=5))) + 
 geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Updated")

collapsed_markers_original_clusters / updated_annotation

(collapsed_markers_original_clusters + updated_annotation) / (collapse_cell_annotation_knn + collapse_cell_annotation)

zm_final_annotation_assignment

write_delim(zm_final_annotation_assignment, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf_svd.knn_100_strict.meta.annotation_V1.txt", 
            col_names = TRUE, quote = "none", delim = "\t")




# Updated Annotation 2  ---------------------------------------------------
#Many of the annotation updates and why they were chosen were 
#faciliated by the heatmap as well as comments from the 
#lab and browswer screenshots. Details can additionally be 
#found on the github issues:https://github.com/Jome0169/Mendieta_et_al_comparative_single_cell/issues/20

zm_final_annotation_assignment.v2 <- zm_final_annotation_assignment %>% 
  mutate(V2_annotation = V1_annotation) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("7"),  "protoderm", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("15"),  "developing_bundle_sheath", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("4"),  "phloem_SE", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("16"),  "sclerenchymaous_bundle_sheath", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("20"),  "xylem", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("21"),  "vasculature", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("17"),  "guard_cell", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("22"),  "vasculature", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("19"),  "guard_mother_cell", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("9"),  "companion_cell", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("20"),  "xylem", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("21"),  "vascular_phloem_like", V2_annotation)) %>% 
  mutate(V2_annotation = if_else(LouvainClusters_t %in% c("22"),  "vascular_parenchyma_like", V2_annotation))
  


cluster_numbers <- zm_final_annotation_assignment.v2 %>% 
    group_by(V2_annotation) %>% 
    summarise(total_cell_count = n())


zm_final_annotation_assignment.v2 <- left_join(zm_final_annotation_assignment.v2, cluster_numbers, by = c("V2_annotation")) %>% 
    mutate(V2_annotation_n = str_c(V2_annotation, total_cell_count, sep = "_ncell_"))

unique(zm_final_annotation_assignment.v2$V2_annotation_n)

updated_annotation_2 <- ggplot(zm_final_annotation_assignment.v2, aes(x=umap1, y = umap2, color = V2_annotation)) +
  scale_colour_manual(values=other)  + guides(colour = guide_legend(override.aes = list(size=5))) + 
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Updated")


collapsed_markers_original_clusters + updated_annotation_2


updated_annotation_2 + theme_bw()
  
write_delim(zm_final_annotation_assignment.v2, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf_svd.knn_100_strict.meta.annotation_V2.txt", 
            col_names = TRUE, quote = "none", delim = "\t")


collapsed_markers_original_clusters + theme_bw()



# Cluster Investigation ---------------------------------------------------


zm_final_annotation_assignment.v2 %>% 
    ungroup() %>% 
    ggplot(., aes(V2_annotation_n,log(total))) + geom_violin() +
    geom_boxplot(width=0.1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


zm_final_annotation_assignment.v2 %>% 
  ungroup() %>% 
  ggplot(., aes(V2_annotation_n,n acrs)) + geom_violin() +
  geom_boxplot(width=0.1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


colnames(zm_final_annotation_assignment.v2)
zm_final_annotation_assignment.v2 %>% 
  ungroup() %>% 
  ggplot(., aes(acrs,log(total))) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


zm_final_annotation_assignment.v2 %>% 
  group_by(V2_annotation_n) %>% 
  summarise(total_reads = sum(total), 
            total_acrs_accessible = mean(acrs)) %>% 
  ungroup() %>% 
  ggplot(., aes(log(total_reads),total_acrs_accessible, label = V2_annotation_n)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  geom_text()


zm_collapsed_markers_original_clusters + updated_annotation_2



zm_final_annotation_assignment.v2 %>% 
  group_by(V2_annotation_n) %>% 
  summarise(total_reads = sum(total), 
            total_acrs_accessible = mean(acrs)) %>% 
  ungroup() %>% 
  ggplot(., aes(V2_annotation_n, log(total_reads))) + geom_violin() +
  geom_boxplot(width=0.1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

zm_final_annotation_assignment.v2




# Maize Annotation V3 -----------------------------------------------------
## This will be my final annotaiton of the zea mays dataset. Not much has changed at all during my analyis.
## many of the clusters have remained static and have not changed.

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


V1_V2_annotation_selected <- zm_final_annotation_assignment.v2 %>% 
    dplyr::select(cellID, V1_annotation, V2_annotation)

v3_collapsed_meta_data<- read.table(here("/Users/feilab/Projects/05.comparative_single_cell/00.data/annotation_meta/zm_v3","Zm_leaf_svd.V3_compressed_markers.meta.txt"))
colnames(v3_collapsed_meta_data)

v3_collapsed_meta_data <- as_tibble(v3_collapsed_meta_data) %>% 
  mutate(LouvainClusters_t = as.character(LouvainClusters))


combined_v2_v3_annotation <- left_join(v3_collapsed_meta_data, V1_V2_annotation_selected)
View(combined_v2_v3_annotation)


generate_v3_color_scale <- generate_color_scheme(combined_v2_v3_annotation)



collapsed_markers_original_clusters <- ggplot(combined_v2_v3_annotation, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


previous_v2_annotation <- ggplot(combined_v2_v3_annotation, aes(x=umap1, y = umap2, color = V2_annotation)) +
  generate_v3_color_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation V2")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

collapse_cell_annotation_knn_v3 <- ggplot(combined_v2_v3_annotation, aes(x=umap1, y = umap2, color = cell_annotation_knn)) +
  generate_v3_color_scale  + geom_point(size = .25, alpha = .8, show.legend = FALSE) + theme_minimal() + ggtitle("Cell Annotation Knn V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapse_cell_annotation_smooth_v3 <- ggplot(combined_v2_v3_annotation, aes(x=umap1, y = umap2, color = cell_annotation_smooth)) +
  generate_v3_color_scale  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation smooth V3")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


collapsed_markers_original_clusters + previous_v2_annotation + collapse_cell_annotation_knn_v3 + collapse_cell_annotation_smooth_v3



# Read in Sub-Cluserting Annotations --------------------------------------------------

zm_vasculature_subcluster_annotations <- (read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf_svd.V3.subclustering.vasculature.annotations.txt", header = TRUE))


'%ni%' <- Negate("%in%")
zm_epidermal_subcluster_annotations <- (read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf_svd.V3.subclustering.epidermal.annotations.txt", header = TRUE)) %>% 
    dplyr::filter(cellID %ni% zm_vasculature_subcluster_annotations$cellID)


### Bind both subclustering annnotations and then left joing them
zm.subcluster_annot <- bind_rows(zm_vasculature_subcluster_annotations, zm_epidermal_subcluster_annotations)

final.v3_annotation <- left_join(combined_v2_v3_annotation, zm.subcluster_annot,
                                 by = c("cellID")) 


zm_final_annotation_assignment.v3 <- final.v3_annotation %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("1", "8", "9", "10"),  "mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("14"),  "developing_mesophyll", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("7", "5", "16", "17"),  "bundle_sheath", V3_final_annnotation)) %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("21"), "unknown", V3_final_annnotation))  %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("19"), "xylem", V3_final_annnotation))  %>% 
  mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("3") & is.na(V3_final_annnotation) == TRUE, "vascular_parenchyma", V3_final_annnotation))  %>%
  mutate(V3_final_annnotation = if_else(is.na(V3_final_annnotation) == TRUE, "unknown", V3_final_annnotation))
  
    

  #mutate(V3_final_annnotation = if_else(LouvainClusters_t %in% c("19"),  "guard_mother_cell", V3_final_annnotation))
  
  
View(zm_final_annotation_assignment.v3)
updated.v3_annotation <- ggplot(zm_final_annotation_assignment.v3, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Final")  +
  guides(colour = guide_legend(override.aes = list(size=5)))


zm_NA_cells <- zm_final_annotation_assignment.v3 %>% 
  filter(is.na(V3_final_annnotation) == TRUE) %>% 
  ggplot(., aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
    scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Final")  +
    guides(colour = guide_legend(override.aes = list(size=5)))

original_v3_clusters <- ggplot(zm_final_annotation_assignment.v3, aes(x=umap1, y = umap2, color = LouvainClusters_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation Old")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

original_v3_clusters + updated.v3_annotation


final_cell_count <- zm_final_annotation_assignment.v3 %>% 
  group_by(V3_final_annnotation) %>% 
  summarise(total_cell_count = n())

zm_final_annotation_assignment.v3.final <- left_join(zm_final_annotation_assignment.v3, final_cell_count, by = c("V3_final_annnotation")) %>% 
  mutate(V3_annotation_n = str_c(V3_final_annnotation, total_cell_count, sep = "_ncell_")) 


#unique(zm_final_annotation_assignment.v3.final$V3_annotation_n)

row.names(zm_final_annotation_assignment.v3.final) <- zm_final_annotation_assignment.v3.final$cellID 
zm_final_annotation_assignment.v3.final %>% 
    ungroup() %>% 
    dplyr::group_by("cellID") %>% 
    dplyr::summarise(counts = n())
    

write_delim(zm_final_annotation_assignment.v3.final, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_leaf.V3_final.txt", 
            col_names = TRUE, quote = "none", delim = "\t")


updated.v3_annotation
# Generate final UMAP for Figure 1 ----------------------------------------
umap_cluster_colors <- c("#DEA940","#A2A763","#D6EEBE","#5894D4","#6AAD51","#6C8686","#FFD92F", "#7f2201","#b53310",
                         "#E7673B","#FB9A99","#E31A1C","#FDBF6F","#CAB2D6","#386CB0","#FFFFB3", "#5dead5","#c11577",
                         "#CCCCCC", "#7570B3","#B15928", "#6A3D9A", "#d87c6a", "#0060e8", "#84c5ff", "#ffd460", "#d81397")

my_pals <- show_col(umap_cluster_colors)

final_annotation_set <- unlist(unique(zm_final_annotation_assignment.v3.final$V3_final_annnotation))
names(umap_cluster_colors)  <- levels(final_annotation_set)
col_scale <- scale_colour_manual(name = "grp", values = umap_cluster_colors)

fig_1_umap_pdf <- ggplot(zm_final_annotation_assignment.v3, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  theme_half_open() +
  col_scale + geom_point(size = .25, alpha = .8) + ggtitle("Zm Cell Annotation Final")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

col_scale

ggsave("zm.cell_annot.pdf", 
       plot = fig_1_umap_pdf, 
       path = "/Users/feilab/Projects/05.comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 10, height = 10)



fig_1_umap_png <- ggplot(zm_final_annotation_assignment.v3, aes(x=umap1, y = umap2, color = V3_final_annnotation)) +
  theme_void() + 
  col_scale + geom_point(size = .5, alpha = .8, show.legend = FALSE) + 
  guides(colour = guide_legend(override.aes = list(size=5)))

fig_1_umap_png

ggsave("zm.cell_annot.png", 
       plot = fig_1_umap_png, 
       path = "/Users/feilab/Projects/05.comparative_single_cell/comparative_single_cell_imgs/figures/figure_1/UMAPs",
       units = "in",
       width = 5, height = 5)
