library("here")
library(devtools)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(RColorBrewer)


nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)


sb.meta_raw <- read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb_leaf.merged_replicates.QC_V2.SVD.full.metadata.txt")
sb.subcluster_results <- read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb_leaf.merged_replicates.QC.subclustering_results.tsv", header = TRUE)

?read.table
head(sb.meta_raw)
head(sb.subcluster_results)
sb.combined_subclustering_results <- left_join(sb.meta_raw, sb.subcluster_results, by =c('cellID'))


subcluster_results <- ggplot(sb.combined_subclustering_results, aes(x=umap1, y = umap2, color = sb_v4_annot)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Current Annotation")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

lc_clusters_global <- ggplot(sb.combined_subclustering_results, aes(x=umap1, y = umap2, color = LouvainCluster_t)) +
   scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Louvain Clustering on ACRs")  +
   guides(colour = guide_legend(override.aes = list(size=5)))
 
 
lc_clusters_global + subcluster_results


sb.combined_subclustering_results.final <- sb.combined_subclustering_results %>% 
  mutate(sb_v4_annot = if_else(LouvainCluster_t %in% c("Louvain_c4", "Louvain_c8", "Louvain_c7"), "mesophyll", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(LouvainCluster_t %in% c("Louvain_c12"), "developing_mesophyll", sb_v4_annot))  %>% 
  mutate(sb_v4_annot = if_else(LouvainCluster_t %in% c("Louvain_c10"), "bundle_sheath", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(LouvainCluster_t %in% c("Louvain_c9","Louvain_c14"), "mesophyll", sb_v4_annot)) 



ggplot(sb.combined_subclustering_results.final, aes(x=umap1, y = umap2, color = sb_v4_annot)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Louvain Clustering on ACRs")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



# Annotate Remaining NAs filtered during sub-clustering -------------------

nas_left <- sb.combined_subclustering_results.final %>% 
  filter(is.na(sb_v4_annot) == TRUE) %>% 
  ggplot(., aes(x=umap1, y = umap2, color = sb_v4_annot)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Current Annotation")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

nas_lef_LC  <- sb.combined_subclustering_results.final %>% 
  filter(is.na(sb_v4_annot) == TRUE) %>% 
  ggplot(., aes(x=umap1, y = umap2, color = LouvainCluster_t)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Louvain Clustering on ACRs")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

nas_left + nas_lef_LC


## Select down the cells to pull 
sb.v4.PCAs <- read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb_leaf.merged_replicates.QC_V2.SVD.full.reduced_dimensions.txt")

sb.na_annots <- sb.combined_subclustering_results.final %>% 
  filter(is.na(sb_v4_annot) == TRUE)
sb.annot.v4 <- sb.combined_subclustering_results.final %>% 
  filter(is.na(sb_v4_annot) != TRUE)

annotate_PCAs <- sb.v4.PCAs[sb.annot.v4$cellID,]
non_annotated_PCAs <- sb.v4.PCAs[sb.na_annots$cellID,]


knn_pred = class::knn(annotate_PCAs, non_annotated_PCAs, sb.annot.v4$sb_v4_annot, k = 5, prob = TRUE)
knn_prob = attributes(knn_pred)$prob
sb.na_annots$sb_v4_annot <- knn_pred


sb.v4.final_annot_meta <- bind_rows(sb.na_annots, sb.annot.v4)


cluster_numbers <- sb.v4.final_annot_meta %>% 
  group_by(sb_v4_annot) %>% 
  summarise(total_cell_count = n())


sb_v4_annot.cell_counts <- left_join(sb.v4.final_annot_meta, cluster_numbers, by = c("sb_v4_annot")) %>% 
  mutate(v4_annotation_n = str_c(sb_v4_annot, total_cell_count, sep = "_ncell_"))


ggplot(sb_v4_annot.cell_counts, aes(x=umap1, y = umap2, color = v4_annotation_n)) +
  scale_colour_manual(values=other)  + geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Sb V4 Annotation")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

write_delim(sb_v4_annot.cell_counts, "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb.leaf_annot.V4.meta.final.txt", 
            col_names = TRUE, quote = "none", delim = "\t")


