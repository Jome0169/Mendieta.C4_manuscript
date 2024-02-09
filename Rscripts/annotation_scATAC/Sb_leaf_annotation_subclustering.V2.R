library("here")
library(devtools)
library(Seurat)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(RColorBrewer)

load_all('/Users/feilab/Programming/Socrates')



# Vasculature SubClustering Annotation ------------------------------------


sb_vasc_annot <- read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb_leaf.merged_replicates.QC_V2.subclustering_vasculature.SVD.full.metadata.txt") %>% 
    mutate(louvain_clusters_sub = str_c("LC", as.character(LouvainClusters), sep = "_"))

nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)



sb.original_vasc_LCs <- ggplot(sb_vasc_annot, aes(x=umap1, y = umap2, color = louvain_clusters_sub)) +
  scale_colour_manual(values=other)  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

sb.original_vasc_LCs


View(sb_vasc_annot)


sb_vasc_annot.updated_annot <- sb_vasc_annot %>% 
  mutate(sb_v4_annot = "NA") %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_3", "LC_5", "LC_12", "LC_13",
                                                          "LC_10","LC_11"),  "bundle_sheath", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_4"), "companion_cell", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_2"), "xylem", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_1"), "mesophyll", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_7", "LC_9","LC_8", "LC_6"), "procambial_meristem", sb_v4_annot))


ggplot(sb_vasc_annot.updated_annot, aes(x=umap1, y = umap2, color = sb_v4_annot)) +
  scale_colour_manual(values=other)  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))




# Epidermal Sub-clustering annotation -------------------------------------



sb_epidermal_annot <- read.table("/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb_leaf.merged_replicates.QC_V2.subclustering_epidermal.SVD.full.metadata.txt") %>% 
  mutate(louvain_clusters_sub = str_c("LC", as.character(LouvainClusters), sep = "_"))

nColor <- 40
other <- randomcoloR::distinctColorPalette(k = 55)


sb_epidermal_annot_LCs <- ggplot(sb_epidermal_annot, aes(x=umap1, y = umap2, color = louvain_clusters_sub)) +
  scale_colour_manual(values=other)  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))
sb_epidermal_annot_LCs


sb_epidermal_annot.updated_annot <- sb_epidermal_annot %>% 
  mutate(sb_v4_annot = "NA") %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_1", "LC_2", "LC_4", "LC_10",
                                                           "LC_10","LC_11"),  "epidermis", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_5", "LC_6", "LC_7", "LC_3",
                                                           "LC_11","LC_13"),  "protoderm", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_8"),  "ground_meristem", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_9"),  "subsidary_cells", sb_v4_annot)) %>% 
  mutate(sb_v4_annot = if_else(louvain_clusters_sub %in% c("LC_12"),"guard_cells", sb_v4_annot))


ggplot(sb_epidermal_annot.updated_annot, aes(x=umap1, y = umap2, color = sb_v4_annot)) +
  scale_colour_manual(values=other)  + geom_point(size = .6, alpha = .8) + theme_minimal() + ggtitle("LouvainClusters")  +
  guides(colour = guide_legend(override.aes = list(size=5)))



sb.v4_updated_annotation.sub_cluster <- bind_rows(sb_epidermal_annot.updated_annot, sb_vasc_annot.updated_annot)

sb.v4_updated_annotation.sub_cluster.annot_ony <- sb.v4_updated_annotation.sub_cluster %>% 
    select(cellID, louvain_clusters_sub, sb_v4_annot)

#rownames(sb.v4_updated_annotation.sub_cluster.annot_ony) <- NULL


write_tsv(sb.v4_updated_annotation.sub_cluster.annot_ony, 
          file = "/Users/feilab/Projects/05.comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4/Sb_leaf.merged_replicates.QC.subclustering_results.tsv",
          col_names = TRUE, quote = "none")
