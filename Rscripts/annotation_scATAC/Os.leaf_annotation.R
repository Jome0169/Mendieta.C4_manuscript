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



os_subcluster_meta_file <- here("/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/os_annot/Os_leaf.merged_replicates.reclustered.SVD.full.metadata.txt")
os_subcluster_leaf <- read.table(os_subcluster_meta_file)




annotation_mapping <- c(
  "Lc_1" = "epidermis",
  "Lc_2" = "companion_cell",
  "Lc_3" = "mesophyll",
  "Lc_4" = "mesophyll",
  "Lc_5" = "mesophyll",
  "Lc_6" = "companion_cell",
  "Lc_7" = "protoderm",
  "Lc_8" = "companion_cell",
  "Lc_9" = "mesophyll",
  "Lc_10" = "bundle_sheath",
  "Lc_11" = "mesophyll",
  "Lc_12" = "bundle_sheath",
  "Lc_13" = "protoderm",
  "Lc_14" = "mesophyll",
  "Lc_15" = "unknown_cells_2",
  "Lc_16" = "unknown_cells_1"
)


# Filter out unwanted subclusters and add annotation
os_subcluster_leaf_annotated <- os_subcluster_leaf %>%
  mutate(subcluster_lc_safe = str_c("Lc", LouvainClusters, sep = "_")) %>%
  mutate(annotation_v1 = annotation_mapping[subcluster_lc_safe])

head(os_subcluster_leaf_annotated)
ggplot(os_subcluster_leaf_annotated, aes(x=umap1, y = umap2, color = as.factor(annotation_v1))) +
  geom_point(size = .25, alpha = .8) + theme_minimal() + ggtitle("Cell Annotation SubClsuter")  +
  guides(colour = guide_legend(override.aes = list(size=5)))

View(os_subcluster_leaf_annotated)


write_delim(os_subcluster_leaf_annotated, "/Users/pablomendieta/Projects/comparative_single_cell/Mendieta_et_al_comparative_single_cell/metrics/annotations/os_annot/os.leaf_annotation.V1.meta.txt", 
            col_names = TRUE, quote = "none", delim = "\t")



