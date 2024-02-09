haidong_os_annotation <- read_delim("/Users/pablomendieta/Desktop/opt_final_comb_annot_ver8.2.txt")

haidong_os_annotation_leaf <- haidong_os_annotation %>% 
  dplyr::filter(grepl("leaf1$|leaf2$", cellID)) %>% 
  dplyr::mutate(cellID = str_replace(cellID, "leaf1$", "1-Os.leaf1")) %>%
  dplyr::mutate(cellID = str_replace(cellID, "leaf2$", "1-Os.leaf2"))


haigong_annotation_leaf_final <- haidong_os_annotation_leaf %>% 
  dplyr::select(cellID, Final_annotation_TCP_up, umap1, umap2) %>% 
  dplyr::rename("haidong_umap1" = umap1) %>% 
  dplyr::rename("haidong_umap2" = umap2) %>% 
  dplyr::mutate(Final_annotation_TCP_up = gsub("leaf\\.", "", Final_annotation_TCP_up),
                Final_annotation_TCP_up = tolower(Final_annotation_TCP_up))

colnames(haidong_os_annotation_leaf)

os.combined_annotations <- full_join(os_subcluster_leaf_annotated, haigong_annotation_leaf_final, by = "cellID")


all_values <- unique(c(os.combined_annotations$annotation_v1, os.combined_annotations$Final_annotation_TCP_up))
color_vec <- c("companion_cell"  = "#1b9e77", 
               "mesophyll"       = "#d95f02", 
               "epidermis"       = "#80b1d3",
               "bundle_sheath"   = "#e7298a",
               "protoderm"       = "#66a61e",
               "subsidiary_cell" = "#e6ab02",
               "vasculature"     = "#a6761d",
               "NA"              = "#666666", # using grey for NA values
               "leaf_un2"        = "#1c9099",
               "bulliform"       = "#756bb1",
               "vascular_cell"   = "#995f0e",
               "guard_cell"      = "#7fc97f")

pm_annot <- ggplot(os.combined_annotations, aes(x=umap1, y = umap2, color = as.factor(annotation_v1))) +
  geom_point(size = .25, alpha = .8) +
  scale_color_manual(values = color_vec) +
  theme_minimal() +
  ggtitle("Pm Cell Annotation") +
  guides(colour = guide_legend(override.aes = list(size=5)))

haidoong_annot <- ggplot(os.combined_annotations, aes(x=umap1, y = umap2, color = as.factor(Final_annotation_TCP_up))) +
  geom_point(size = .25, alpha = .8) +
  scale_color_manual(values = color_vec) +
  theme_minimal() +
  ggtitle("Cell Annotation SubClsuter") +
  guides(colour = guide_legend(override.aes = list(size=5)))


library(patchwork)
pm_annot + haidoong_annot



pm_annot_haidong_umap <- ggplot(os.combined_annotations, aes(x=haidong_umap1, y = haidong_umap2, color = as.factor(annotation_v1))) +
  geom_point(size = .25, alpha = .8) +
  scale_color_manual(values = color_vec) +
  theme_minimal() +
  ggtitle("Pm Cell Annotation") +
  guides(colour = guide_legend(override.aes = list(size=5)))

haidoong_annot_haidong_umap <- ggplot(os.combined_annotations, aes(x=haidong_umap1, y = haidong_umap2, color = as.factor(Final_annotation_TCP_up))) +
  geom_point(size = .25, alpha = .8) +
  scale_color_manual(values = color_vec) +
  theme_minimal() +
  ggtitle("Cell Annotation SubClsuter") +
  guides(colour = guide_legend(override.aes = list(size=5)))




pm_annot_haidong_umap + haidoong_annot_haidong_umap





final <- os.combined_annotations %>% 
  group_by(annotation_v1, Final_annotation_TCP_up) %>% 
  summarise(counts = n()) %>% 
  mutate(score = case_when(annotation_v1 == Final_annotation_TCP_up ~ counts,
                           annotation_v1 != Final_annotation_TCP_up ~ 0))


dim(haigong_annotation_leaf_final)
