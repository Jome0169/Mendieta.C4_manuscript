library(tidyverse)
library(here)

sb_v4_file_location <- here("Mendieta_et_al_comparative_single_cell/metrics/annotations/sb_annot_v4", "Sb.leaf_annot.V4.meta.final.2022-11-14.txt")
sb_v4_annotation <- read_delim(sb_v4_file_location, col_names = TRUE, delim = "\t")




dropping_annotations <- c("phloem", "unknown.7", "proto_xylem")

'%ni%' <- Negate("%in%")
filtered_cell_types <- sb_v4_annotation %>% 
  dplyr::filter(FRiP > .2) %>% 
  dplyr::filter(sb_v4_annot %ni% dropping_annotations) %>% 
  dplyr::mutate(final_annotation_safe = case_when(sb_v4_annot == "developing_mesophyll" ~ "mesophyll_desc_developing",
                                             sb_v4_annot == "ground_meristem" ~ "unknown_desc_1",
                                             sb_v4_annot ==  "procambium" ~ "procambial_meristem",
                                             TRUE ~ sb_v4_annot)) %>% 
  
  dplyr::mutate(final_annotation = case_when(sb_v4_annot == "developing_mesophyll" ~ "mesophyll;developing",
                                             sb_v4_annot == "ground_meristem" ~ "unknown;1",
                                             sb_v4_annot ==  "procambium" ~ "procambial_meristem",
                                             TRUE ~ sb_v4_annot)) %>% 
  dplyr::mutate(reduce_resolution_annotation = case_when(final_annotation == "mesophyll;developing" ~ "mesophyll",
                                                         final_annotation == "unknown;1" ~ "unknown",
                                                         TRUE ~ final_annotation)) %>% 
  dplyr::select(-total_cell_count)

unique(filtered_cell_types$reduce_resolution_annotation)

filtered_cell_types_count_n <- filtered_cell_types %>%
  group_by(final_annotation_safe) %>% 
  summarise(count_cells = n())

finalized_annotation <- left_join(filtered_cell_types, filtered_cell_types_count_n, by = c("final_annotation_safe")) %>% 
  mutate(final_annotation_n = str_c(final_annotation_safe, count_cells, sep = "_n_cell_")) %>% 
  dplyr::select(-count_cells)


unique(finalized_annotation$reduce_resolution_annotation)


write_delim(finalized_annotation, 
            file = "Mendieta_et_al_comparative_single_cell/metrics/annotations/Sb_annot_final/Sb.leaf_annot.V5.meta.frozen.txt",
            delim = "\t", 
            col_names = TRUE )

## Write these outputs to test switching of cell-type specific ACRs in a single genome
finalized_annotation %>% 
  group_by(sampleID) %>% 
  group_walk(~ write_delim(.x, paste0("Mendieta_et_al_comparative_single_cell/metrics/annotations/Sb_annot_final/Sb.leaf_annot.V5.meta.frozen.", .y$sampleID, ".txt"), delim = "\t", col_names = TRUE))




