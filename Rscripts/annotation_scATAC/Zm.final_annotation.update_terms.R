library(tidyverse)
library(here)

zm_v4_file_location <- here("Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_annot_v4", "Zm.leaf_annot.V4.meta.final.txt")
zm_v4_annotation <- read_delim(zm_v4_file_location, col_names = TRUE, delim = "\t")

unique(zm_v4_annotation$Zm_v4_annot)


dropping_annotations <- c("ground_meristem", "subsidiary_mother_cells",
                          "stomatal_precursor", "subsidiary_cells", "unknown")

'%ni%' <- Negate("%in%")
filtered_cell_types <- zm_v4_annotation %>% 
  dplyr::filter(FRiP > .2) %>% 
  dplyr::filter(Zm_v4_annot %ni% dropping_annotations) %>% 
  dplyr::mutate(final_annotation_safe = case_when(Zm_v4_annot == "developing_mesophyll" ~ "mesophyll_desc_developing",
                                             Zm_v4_annot == "sieve_elements" ~ "phloem_sieve_elements",
                                             Zm_v4_annot ==  "procambium" ~ "procambial_meristem",
                                             TRUE ~ Zm_v4_annot)) %>% 
  dplyr::mutate(final_annotation = case_when(Zm_v4_annot == "developing_mesophyll" ~ "mesophyll;developing",
                                     Zm_v4_annot == "sieve_elements" ~ "phloem_sieve_elements",
                                     Zm_v4_annot ==  "procambium" ~ "procambial_meristem",
                                     TRUE ~ Zm_v4_annot)) %>% 
  dplyr::mutate(reduce_resolution_annotation = case_when(final_annotation == "mesophyll;developing" ~ "mesophyll",
                                                         final_annotation == "sieve_elements" ~ "companion_cells_sieve_elements",
                                                         final_annotation == "companion_cells" ~ "companion_cells_sieve_elements",
                                                         final_annotation == "phloem_sieve_elements" ~ "companion_cells_sieve_elements",
                                                        TRUE ~ final_annotation)) %>% 
  dplyr::select(-total_cell_count)
  


filtered_cell_types_count_n <- filtered_cell_types %>%
    group_by(final_annotation_safe) %>% 
    summarise(count_cells = n())

finalized_annotation <- left_join(filtered_cell_types, filtered_cell_types_count_n, by = c("final_annotation_safe")) %>% 
    mutate(final_annotation_n = str_c(final_annotation_safe, count_cells, sep = "_n_cell_")) %>% 
    dplyr::select(-count_cells)


unique(finalized_annotation$reduce_resolution_annotation)

write_delim(finalized_annotation, 
            file = "Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_annot_final/Zm.leaf_annot.V5.meta.frozen.txt",
            delim = "\t", 
            col_names = TRUE )


## Write these outputs to test switching of cell-type specific ACRs in a single genome
finalized_annotation %>% 
  group_by(sampleID) %>% 
  group_walk(~ write_delim(.x, paste0("Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_annot_final/Zm.leaf_annot.V5.meta.frozen.", .y$sampleID, ".txt"), delim = "\t", col_names = TRUE))


unique(finalized_annotation$final_annotation_n)
