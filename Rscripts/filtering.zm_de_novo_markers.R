library(tidyverse)

read_de_novo_files <- function(file_paths) {
  # Use map_dfr from purrr to read each file and bind them into one data frame
  all_data <- map_dfr(file_paths, function(file_path) {
    # Extract cell-type from the file name
    cell_type <- str_extract(file_path, "(?<=de_novo\\.)[\\w_]+(?=\\.upregulated)")
    
    # Read the file
    data <- read_tsv(file_path)
    
    # Add a column for cell-type and file name
    data <- data %>% 
      mutate(cell_type = cell_type,
             file_name = basename(file_path))
    
    return(data)
  })
  
  return(all_data)
}

# Define the full paths of the files
file_paths <- c(
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.bundle_sheath.upregulated_genes.deseq2_output.tsv",
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.companion_cells_sieve_elements.upregulated_genes.deseq2_output.tsv",
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.epidermis.upregulated_genes.deseq2_output.tsv",
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.mesophyll.upregulated_genes.deseq2_output.tsv",
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.procambial_meristem.upregulated_genes.deseq2_output.tsv",
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.protoderm.upregulated_genes.deseq2_output.tsv"
)

`%nin%` = Negate(`%in%`)
# Read the files and combine the data
data <- read_de_novo_files(file_paths)
colnames(data)


filtered_de_novo_markers <- data %>% 
    dplyr::filter(log2FoldChange > 3 & lfcSE < .6) %>% 
    select(-file_name)

glimpse(filtered_de_novo_markers)
    
filtered_de_novo_markers_counts <- filtered_de_novo_markers %>% 
    group_by(gene_name) %>% 
    dplyr::summarise(counts = n()) %>% 
    dplyr::filter(counts > 1)

#filtered_de_novo_markers %>% 
#    dplyr::filter(gene_name %in% filtered_de_novo_markers_counts$gene_name) %>% 
#    View()



read_bed_files <- function(file_paths) {
  library(purrr)
  library(readr)
  library(stringr)
  library(dplyr)
  
  all_data <- map_dfr(file_paths, function(file_path) {
    # Extract marker type or cell type from the file name
    marker_type <- str_extract(file_path, "(?<=Zm\\.de_novo\\.)[^\\.]+")
    
    # Read the BED file with explicit column names
    data <- read_tsv(file_path, col_names = c("chr", "start", "end", "geneID", "name", "type"), col_types = cols(.default = "c"))
    
    # Add a column for marker/cell type and file name
    data <- data %>%
      mutate(marker_type = marker_type,
             file_name = basename(file_path))
    
    return(data)
  })
  
  return(all_data)
}

# Define the full paths of the BED files
bed_files <- c(
  "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/de_novo_markers/Zm.de_novo.all_combined_cell_types_upregulated.bed"
)

# Read the files and combine the data
bed_files_markers <- read_bed_files(bed_files)

allowed_cell_types <- c("bundle_sheath,procambial_meristem", 
                        "epidermis,protoderm",
                        "companion_cells_sieve_elements,bundle_sheath",
                        "companion_cells_sieve_elements",
                        "epidermis",
                        "protoderm",
                        "mesophyll",
                        "bundle_sheath",
                        "procambial_meristem"
)
View(bed_files_markers)
# View the first few rows of the combined data frame
filtered_bed_markers <- bed_files_markers %>% 
    dplyr::filter(geneID %in% filtered_de_novo_markers$gene_name) %>% 
    dplyr::filter(type %in% allowed_cell_types) %>% 
    dplyr::select(-file_name)

#colnames(filtered_de_novo_markers)
final_de_novo_tsv_save <- filtered_de_novo_markers %>% 
    dplyr::filter(gene_name %in% filtered_bed_markers$geneID) %>% 
    dplyr::distinct()

View(filtered_bed_markers)

write_tsv(final_de_novo_tsv_save, file = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/Zm.de_novo_markers.final.tsv")
write_tsv(filtered_bed_markers, file = "/Users/pablomendieta/Projects/comparative_single_cell/comparative_single_cell_imgs/tables/Zm.de_novo_markers.final.bed")
