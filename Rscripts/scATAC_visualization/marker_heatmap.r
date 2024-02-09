# load arguments
args <- commandArgs(T)
#if(length(args)!=5){stop("Rscript normGBA.R <gene.sparse> <meta> <Zea_mays.AGPv4.36.Allgene.nuclear.bed> <prefix> <F>")}
input_data <- as.character(args[1])
meta <- as.character(args[2])
gene <- as.character(args[3])
prefix <- as.character(args[4])

# load libraries
library("heatmaply")
library(plotly)
library(dplyr)
library(edgeR)
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(tidyverse)
library(ComplexHeatmap)

#gene_bodysorghum_bicolor_tis_leaf_nmf_step_2_knn_27.out.de_novo.rds
#input <- here("/home/jpm73279/r_script_dev/heat_map_plotting","zea_mays_tis_leaf_svd_knn_100_strict_step_2_knn_27.out.rds")
#input_2 <- here("/home/jpm73279/r_script_dev/lw_plotting","zea_mays.normalized_gene_acc_scores.leaf_svd_knn_100_strict.GBaccessibility.sparse")
#meta <- here("/home/jpm73279/r_script_dev/heat_map_plotting","Zm_leaf.merged_replicates.knn_100.stict_peak_filtering.SVD.full.metadata.txt")
#gene <- here("/home/jpm73279/r_script_dev/heat_map_plotting","ZM.Markers.unique_counter.leaf.final.bed")
#gene_DA <- here("/home/jpm73279/r_script_dev/lw_plotting","sorghum_bicolor.tis_leaf_nmf.cluster.DA_genes.merged.bed")
#prefix <- "TEST_SORGHUM_TEST"

meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)  
gene_markers <- gene_markers  %>% 
    arrange(type)

all_markers <- gene_markers$geneID

raw_cpm_counts_all_genes <- read_delim(input_data, delim=" ", col_names = c("gene_name", "barcode", "accessability")) %>%
    dplyr::mutate(cellID = barcode)  %>% 
    dplyr::mutate(geneID = gene_name)



merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))  %>% 
    mutate(safe_cluster_name = str_c("Louvain_C", LouvainClusters, sep ="_"))  %>% 
    dplyr::select(-LouvainClusters)  %>% 
    group_by(safe_cluster_name, geneID)  %>% 
    summarise(counts = sum(accessability, na.rm = TRUE))

### Alt CPM Calc
merged_meta_cpm_information_copied <- merged_meta_cpm_information 
catch <- merged_meta_cpm_information_copied  %>% 
    group_by(safe_cluster_name) %>% 
    group_map(~(cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>% 
    unlist()

caught_values <- as_tibble(catch)
merged_meta_cpm_information_copied$grouped_CPM <- caught_values

merged_meta_cpm_information_markers <- left_join(merged_meta_cpm_information_copied, gene_markers, by = c("geneID"))  %>% 
    filter(geneID %in% gene_markers$geneID)

merged_meta_cpm_information_markers_before_matrix <- dplyr::select(merged_meta_cpm_information_markers, safe_cluster_name, name, grouped_CPM)


generate_matrix <- as.matrix(pivot_wider(merged_meta_cpm_information_markers_before_matrix, names_from = "name", values_from = grouped_CPM))
generate_matrix[is.na(generate_matrix)] <- 0
generate_matrix[is.nan(generate_matrix)] <- 0

generate_matrix_2 <- generate_matrix[,-1]
rownames(generate_matrix_2) <- generate_matrix[,1]
class(generate_matrix_2) <- "numeric"

test <- scale(generate_matrix_2)
test[is.na(test)] <- 0
test[is.nan(test)] <- 0


#mtscaled <- as.matrix(scale(quick))
grab_map <- heatmaply(scale(generate_matrix_2),  
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red"))

output_name <- paste0(prefix,".html")

htmlwidgets::saveWidget(as_widget(grab_map), output_name)
