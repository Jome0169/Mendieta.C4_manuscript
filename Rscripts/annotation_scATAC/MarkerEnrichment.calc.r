# load libraries
library(dplyr)
library("edgeR")
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(tidyverse)
library(ComplexHeatmap)
library(here)
library(modelr)
library(preprocessCore)
library(parallel)

# load arguments
args <- commandArgs(trailingOnly=T)

#zm_working_dir <- "/scratch/jpm73279/comparative_single_cell/08.annotation_figures/zea_mays"
zm_input <- as.character(args[1])
zm_meta <- as.character(args[2])
zm_gene <- as.character(args[3])
column_name <- as.character(args[4])
output_name <- as.character(args[5])

generate_sparse_matrix <- function(raw_counts_file){
    
    zm_raw_cpm_counts_all_genes <- read_delim(zm_input, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
        dplyr::mutate(cellID = barcode)  %>%
        dplyr::mutate(geneID = gene_name)

    
    combined_large_w_sparse <- zm_raw_cpm_counts_all_genes  %>% 
        filter(gene_name != "Annotation")  %>% 
        dplyr::select(gene_name, barcode, accessability)  %>% 
        dplyr::filter(barcode %in% zm_meta_data$cellID)
    

    combined_large_w_sparse <- combined_large_w_sparse  %>% 
        dplyr::select(gene_name, barcode, accessability)  %>% 
        mutate(across(accessability, as.numeric))
    
    
    gene_names <- unique(combined_large_w_sparse$gene_name)
    barcodes <- unique(combined_large_w_sparse$barcode)

    combined_large_w_sparse$row <- match(combined_large_w_sparse$gene_name, gene_names)
    combined_large_w_sparse$col <- match(combined_large_w_sparse$barcode, barcodes)


    UIMatrix <- sparseMatrix(i = combined_large_w_sparse$row,
        j = combined_large_w_sparse$col,
        x = combined_large_w_sparse$accessability,
        dimnames=list(gene_names, barcodes))
    
    return(UIMatrix)
}

random_sample_calculate_mean_sd <- function(x, zm_meta_all, n, chosen_cell_ids, sparse_matrix){
    
    `%ni%` <- Negate(`%in%`)
    non_subset_cells <- zm_meta_all  %>% 
        dplyr::filter(cellID %ni% chosen_cell_ids)
    
    
    random_cells <- sample_n(non_subset_cells, n, replace = FALSE)
    
    random.cell.mat <- sparse_matrix[, random_cells$cellID]
    random.cell.mat.normalized <- random.cell.mat %*% diag(1e4/colSums(random.cell.mat))
    
    random.cell.mat.log_mean <- (rowMeans(random.cell.mat.normalized))
    random_cells.marker_values.log_std_dev <- apply(random.cell.mat.normalized, 1, sd)
    random_cells.merged_values <- rbind(random.cell.mat.log_mean, random_cells.marker_values.log_std_dev)

    return(random_cells.merged_values)
}

louvain_cluster_permutaion <- function(x, meta_data, n_permutations, col_name, sparse_matrix_markers) {
    
    col_var <- c(col_name)
    print(col_var)
    selected.cells.mat <- meta_data  %>% 
        mutate(across(!!sym(col_var), as.character))  %>% 
        dplyr::filter((!!sym(col_var) == (x)))  %>% 
        dplyr::filter(cellID %in% colnames(sparse_matrix_markers))
    
    print(dim(selected.cells.mat))

    selected_cells.sparse <- sparse_matrix_markers[, selected.cells.mat$cellID]
    selected_cells.sparse.normalized <- selected_cells.sparse %*% diag(1e4/colSums(selected_cells.sparse))
    selected_cells.sparse.mean <- (rowMeans(selected_cells.sparse.normalized))
    #bs.marker_values.log_std_dev <- (apply(bs.marker_values.mat.normalized, 1, sd))
    selected.cells.marker_mean <- rbind(selected_cells.sparse.mean)
    
    #Generate permutations. This could be altered in the future as well to calculate actual P values.
    #But for now relative enrichment values works.
    x <- seq(1:1000)
    permutations.test <- mclapply(x, random_sample_calculate_mean_sd, meta_data, dim(selected.cells.mat)[1], selected.cells.mat$cellID, sparse_matrix_markers, mc.cores = 15)
    permutations.test <- setNames(permutations.test,x)
    permutations.test.output <- as.data.frame(do.call(rbind, permutations.test))
    
    
    z <- permutations.test.output  %>% 
        t()  %>% 
        as_tibble()  %>% 
        dplyr::select(contains("log_mean"))  %>% 
        t()

    mean_vals <- colMeans(z)

    q <- permutations.test.output  %>% 
        t()  %>% 
        as_tibble()  %>% 
        dplyr::select(contains("log_std"))  %>% 
        t()

    std_vals <- colMeans(q)

    enrich_values <- (selected.cells.marker_mean[1,] - mean_vals)/std_vals

    return(enrich_values)
}

#zm_working_dir <- "/scratch/jpm73279/comparative_single_cell/08.annotation_figures/zea_mays"
#gene_bodysorghum_bicolor_tis_leaf_nmf_step_2_knn_27.out.de_novo.rds
#zm_input <- here(zm_working_dir,"00.data/zea_mays.gene_body_acc_zea_may.v4_annot.counts.txt")
#input_2 <- here("/home/jpm73279/r_script_dev/lw_plotting","zea_mays.normalized_gene_acc_scores.leaf_svd_knn_100_strict.GBaccessibility.sparse")
#zm_meta <- here(zm_working_dir,"00.data/Zm.leaf_annot.V4.meta.final.txt")
#zm_gene <- here(zm_working_dir,"00.data/Zm.markers.leaf.txt")
#zm_gene_DA <- here(zm_working_dir,"00.data/Zm-B73-REFERENCE-NAM_Zm00001eb.1.genes.bed")

#Read in input values
zm_meta_data <- read.delim(zm_meta)
zm_gene_markers <- read.delim(zm_gene)
zm_gene_markers <- zm_gene_markers  %>%
    arrange(type)

zm_all_markers <- zm_gene_markers$geneID

#Read in the matrix
UIMatrix <- generate_sparse_matrix(zm_input)

#Filter the matrix as well as the markers for those present in the matrix
zm_gene_markers.filtered <- zm_gene_markers  %>% 
    dplyr::filter(geneID %in% rownames(UIMatrix))

only_marker_matrix <- UIMatrix[zm_gene_markers.filtered$geneID, ]


#column_name <- "LouvainClusters"
#Testing the above function on a single Louvain cluster
#x <- louvain_cluster_permutaion("1", zm_meta_data, n=1000, column_name, only_marker_matrix)


col_unique_values <- as.character(unique(zm_meta_data[[column_name]]))

louvain_cluster_enrichment <- lapply(col_unique_values, louvain_cluster_permutaion, zm_meta_data, n=1000, column_name, only_marker_matrix)

names(louvain_cluster_enrichment) <- col_unique_values
names(louvain_cluster_enrichment)

marker_names <- as.vector(zm_gene_markers.filtered$name)
marker_gene_IDs <- as.vector(zm_gene_markers.filtered$geneID)


flattened_louvain_enrichment <- bind_rows(louvain_cluster_enrichment, .id = 'cluster')

#For testing.
#output_name <- "marker_enrichment"

all_marker_output <- str_c(output_name, ".all_marker_enrichment.tsv")
write_delim(flattened_louvain_enrichment, file = all_marker_output, delim = "\t", quote = c("none"))

top_10_enrich <- str_c(output_name, ".top_10_markers_enrichment.tsv")
pivoted_names <- pivot_longer(flattened_louvain_enrichment, cols=-cluster, names_to = "geneID", values_to = "enrichment")
                       
top_10_markers_per_cluster <- pivoted_names  %>% 
    group_by(cluster)  %>% 
    arrange(cluster, desc(enrichment))  %>% 
    top_n(10)  %>% 
    left_join(., zm_gene_markers.filtered)  %>% 
    dplyr::select(-chr, -start, -end) 
                           
write_delim(top_10_markers_per_cluster, file = top_10_enrich, delim = "\t", quote = c("none"))
