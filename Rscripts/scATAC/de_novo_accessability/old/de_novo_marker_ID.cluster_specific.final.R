library(dplyr)
library(here)
library(DESeq2)
library(tidyverse)
library(rlang)

# load arguments
args <- commandArgs(trailingOnly=T)

#args    
meta_data_file <- as.character(args[1])
gene_accessability_file <- as.character(args[2])
marker_gene_file <- as.character(args[3])
all_genes_bed <- as.character(args[4])
column_name_inp <- as.character(args[5])
species <- as.character(args[6])
output_base <- as.character(args[7])
output_locaiton <- as.character(args[8])


#Local Development Files
#data_path <- "/Users/feilab/Projects/05.comparative_single_cell/00.data/test_data_log2fc"
#meta_data_file <- here(data_path, "sb_leaf_nmf_compressed_markers.meta.txt")
#gene_accessability_file <- here(data_path, "sorghum_bicolor.normalized_gene_acc_scores.leaf_nmf.sctGBAcounts.sparse")
#1marker_gene_file <- here(data_path, "Sb_leaf.maize_markers.ortho.visualize.bed")


#Load Marker Genes

message("Loading Meta Data")
marker_genes <- read.table(marker_gene_file, header =TRUE)
#Load Meta Data

message("LOading markers Genes")
loaded_meta_data <- read.table(meta_data_file, header = TRUE) #%>% 
  #mutate(Louvain_cluster_name = str_c("LouvainC", LouvainClusters_t, sep = "_"))


if (species == "sorghum_bicolor") {
  header_names <- c("chr", "start", "end", "geneID", "val", "strand")
  #all_gene_File <- here(data_path, "Sbicolorv5.1.gene.bed")
  all_genes <- read.table(all_genes_bed, header =FALSE) %>% 
      dplyr::select(V1:V6) %>% 
      setNames(header_names) %>%  
      mutate(geneID = str_remove_all(geneID,".g")) 
  #Load Sparse Matrix
  loaded_sparse_matric <- read_delim(gene_accessability_file, delim='\t', col_names = c("gene_name", "barcode", "accessability")) %>% 
      mutate(gene_name = str_remove_all(gene_name,".g"))
  
} else {
    
  header_names <- c("chr", "start", "end", "geneID", "val", "strand")
  #all_gene_File <- here(data_path, "Sbicolorv5.1.gene.bed")
  all_genes <- read.table(all_genes_bed, header =FALSE) %>% 
    dplyr::select(V1:V6) %>% 
    setNames(header_names)
  #Load Sparse Matrix
  loaded_sparse_matric <- read_delim(gene_accessability_file, delim='\t', col_names = c("gene_name", "barcode", "accessability"))
  
  }
    

  

# Generate Function  ------------------------------------------------------
get_downregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-.1)], rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}
get_upregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange>=.1)], rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}
generate_null_sample <- function(meta_data, sparse_matrix, slot_name, subsampled_cluster) {
  
  slot_var <- c(slot_name)
    
  subsampled_cluster_metrics <- meta_data  %>% 
    filter(!!sym(slot_var) == subsampled_cluster)  %>% 
    ungroup()  %>% 
    dplyr::summarise(total_cells = n(), 
             total_tn5 = sum(total))
    
    
    prcnt_10 <- subsampled_cluster_metrics$total_tn5 * .1
    up_range <- subsampled_cluster_metrics$total_tn5 + prcnt_10
    down_range <- subsampled_cluster_metrics$total_tn5 - prcnt_10

    passing_null_sample <- FALSE
  if (passing_null_sample == FALSE) {
      
      print("generating Null Sample")
      generated_subsample <- meta_data %>% 
        filter(!!sym(slot_var) != subsampled_cluster)  %>% 
        sample_n(nrow(.))  %>% 
        sample_n(subsampled_cluster_metrics$total_cells) %>% 
        mutate(sampled_slot = "ZZZZ_cell_pop") %>% 
        sample_n(nrow(.))
    
      generated_subsample_tn5 <- generated_subsample %>% 
        ungroup()  %>% 
        dplyr::summarise(total_tn5 = sum(total))
      
      if (between(generated_subsample_tn5, down_range, up_range) == TRUE){
        passing_null_sample <- TRUE
          } else {
        passing_null_sample <- FALSE
      }
      
    } else {
      print("Null Sample Generated")
  }
   
    
    
  
  sub_sample_group_1 <- generated_subsample %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, "QQQQQ", sep = "."))
  
  sub_sample_group_2 <- generated_subsample %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(sampled_slot, "CCCCC", sep = "."))
  
  combined_null <- bind_rows(sub_sample_group_1, sub_sample_group_2)
  
  subsampled_sparse <- filter(sparse_matrix, barcode %in% combined_null$cellID)
  subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
  
  subsample_combined_counts <- subsample_joined %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-sampled_slot) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  
  subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)

  return(subsample_accessability_counts_wide)
  
}

sample_cluster <- function(meta_data, sparse_matrix, slot_name, cluster){
  
    slot_var <- c(slot_name)
    
  subsampled_meta_data <- meta_data %>% 
    dplyr::select(cellID, !!sym(slot_var)) %>% 
    dplyr::filter(!!sym(slot_var) == cluster) %>%
    sample_n(nrow(.))
  
  group_1 <- subsampled_meta_data %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), "QQQQQ", sep = "."))
  
  group_2 <- subsampled_meta_data %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(!!sym(slot_var), "CCCCC", sep = "."))
  
  
  combined_clustering <- bind_rows(group_1, group_2)
  
  test_set_sparse <- dplyr::filter(sparse_matrix, barcode %in% combined_clustering$cellID)
  testing_join <- left_join(test_set_sparse, combined_clustering, by = c("barcode" = "cellID"))
  
  colnames(testing_join)
  
  combined_counts <- testing_join %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-!!slot_name) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  accessability_counts_wide <- pivot_wider(combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0) %>% 
    rowwise %>% 
    mutate(row.sum = sum(c_across(where(is.numeric)))) %>% 
    filter(row.sum > 5) %>% 
    select(-row.sum)

  return(accessability_counts_wide)
}
combine_clusters_prepare <- function(sample_cluster_sparse, null_cluster_sparse){
  
  
  combined_sample_real <- left_join(sample_cluster_sparse, null_cluster_sparse, by = c("gene_name")) %>% 
    replace(is.na(.), 0)
  
  final_count_data <- as.data.frame(combined_sample_real)
  rownames(final_count_data) <- final_count_data[,1]
  final_count_data[,1] <- NULL
  
  return(final_count_data)
  
}
generate_sample_matrix <- function(combined_sparse){
  
  colnames(combined_sparse)
  generate_sample_matrix <-colnames(combined_sparse)
  one_replace <- str_replace(generate_sample_matrix, "\\.QQQQQ", "")
  two_replace <- str_replace(one_replace, "\\.CCCCC", "")
  generate_sample_matrix <- two_replace 
  
  
  test <- unlist(generate_sample_matrix)
  sample_df <- as.data.frame(test)
  colnames(sample_df) <- "sample_type"
    
  take_reference_factor <- "ZZZZ_cell_pop"
  alternative_factor <- as.character(unique(sample_df[sample_df$sample_type != "ZZZZ_cell_pop",]))
    
  sample_df$sample_type <- factor(sample_df$sample_type, levels = c(take_reference_factor, alternative_factor))
    
    print(sample_df)
  
  return(sample_df)
}
run_de_seq_2 <- function(counts_matrix, sample_matrix,output_base){
  
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = sample_matrix,
                                design = ~ sample_type)
  
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  result_vector <- resultsNames(dds)
  grab_comparison_vector <- as.character(result_vector[2])
  
  
  output_name <-str_replace(grab_comparison_vector, "sample_type", "")

  res <- results(dds, name=grab_comparison_vector)
  resOrdered <- res[order(res$padj),]
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(resOrdered, rownames = "gene_name")
  
  final_output_name <- paste0(output_base,grab_comparison_vector, ".tsv")
  readr::write_tsv(as.data.frame(resOrdered.final),file=final_output_name)
  
  return(resOrdered)

}
call_cluster_specific_genes <- function(DE_seq_2_output,cluster_name,output_base){
  
  result_vector <- resultsNames(DE_seq_2_output)
  grab_comparison_vector <- as.character(result_vector[2])
  
  #For Saving Tsv
  resOrdered.final <- as_tibble(DE_seq_2_output, rownames = "gene_name")
  down_regulated <- get_upregulated(resOrdered.final)
  final_output_name_2 <- paste0(output_base, ".", cluster_name, ".",grab_comparison_vector, ".deseq2_output.tsv")
  readr::write_tsv(as.data.frame(down_regulated),file=final_output_name_2)
  
  return(resOrdered.final)
  
  
}
generate_marker_list <-function(DE_seq_2_output, marker_gene_bed, total_gene_bed, cluster, output_base){
  
  
  
  cluster_specific_genes.DA.2 <- get_upregulated(DE_seq_2_output) 
  cluster_specific_genes.DA <- cluster_specific_genes.DA.2
  
  #marker_in_de_genes <- marker_gene_bed %>% 
  #  dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
  #  left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
  #  filter(is.na(padj) != TRUE & padj < .1 ) %>% 
  #  arrange(padj) %>% 
  #  mutate(name = str_c(name, "pval", round(padj,4), cluster, sep = "_")) %>% 
  #  select(chr:geneID,name,type)
  
  
  
  `%ni%` <- Negate(`%in%`)
  clust_name <- c(cluster)
  bed_in_de_genes <- total_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
    #dplyr::filter(geneID %ni% as.character(marker_in_de_genes$geneID)) %>% 
    left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
    filter(is.na(padj) != TRUE) %>% 
    filter(padj < .1) %>% 
    select(chr:strand,padj) %>% 
    select(-val) %>% 
    arrange(padj) %>% 
    mutate(count = row_number()) %>% 
    mutate(name = str_c(geneID, "pval", round(padj,4), clust_name, sep = "_")) %>% 
    mutate(type = clust_name) %>%
    dplyr::select(-count, -strand, -padj) %>% 
    select(chr:geneID,name,type)
  
  
  combine_markers_new_denovo <- bind_rows(bed_in_de_genes)
  
  output_file = paste0(output_base,".",cluster, ".markers_de_novo.visualize.bed")
  readr::write_tsv(as.data.frame(combine_markers_new_denovo),file=output_file, col_names = TRUE)
  
  
}

grab_up_regulated_genes <- function(deseq2_output, cluster_name){
    
    result_vector <- resultsNames(deseq2_output)
    grab_comparison_vector <- as.character(result_vector[2])
  
    #For Saving Tsv
    resOrdered.final <- as_tibble(deseq2_output, rownames = "gene_name")
    cluster_specific <- get_upregulated(resOrdered.final) %>% 
        dplyr::mutate(origin = cluster_name)
    
    generate_output_name <- paste0("Cluster_", cluster_name, ".DA_list.txt")
    
    write_delim(cluster_specific, generate_output_name, 
        col_names = TRUE, quote = "none", delim = "\t")
    
    return(cluster_specific)
}


wrapper_function_updated <- function(meta_data, sparse_matrix, column_name,
                                    cluster_name,marker_bed, all_bed, 
                                     acr_file, output_base, output_location) {
  
  message(paste("Working on cluster: ",cluster_name))
  sample_cluster_test <- sample_cluster(meta_data, sparse_matrix, column_name, cluster_name)
  function_null <- generate_null_sample(meta_data, sparse_matrix, column_name,cluster_name)
  final_test <- combine_clusters_prepare(sample_cluster_test, function_null)
  deseq2_sample_matrix <- generate_sample_matrix(final_test)
  caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, output_base)
    message("Finding Upregulated genes") 
  up_reg_genes <- grab_up_regulated_genes(caught_output,cluster_name)
  #DA_marker_intersect <- intersect_ACRs_with_markers(up_reg_genes, acr_file, marker_bed)
  cluster_specific_genes <- call_cluster_specific_genes(caught_output, cluster_name, output_base)
  generate_marker_list(cluster_specific_genes, marker_bed, all_bed, cluster_name, output_base)
  return(up_reg_genes)
  
}


# Run on All --------------------------------------------------------------
print(colnames(loaded_meta_data))
final_cluster_list <- (unique(loaded_meta_data$V3_annotation_n))
message("Finding DA genes for the Following Clusters:")
print(final_cluster_list)
test <- lapply(final_cluster_list, wrapper_function_updated, meta_data = loaded_meta_data, sparse_matrix=loaded_sparse_matric,  marker_bed = marker_genes, 
               column_name = "V3_annotation_n", all_bed = all_genes, output_base = output_base, output_location = output_location)



