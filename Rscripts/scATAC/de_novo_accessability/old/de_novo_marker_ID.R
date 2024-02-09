library(dplyr)
library(here)
library(DESeq2)
library(tidyverse)

# load arguments
args <- commandArgs(trailingOnly=T)

#args    
meta_data_file <- as.character(args[1])
gene_accessability_file <- as.character(args[2])
marker_gene_file <- as.character(args[3])
all_genes_bed <- as.character(args[4])
species <- as.character(args[5])
output_base <- as.character(args[6])

#Local Development Files
#data_path <- "/Users/feilab/Projects/05.comparative_single_cell/00.data/test_data_log2fc"
#meta_data_file <- here(data_path, "sb_leaf_nmf_compressed_markers.meta.txt")
#gene_accessability_file <- here(data_path, "sorghum_bicolor.normalized_gene_acc_scores.leaf_nmf.sctGBAcounts.sparse")
#1marker_gene_file <- here(data_path, "Sb_leaf.maize_markers.ortho.visualize.bed")


#Load Marker Genes
marker_genes <- read.table(marker_gene_file, header =TRUE)
#Load Meta Data
loaded_meta_data <- read.table(meta_data_file) %>% 
  mutate(Louvain_cluster_name = str_c("LouvainC", LouvainClusters, sep = "_"))


if (species == "sorghum_bicolor") {
  header_names <- c("chr", "star", "end", "geneID", "val", "strand")
  #all_gene_File <- here(data_path, "Sbicolorv5.1.gene.bed")
  all_genes <- read.table(all_genes_bed, header =FALSE) %>% 
      dplyr::select(V1:V6) %>% 
      setNames(header_names) %>%  
      mutate(geneID = str_remove_all(geneID,".g")) 
  #Load Sparse Matrix
  loaded_sparse_matric <- read_delim(gene_accessability_file, delim='\t', col_names = c("gene_name", "barcode", "accessability")) %>% 
      mutate(gene_name = str_remove_all(gene_name,".g"))
  
} else {
    
  header_names <- c("chr", "star", "end", "geneID", "val", "strand")
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

# Older version of null sample which didn't take into consideration the numbe
# of cells
generate_null_sample <- function(meta_data, sparse_matrix) {
  
  generated_subsample <- meta_data %>% 
    group_by(Louvain_cluster_name) %>% 
    slice_sample(prop=.25) %>% 
    ungroup() %>% 
    mutate(Louvain_cluster_name = "mixed_cell_pop") %>% 
    sample_n(nrow(.))
  
  sub_sample_group_1 <- generated_subsample %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "1", sep = "."))
  
  sub_sample_group_2 <- generated_subsample %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "2", sep = "."))
  
  combined_null <- bind_rows(sub_sample_group_1, sub_sample_group_2)
  
  subsampled_sparse <- filter(sparse_matrix, barcode %in% combined_null$cellID)
  subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
  
  subsample_combined_counts <- subsample_joined %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-Louvain_cluster_name) %>% 
    group_by(gene_name, louvain_grouping_sample) %>% 
    summarise(total_accessability = sum(accessability))
  
  
  
  subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)

  return(subsample_accessability_counts_wide)
  
}

#generate_null_sample <- function(meta_data, sparse_matrix, column_name, n_value) {
#  
#  var <- c(column_name)
#  print(var)
#  print("SubSampling for Null")
#  generated_subsample <- meta_data %>% 
#    sample_n(n_value) %>% 
#    #slice_sample(prop=.50) %>% 
#    ungroup() %>%
#    mutate(Louvain_cluster_name = "mixed_cell_pop") %>% 
#    sample_n(nrow(.))
#  
#  print("Splitting Null into Two groups")
#  sub_sample_group_1 <- generated_subsample %>% 
#    slice_head(prop = .5) %>% 
#    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "1", sep = "."))
#  
#  sub_sample_group_2 <- generated_subsample %>% 
#    slice_tail(prop = .5) %>% 
#    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "2", sep = "."))
#  
#  combined_null <- bind_rows(sub_sample_group_1, sub_sample_group_2)
#  
#  subsampled_sparse <- filter(sparse_matrix, barcode %in% combined_null$cellID)
#  subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
#  
#  subsample_combined_counts <- subsample_joined %>% 
#    dplyr::select(-barcode) %>%  
#    ungroup() %>% 
#    dplyr::select(-(!!column_name)) %>% 
#    group_by(gene_name, louvain_grouping_sample) %>% 
#    summarise(total_accessability = sum(accessability))
#  
#  
#  
#  subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, 
#                                                     names_from = "louvain_grouping_sample", 
#                                                     values_from = "total_accessability",
#                                                     values_fill = 0)
#
#  return(subsample_accessability_counts_wide)
#  
#}


sample_cluster <- function(meta_data, sparse_matrix, cluster){
  
  subsampled_meta_data <- meta_data %>% 
    dplyr::select(cellID, Louvain_cluster_name) %>% 
    dplyr::filter(Louvain_cluster_name == cluster) %>%
    sample_n(nrow(.))
  
  group_1 <- subsampled_meta_data %>% 
    slice_head(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "1", sep = "."))
  
  group_2 <- subsampled_meta_data %>% 
    slice_tail(prop = .5) %>% 
    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "2", sep = "."))
  
  
  combined_clustering <- bind_rows(group_1, group_2)
  
  test_set_sparse <- filter(sparse_matrix, barcode %in% combined_clustering$cellID)
  testing_join <- left_join(test_set_sparse, combined_clustering, by = c("barcode" = "cellID"))
  
  colnames(testing_join)
  
  combined_counts <- testing_join %>% 
    dplyr::select(-barcode) %>%  
    ungroup() %>% 
    dplyr::select(-Louvain_cluster_name) %>% 
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
  one_replace <- str_replace(generate_sample_matrix, "\\.1", "")
  two_replace <- str_replace(one_replace, "\\.2", "")
  generate_sample_matrix <- two_replace 
  
  
  test <- unlist(generate_sample_matrix)
  sample_df <- as.data.frame(test)
  colnames(sample_df) <- "sample_type"
  #rownames(sample_df) <- generate_sample_matrix
  sample_df$sample_type <- factor(sample_df$sample_type)
  
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
  down_regulated <- get_downregulated(resOrdered.final)
  final_output_name_2 <- paste0(output_base, ".", cluster_name, ".",grab_comparison_vector, ".cluster_specific.tsv")
  readr::write_tsv(as.data.frame(down_regulated),file=final_output_name_2)
  
  return(resOrdered.final)
  
  
}
generate_marker_list <-function(DE_seq_2_output, marker_gene_bed, total_gene_bed, cluster, output_base){
  
  
  
  cluster_specific_genes.DA.2 <- get_downregulated(DE_seq_2_output) 
  cluster_specific_genes.DA <- cluster_specific_genes.DA.2
  
  marker_in_de_genes <- marker_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
    left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
    filter(is.na(padj) != TRUE & padj < .1 ) %>% 
    arrange(padj) %>% 
    mutate(name = str_c(name, "pval", round(padj,4), cluster, sep = "_")) %>% 
    select(chr:geneID,name,type)
  
  
  
  `%ni%` <- Negate(`%in%`)
  clust_name <- c(cluster)
  bed_in_de_genes <- total_gene_bed %>% 
    dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
    dplyr::filter(geneID %ni% as.character(marker_in_de_genes$geneID)) %>% 
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
  
  
  combine_markers_new_denovo <- bind_rows(marker_in_de_genes, bed_in_de_genes)
  
  output_file = paste0(output_base,".",cluster, ".markers_de_novo.visualize.bed")
  readr::write_tsv(as.data.frame(combine_markers_new_denovo),file=output_file, col_names = TRUE)
  
  
}
wrapper_function <- function(meta_data, sparse_matrix, marker_bed, all_bed, cluster_name, output_base) {
  
  sample_cluster_test <- sample_cluster(meta_data, sparse_matrix, cluster_name)
  function_null <- generate_null_sample(meta_data, sparse_matrix)
  final_test <- combine_clusters_prepare(sample_cluster_test, function_null)
  deseq2_sample_matrix <- generate_sample_matrix(final_test)
  caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, output_base)
  cluster_specific_genes <- call_cluster_specific_genes(caught_output, cluster_name, output_base)
  generate_marker_list(cluster_specific_genes, marker_bed, all_bed, cluster_name, output_base)
  return(cluster_specific_genes)
  
}


# Run on All --------------------------------------------------------------
final_cluster_list <- (unique(loaded_meta_data$Louvain_cluster_name))
test <- lapply(final_cluster_list, wrapper_function, meta_data = loaded_meta_data, sparse_matrix=loaded_sparse_matric,  marker_bed = marker_genes, all_bed = all_genes, output_base = output_base)


# TESTING FUNCTIONS -------------------------------------------------------

#function_null <- generate_null_sample(loaded_meta_data, loaded_sparse_matric)
#sample_cluster_test <- sample_cluster(loaded_meta_data, loaded_sparse_matric, "LouvainC_9")
#final_test <- combine_clusters_prepare(sample_cluster_test, function_null)
#deseq2_sample_matrix <- generate_sample_matrix(final_test)
#caught_output <- run_de_seq_2(final_test, deseq2_sample_matrix, "TEST.LouvainC_0.")
#cluster_specific_genes <- call_cluster_specific_genes(caught_output, "TEST.LouvainC_0.")
#de.novo.markers <- generate_marker_list(cluster_specific_genes, marker_genes, all_genes, "LouvainC_9", "TEST.LouvainC_9")

#summary(caught_output)
#Function to pull out the genes which are DE in the general cluster
#final_de_set <- get_downregulated(caught_output)
#write.csv(as.data.frame(resOrdered),file=final_output_name)


# Graphing  -----------------------------------------------------------


#final_test %>% 
#  pivot_longer(., c("LouvainC_9.y","LouvainC_9.x","mixed_cell_pop.1","mixed_cell_pop.2"), names_to = "type") %>% 
#  ggplot(., aes(x = value,)) + geom_histogram(binwidth = 2) + facet_grid(type~.)
#
#final_test %>% 
#  pivot_longer(., c("LouvainC_9.y","LouvainC_9.x","mixed_cell_pop.1","mixed_cell_pop.2"), names_to = "type") %>% 
#  ggplot(., aes(x = value, y = type)) + geom_violin()
#
#
#final_test %>% 
#  pivot_longer(., c("LouvainC_9.y","LouvainC_9.x","mixed_cell_pop.1","mixed_cell_pop.2"), names_to = "type") %>% 
#  filter(type %in% c("LouvainC_9.x", "LouvainC_9.y")) %>% 
#  group_by(type) %>% 
#  summarise(x = quantile(value, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75))
#  
#
#cluster_specific_genes.DA.2 <- get_downregulated(cluster_specific_genes) 
#cluster_specific_genes.DA <- cluster_specific_genes.DA.2 %>% 
#  mutate(gene_name = str_remove_all(gene_name,".g"))
#
#colnames(marker_genes)
#
#marker_in_de_genes <- marker_genes %>% 
#  dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
#  left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
#  mutate(name = str_c(name, "pval", round(padj,4), sep = ".")) %>% 
#  select(chr:geneID,name,type)
#
#
#
#`%ni%` <- Negate(`%in%`)
#clust_name <- c("LouvainC_9")
#bed_in_de_genes <- all_genes %>% 
#  mutate(gene_name = str_remove_all(geneID,".g")) %>% 
#  dplyr::filter(geneID %in% as.character(cluster_specific_genes.DA$gene_name)) %>% 
#  dplyr::filter(geneID %ni% as.character(marker_in_de_genes$geneID)) %>% 
#  left_join(.,  cluster_specific_genes.DA, by = c("geneID" = "gene_name")) %>% 
#  filter(is.na(padj) != TRUE) %>% 
#  filter(padj < .05) %>% 
#  select(chr:strand,padj) %>% 
#  select(-val) %>% 
#  arrange(padj) %>% 
#  mutate(count = row_number()) %>% 
#  mutate(name = str_c("DA.gene", count, "pval", round(padj,4), sep = ".")) %>% 
#  mutate(type = clust_name) %>%
#  dplyr::select(-count, -strand, -padj) %>% 
#  select(chr:geneID,name,type)
#  
#
#View((bed_in_de_genes))
#
#
#combine_markers_new_denovo <- bind_rows(marker_in_de_genes, bed_in_de_genes)
#View(combine_markers_new_denovo)
#
#readr::write_tsv(as.data.frame(combine_markers_new_denovo),file="first_output_test.tsv", col_names = TRUE)




# Develop Generate Subsample ------------------------------------------------------
#generated_subsample <- loaded_meta_data %>% 
#    group_by(Louvain_cluster_name) %>% 
#    slice_sample(n=100) %>% 
#    ungroup() %>% 
#    mutate(Louvain_cluster_name = "mixed_cell_pop") %>% 
#    sample_n(nrow(.))
#
#sub_sample_group_1 <- generated_subsample %>% 
#  slice_head(prop = .5) %>% 
#  mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "1", sep = "."))
#
#sub_sample_group_2 <- generated_subsample %>% 
#  slice_tail(prop = .5) %>% 
#  mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "2", sep = "."))
#
#combined_null <- bind_rows(sub_sample_group_1, sub_sample_group_2)
#
#subsampled_sparse <- filter(loaded_sparse_matric, barcode %in% combined_null$cellID)
#subsample_joined <- left_join(subsampled_sparse, combined_null, by = c("barcode" = "cellID"))
#
#
#subsample_combined_counts <- subsample_joined %>% 
#  dplyr::select(-barcode) %>%  
#  ungroup() %>% 
#  dplyr::select(-Louvain_cluster_name) %>% 
#  group_by(gene_name, louvain_grouping_sample) %>% 
#  summarise(total_accessability = sum(accessability))
#
#
#subsample_accessability_counts_wide <- pivot_wider(subsample_combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)


# Develop Generate Test Set -------------------------------------------------------
#loaded_meta_data %>% 
#    group_by(Louvain_cluster_name) %>% 
#    summarise(n()) %>% 
#    View()
#
#subsampled_meta_data <- loaded_meta_data %>% 
#  dplyr::select(cellID, Louvain_cluster_name) %>% 
#  dplyr::filter(Louvain_cluster_name == "LouvainC.6") %>%
#  sample_n(nrow(.))
#  
#group_1 <- subsampled_meta_data %>% 
#    slice_head(prop = .5) %>% 
#    mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "1", sep = "."))
#
#group_2 <- subsampled_meta_data %>% 
#  slice_tail(prop = .5) %>% 
#  mutate(louvain_grouping_sample = str_c(Louvain_cluster_name, "2", sep = "."))
#
#
#combined_clustering <- bind_rows(group_1, group_2)
#
#test_set_sparse <- filter(loaded_sparse_matric, barcode %in% combined_clustering$cellID)
#testing_join <- left_join(test_set_sparse, combined_clustering, by = c("barcode" = "cellID"))
#
#colnames(testing_join)
#
#combined_counts <- testing_join %>% 
#    dplyr::select(-barcode) %>%  
#    ungroup() %>% 
#    dplyr::select(-Louvain_cluster_name) %>% 
#    group_by(gene_name, louvain_grouping_sample) %>% 
#    summarise(total_accessability = sum(accessability))
#  
#  
#accessability_counts_wide <- pivot_wider(combined_counts, names_from = "louvain_grouping_sample", values_from = "total_accessability",   values_fill = 0)



# Develop Merge Data - Run Deseq2 -------------------------------------------------

#combined_sample_real <- full_join(subsample_accessability_counts_wide, accessability_counts_wide, by = c("gene_name")) %>% 
#  replace(is.na(.), 0)
#
#final_count_data <- as.data.frame(combined_sample_real)
#rownames(final_count_data) <- final_count_data[,1]
#final_count_data[,1] <- NULL
#
#colnames(combined_sample_real)
#generate_sample_matrix <-colnames(combined_sample_real)[2:length(colnames(combined_sample_real))]
#one_replace <- str_replace(generate_sample_matrix, ".1", "")
#two_replace <- str_replace(one_replace, ".2", "")
#generate_sample_matrix <- two_replace 
#
#
#test <- unlist(generate_sample_matrix)
#sample_df <- as.data.frame(test)
#colnames(sample_df) <- "sample_test"
##rownames(sample_df) <- generate_sample_matrix
#sample_df$sample_test <- factor(sample_df$sample_test)



#  Develop Run DEseq2 --------------------------------------------------------------
#dds <- DESeqDataSetFromMatrix(countData = final_count_data,
#                              colData = sample_df,
#                              design = ~ sample_test)
##keep <- rowSums(counts(dds)) >= 1
##dds <- dds[keep,]
#
#dds <- DESeq(dds, useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)
#res <- View(results(dds))
#View(res)
#
#test <- resultsNames(dds)
#
#resLFC <- lfcShrink(dds, coef=test, type="apeglm")
#View(resLFC)
#res
#
#quick_make <- as.character(test[2])
#quick_make
#
#res <- results(dds, name=(quick_make))
#resOrdered <- res[order(res$padj),]
#summary(res)
#sum(res$padj < 0.1, na.rm=TRUE)
#plotMA(res, ylim=c(-2,2))
#
#write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")
#
#stored_results <- as.data.frame(resOrdered)
#
#resOrdered[resOrdered$padj < .05, na.rm=TRUE]
#
#
#
#caught_output$lfcSE
##summary(resOrdered, altHypothesis = greater)
#quick <- results(dds, alpha = 0.05, lfcThreshold = 0.58)
#quick
#final <- quick[quick$padj < .05,]
#head(final)

