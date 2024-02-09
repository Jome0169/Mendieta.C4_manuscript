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


#working_dir <- "/scratch/jpm73279/comparative_single_cell/08.annotation_figures/zea_mays"

# load arguments
args <- commandArgs(T)
#if(length(args)!=5){stop("Rscript normGBA.R <gene.sparse> <meta> <Zea_mays.AGPv4.36.Allgene.nuclear.bed> <prefix> <F>")}
input <- as.character(args[1])
meta <- as.character(args[2])
gene <- as.character(args[3])
annot_col <- as.character(args[4])
prefix <- as.character(args[5])


slot_var <- c(annot_col)

#gene_bodysorghum_bicolor_tis_leaf_nmf_step_2_knn_27.out.de_novo.rds
#input <- here(working_dir,"00.data/zea_mays.gene_body_acc_leaf_V3_clustering_strict.counts.tab.txt")
#input_2 <- here("/home/jpm73279/r_script_dev/lw_plotting","zea_mays.normalized_gene_acc_scores.leaf_svd_knn_100_strict.GBaccessibility.sparse")
#meta <- here(working_dir,"00.data/Zm_leaf.V3_final.txt")
#gene <- here(working_dir,"00.data/Zm.markers.leaf.txt")
#gene_DA <- here(working_dir,"00.data/Zm-B73-REFERENCE-NAM_Zm00001eb.1.genes.bed")
#prefix <- "TEST_SORGHUM_TEST"


message("Loading Data....")
meta_data <- read.delim(meta)
gene_markers <- read.delim(gene)
gene_markers <- gene_markers  %>%
    arrange(type)

all_markers <- gene_markers$geneID

raw_cpm_counts_all_genes <- read_delim(input, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
    dplyr::mutate(cellID = barcode)  %>%
    dplyr::mutate(geneID = gene_name)

colnames(meta_data)


message("Calculating CPM Values....")
merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))  %>%
    group_by(!!sym(slot_var), geneID)  %>%
    summarise(counts = sum(accessability, na.rm = TRUE))

### Alt CPM Calc
merged_meta_cpm_information_copied <- merged_meta_cpm_information
catch <- merged_meta_cpm_information_copied  %>%
    group_by(!!sym(slot_var)) %>%
    group_map(~(cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>%
    unlist()

caught_values <- as_tibble(catch)
see <- ungroup(merged_meta_cpm_information_copied)
merged_meta_cpm_information_copied<- bind_cols(merged_meta_cpm_information_copied,caught_values)  %>% 
    rename(grouped_CPM = value)


message("Calculating ZScore Approximations...")
head(merged_meta_cpm_information_copied)
altered_deseq2 <- merged_meta_cpm_information_copied %>% 
    dplyr::select(-counts) %>% 
    pivot_wider(names_from = geneID, values_from = grouped_CPM, values_fill = 0) %>% 
    pivot_longer(cols = -!!sym(slot_var), names_to = "geneID", values_to = "grouped_CPM") %>% 
    group_by(geneID) %>% 
    mutate(Zscore = scale(grouped_CPM)) %>% 
    ungroup()  %>% 
    #mutate(relative_accessability = rescale(Zscore, to = c(0,1))) %>% 
    group_by(!!sym(slot_var))  %>% 
    mutate(Zscore_group = scale(Zscore))

message("Calculating Proportion of cells marker is Accessible IN...")
# Create Proportion Cells Accessible Metrics ------------------------------
merged_meta_cellID_values <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))
take_unq_genes <- unique(merged_meta_cellID_values$geneID)


merged_meta_cellID_values_all_genes <- merged_meta_cellID_values %>% 
    select(cellID, !!sym(slot_var), accessability, geneID) 


wider_all_genes_altered <- merged_meta_cellID_values_all_genes %>% 
                    pivot_wider(names_from = geneID, 
                                values_from = accessability,  
                                values_fill = 0) %>% 
                    pivot_longer(cols = c(-!!sym(slot_var), -cellID), 
                                 names_to = "geneID", 
                                 values_to = "accessability") %>% 
                    mutate(expression_bool = case_when(accessability < 1 ~ 0,
                                                       accessability >= 1 ~ 1)) %>% 
                    group_by(!!sym(slot_var), geneID) %>% 
                    summarise(total_cells = n(), 
                              proportion_expressing = (sum(expression_bool)/total_cells * 100))


message("Heirarchcial Clustering")
generate_vec <- c((slot_var), "geneID")
print(generate_vec)

marker_final_plotting <- left_join(altered_deseq2, wider_all_genes_altered, by = generate_vec) %>% 
        left_join(., gene_markers, by = c("geneID"))


message("Generating Final Matrix")
test_marker_clust <- marker_final_plotting  %>% 
    filter(geneID %in% gene_markers$geneID)  %>% 
    ungroup() %>% 
    dplyr::select(!!sym(slot_var), Zscore, name)  %>% 
    pivot_wider(names_from = !!sym(slot_var), values_from = Zscore) %>% 
      data.frame() # make df as tibbles -> matrix annoying

row.names(test_marker_clust) <- test_marker_clust$name  # put gene in `row`
test_marker_clust <- test_marker_clust[,-1] #drop gene column as now in rows
clust <- hclust(dist(test_marker_clust %>% as.matrix())) # hclust with distance matrix

message("Ordering Results")
mat <- marker_final_plotting %>% 
    filter(geneID %in% gene_markers$geneID)  %>% 
    ungroup() %>% 
    dplyr::select(!!sym(slot_var), Zscore, name)  %>% 
    pivot_wider(names_from = !!sym(slot_var), values_from = Zscore) %>% 
      data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$name  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix


message("Generating Dotplot...")
message("I hope it's pretty :) !")
options(repr.plot.width=12, repr.plot.height=5)
supplamental_marker_dotplot_final <- marker_final_plotting  %>% 
    mutate(Gene_name = factor(name, clust$labels[clust$order]),
          Cluster_name = factor(!!sym(slot_var), levels = v_clust$labels[v_clust$order])) %>% 
    filter(is.na(Gene_name) != TRUE) %>% 
    ggplot(., aes(y=Cluster_name, x = Gene_name, 
                               color = Zscore, size = proportion_expressing)) + 
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_colour_gradient2() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line.x = element_line(color="black", size = 1), 
          axis.line.y = element_line(color="black", size = 1)) +
    ggtitle("Z Score of Markers - Across Annotations")



generate_output_name <- paste0(prefix, ".dotmarkers.pdf")
ggsave(generate_output_name, plot = supplamental_marker_dotplot_final,
    width = 10, height = 5,
    units = c('in'), limitsize = FALSE, dpi = 300)



