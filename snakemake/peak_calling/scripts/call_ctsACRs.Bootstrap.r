#Please comment thiis out if you're not using my Dir and Notebook. This is hacky to make it run in snakemake :/ 
.libPaths(c("/home/jpm73279/R/x86_64-conda-linux-gnu-library/4.1","/home/jpm73279/.conda/envs/R_final_install/lib/R/library", .libPaths()))


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
library(rsample)
library(purrr)
library(data.table)
library(future)
library(furrr)
library(parallel)
library(argparse)

# load arguments
args <- commandArgs(T)

calculating_specificity <- function(x, threads=30){
    
    # add pseudo-count
    x <- x+1
    
    # convert to probability distribution
    p <- t(apply(x, 1, function(z){
        z/sum(z)
    }))
    hp <- apply(p, 1, function(z){
        z <- z[z > 0]
        -1*sum(z*(log2(z)))
        
    })
    sp <- apply(p, 2, function(z){
        hp - log2(z)
    })
    
    return(sp)
}

convert_to_sparse_matrix <- function(three_col_tribble, meta_slot_var) {
    
    CPM_matrix_prep <- three_col_tribble  %>% 
        dplyr::select(!!sym(meta_slot_var), geneID, grouped_CPM)

    three_col_prep <- CPM_matrix_prep  %>% 
        dplyr::rename("V2" = "geneID")  %>% 
        dplyr::rename("V1" = !!sym(meta_slot_var))


    # make sure bins/cells are factors
    three_col_prep$V1 <- factor(three_col_prep$V1)
    three_col_prep$V2 <- factor(three_col_prep$V2)


    # convert to sparseMatrix format
    sparse_count_matrix <- Matrix::sparseMatrix(i=as.numeric(three_col_prep$V1),
        j=as.numeric(three_col_prep$V2),
        x=as.numeric(three_col_prep$grouped_CPM),
        dimnames=list(levels(three_col_prep$V1),levels(three_col_prep$V2)))


    return(sparse_count_matrix)
    
}

# Function to generate null distribution
generate_null_distribution <- function(data, col_name) {
  # Determine the number of rows per class (equal for all classes)
  n_rows_per_class <- min(250)

  ## Replicate splitting is old and being abandoned in place of bootstrapping
  # Split the data into replicate1 and replicate2 groups based on the rep_values column
  #data_rep1 <- data %>% filter(data[[rep_values]] == "rep1")
  #data_rep2 <- data %>% filter(data[[rep_values]] == "rep2")

  # Group data by the specified column and select n_rows_per_class for each class 
  sampled_data <- data %>%
    group_by_at(col_name) %>%
    sample_n(n_rows_per_class,replace = TRUE) %>%
    ungroup()

  # Shuffle the column and make sure less than 20% of cells retain their original cluster ID
  #shuffled_column <- shuffle_and_assign_v3(sampled_data, col_name)
  #sampled_data[[col_name]] <- shuffled_column
   sampled_data[[col_name]] <- sampled_data[[col_name]][sample(nrow(sampled_data))]

  return(sampled_data)
}


generate_null_dist_values <- function(meta_data, meta_slot_var, raw_cpm_counts_all_genes){
#    meta_slot_var <- c("final_annotation")
    options(dplyr.summarise.inform = FALSE)
    merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"), relationship = "many-to-many")  %>%
        group_by(!!sym(meta_slot_var), geneID)  %>%
        summarise(counts = sum(accessability, na.rm = TRUE))

    ### Alt CPM Calc
    merged_meta_cpm_information_copied <- merged_meta_cpm_information
    catch <- merged_meta_cpm_information_copied  %>%
        group_by(!!sym(meta_slot_var)) %>%
        group_map(~(edgeR::cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>%
        unlist()

    caught_values <- as_tibble(catch)
    see <- ungroup(merged_meta_cpm_information_copied)
    merged_meta_cpm_information_copied<- bind_cols(merged_meta_cpm_information_copied,caught_values)  %>% 
        rename(grouped_CPM = value)
    
    sparse_null_dist <- convert_to_sparse_matrix(merged_meta_cpm_information_copied, meta_slot_var)
    transposed_ACRs_by_ct <- as.matrix(t(sparse_null_dist))

    message("Generating Null Distribution ...")

    calculate_specificity <- calculating_specificity(transposed_ACRs_by_ct)
       
    return(calculate_specificity)
}

generate_null_dist_values_optimized <- function(meta_data, meta_slot_var, raw_cpm_counts_all_genes) {

  #message(paste0("Working on bootstrap ", counter))
    print(paste("Processing permutation:", Sys.getpid()))

  # Convert dataframes to data.tables
  setDT(meta_data)
  setDT(raw_cpm_counts_all_genes)
  
  # Replace left_join with merge
  merged_meta_cpm_information <- merge(meta_data, raw_cpm_counts_all_genes, by = "cellID", all.x = TRUE, allow.cartesian=TRUE)
  
  # Group by and summarise
  merged_meta_cpm_information[, counts := sum(accessability, na.rm = TRUE), by = c(meta_slot_var, "geneID")]

  # Calculate CPM
  merged_meta_cpm_information[, grouped_CPM := edgeR::cpm(counts, log = FALSE, group = get(meta_slot_var)), by = meta_slot_var]
  
  sparse_null_dist <- convert_to_sparse_matrix(merged_meta_cpm_information, meta_slot_var)
  transposed_ACRs_by_ct <- as.matrix(t(sparse_null_dist))

  message("Generating Null Distribution ...")

  calculate_specificity <- calculating_specificity(transposed_ACRs_by_ct)
       
  return(calculate_specificity)
}


generate_pvalues_bootstraps_fast <- function(meta_data, raw_cpm_counts_all_genes, meta_slot_var) {
  
  message(paste0("Working on bootstrap..."))
  #counter <<- counter +1 

  setDT(meta_data)
  setDT(raw_cpm_counts_all_genes)
  
  merged_meta_cpm_information <- merge(meta_data, raw_cpm_counts_all_genes, by = "cellID", all.x = TRUE, allow.cartesian=TRUE)
  
  merged_meta_cpm_information[, counts := sum(accessability, na.rm = TRUE), by = c(meta_slot_var, "geneID")]
  
  #message("generating the CPM values")
  
  merged_meta_cpm_information[, grouped_CPM := edgeR::cpm(counts, log = FALSE, group = get(meta_slot_var)), by = meta_slot_var]
  
  sparse_null_dist <- convert_to_sparse_matrix(merged_meta_cpm_information, meta_slot_var)
  transposed_ACRs_by_ct <- as.matrix(t(sparse_null_dist))
  
  #message("Generating P values ...")
  
  calculate_specificity <- calculating_specificity(transposed_ACRs_by_ct)
  return(calculate_specificity)
}

parser <- ArgumentParser(description = "Command Line Args....")

# # Define the arguments
parser$add_argument("--input_data", required=TRUE, help = "The input data file.")
parser$add_argument("--meta", required=TRUE, help = "The meta data file.")
parser$add_argument("--peak_file", required=TRUE, help = "The peak file.")
parser$add_argument("--meta_slot", required=TRUE, help = "The meta slot.")
parser$add_argument("--prefix", required=TRUE, help = "The prefix for output files.")

# Add new arguments
parser$add_argument("--threshold", type="numeric", required=TRUE, help="The threshold value for statistical testing.")
parser$add_argument("--stat_test", default="perm", choices=c("perm", "pnorm"), help="The statistical test to use. Choose 'perm' for permutation test, or 'pnorm' for normal distribution test. Defaults to 'perm'.")
parser$add_argument("--null_permutations", type="integer", required=FALSE, default=5000, help="The number of null permutations to generate. Default is 1000.")
parser$add_argument("--entropy_bootstraps", type="integer", required=FALSE, default=1000, help="The number of bootstraps for the entropy metric to generate. Default is 1000.")


args <- parser$parse_args()

# Access the values with args$<name>
input_data <- args$input_data
meta <- args$meta
peak_file <- args$peak_file
meta_slot <- args$meta_slot
prefix <- args$prefix

# Access new arguments
threshold <- args$threshold
stat_test <- args$stat_test
null_permutations <- args$null_permutations
entropy_bootstraps <- args$entropy_bootstraps



# ## Testing Data Original 
# input_data <- "/scratch/jpm73279/comparative_single_cell/07.call.ACRs/replicate_analysis_one_off/zm/zm.peaks_accessability.txt"
# meta <- "/scratch/jpm73279/comparative_single_cell/07.call.ACRs/replicate_analysis_one_off/zm/Zm.leaf_annot.V5.meta.frozen.txt"
# peak_file <- "/scratch/jpm73279/comparative_single_cell/07.call.ACRs/replicate_analysis_one_off/zm/zm.peaks.500bp_peaks.bed"
# meta_slot <- "final_annotation_n"
# prefix <- "finalizing_cts_ACR_calling"
# threshold <- .01
# stat_test <- "perm"
# null_permutations <- 40
# entropy_bootstraps <- 10



## Read Inputs 
input <- input_data
message("Reading Peak File")
bed_file_read <- read_delim(peak_file, col_names = c("chrom", "start", "stop", "acr_number", "accessability"))

message("Reading Meta Data...")
meta_data <- read.delim(meta)

## Use the column for meta_data for cell type ACR calling 
meta_slot_var <- c(meta_slot)


message("Reading Input Data...")
raw_cpm_counts_all_genes <- read_delim(input, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
    dplyr::mutate(cellID = barcode)  %>%
    dplyr::mutate(geneID = gene_name)


## Generate the Null distribution of values... 
null_distributions <- replicate(null_permutations, generate_null_distribution(meta_data, meta_slot_var), simplify = FALSE)
#Older slower implementation 
#null_dist_values <- lapply(null_distributions, generate_null_dist_values_optimized, "final_annotation_n", raw_cpm_counts_all_genes)


counter <- 1
null_dist_values <- mclapply(null_distributions, generate_null_dist_values_optimized, 
         meta_slot_var, raw_cpm_counts_all_genes,
         mc.preschedule = FALSE, mc.set.seed = TRUE,
         mc.silent = FALSE, mc.cores = 25,
         mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL)

is_valid_matrix <- function(mat) {
  # Check if it's a matrix
  if(!is.matrix(mat)) return(FALSE)
  
  # Check if it's empty
  if(any(dim(mat) == 0)) return(FALSE)
  
  # Check for negative values
  if(any(mat < 0, na.rm = TRUE)) return(FALSE)
  
  # If passed all checks
  return(TRUE)
}

# Apply the function to the list of matrices
valid_indices <- sapply(null_dist_values, is_valid_matrix)
# Filter out invalid matrices
null_dist_values <- null_dist_values[valid_indices]










message("Generating 1000 Bootstraps")
y <- bootstraps(meta_data, times = entropy_bootstraps, strata = !!sym(meta_slot_var))

# Set up parallel processing
message("Running Bootstraps...")

## Older Slower Version
# counter <- 1
# tic()
# results <- y %>%
#   mutate(p_values = map(splits, ~ generate_pvalues_bootstraps_fast(analysis(.x), raw_cpm_counts_all_genes, "final_annotation_n"))) %>%
#   select(id, p_values)
# toc()

#Set to 3gb
options(future.globals.maxSize= 9e+9)
## Increase speed by threading...

generate_pvalues_bootstraps_fast_safe_function <- purrr::safely(generate_pvalues_bootstraps_fast)
counter <- 1
plan(multicore)
results <- y %>%
 mutate(safe_output = future_map(splits, ~ generate_pvalues_bootstraps_fast_safe_function(analysis(.x), raw_cpm_counts_all_genes, meta_slot_var))) %>%
 mutate(has_error = map_lgl(safe_output, ~ !is.null(.x$error))) %>%
 # Filter out rows with errors
 filter(!has_error) %>%
 # Extract p_values from the safe_output
 mutate(p_values = map(safe_output, "result")) %>%
 select(id, p_values)

expaneded_bootstraps <- results %>%
  mutate(p_values = map(p_values, ~ as_tibble(.x, rownames = "ACR_values"))) %>% # Convert matrices to tibbles and include row names
  unnest(cols = p_values) %>%                          # Unnest the tibbles
  pivot_longer(cols = -c(id, ACR_values),              # Keep 'id' and 'ACR_values' fixed
               names_to = "cell_type",                 # Assign the column names to 'cell_type'
               values_to = "value") %>%                # Assign the values to 'value'
  rename(BootstrapID = id)                             # Rename the columns to the desired names

nested_data <- expaneded_bootstraps %>% 
    select(-BootstrapID) %>%
    group_by(ACR_values, cell_type) %>%
    nest() %>% 
    rename(distribution = "data")


message("Merging Bootstraps and Null...")

melted_nested_nulls <- imap_dfr(null_dist_values, function(matrix, index) {
  df_matrix <- reshape2::melt(matrix)
  df_matrix <- df_matrix %>%
    rename(
      row_name = Var1,
      column_name = Var2,
      value = value
    ) %>%
    mutate(matrix_index = index)
  return(df_matrix)
})


rm(null_dist_values)
gc()

null_dist_generation <- melted_nested_nulls %>% 
    rename(ACR_values = row_name) %>%
    select(ACR_values, value) %>%
    group_by(ACR_values) %>%
    nest() %>%
    rename(null_dist = data)

merged_bootstraps_nulls <- left_join(nested_data, null_dist_generation, by = c("ACR_values"))


# ## Generate the same plot looking at ACRs associated with marker genes 
# options(repr.plot.width=15, repr.plot.height=15)
# #PEPC1 ME3 GL1 LRD3
# look_group <- c("scACR_51864", "scACR_18990", "scACR_43190", "scACR_53777")

# plot_acr_null_real <- merged_bootstraps_nulls %>% 
#     ungroup() %>% 
#     dplyr::filter(ACR_values %in% look_group) %>% 
#     mutate(ACR_values = factor(ACR_values, levels = look_group)) %>%  # Reorder the levels of ACR_values
#     unnest(distribution) %>% 
#     rename(real_value = value) %>% 
#     unnest(null_dist)%>% 
#     rename(null_value = value) %>% 
#     pivot_longer(c(real_value,null_value), names_to = "class", values_to = "val") 

# acr_meds <- plot_acr_null_real %>% 
#     group_by(class, cell_type, ACR_values) %>% 
#     summarise(median_val = mean(val))
 
# ggplot(plot_acr_null_real, aes(val, color = class)) + geom_density() + facet_grid(cell_type~ACR_values, scales="free_y") +
#       geom_vline(data=acr_meds, aes(xintercept=median_val, color=class),
#              linetype="dashed")

calc_pvals <- function(qp, mean_val, sd) {
  obs <- qp[is.finite(qp)]
  #ave <- mean(all_null_values_array, na.rm = TRUE)
  #sd <- sd(all_null_values_array, na.rm = TRUE)
  pvals <- pnorm(obs, mean = mean_val, sd = sd, lower.tail = TRUE)
  return(pvals)
}

# calculate_pvalues_per_ACR <- function(ACR_tribble_list) {
#     results <- ACR_tribble_list %>%
#         ungroup() %>%
#         rowwise() %>%
#         mutate(
#            list_len = lengths(null_dist) + 1,
           
#            median_val = mean(unlist(distribution), na.rm = TRUE),
#            lower_tail_test = (sum(unlist(distribution) > unlist(null_dist))),
#            perm_pval = lower_tail_test/list_len,
            
#            null_dist_mean = mean(unlist(null_dist), na.rm = TRUE),
#            null_dist_sd = sd(unlist(null_dist), na.rm = TRUE),
#            pnorm_pval = map(null_dist, ~calc_pvals(median_val, null_dist_mean, null_dist_sd)),
#            z_score = (median_val - null_dist_mean) / null_dist_sd) %>% # calculate Z score %>% 
#     mutate(pnorm_pval = map_dbl(pnorm_pval, ~ .x[[1]])) %>% 
#     ungroup()
#     return(results)
# }

calculate_pvalues_per_ACR <- function(ACR_tribble_list) {
    results <- ACR_tribble_list %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
           list_len = lengths(null_dist) + 1,
           median_val = mean(unlist(distribution), na.rm = TRUE),
            conf_interval_lower = quantile(unlist(distribution), probs = c(0.025), na.rm = TRUE),           
            conf_interval_upper = quantile(unlist(distribution), probs = c(0.975), na.rm = TRUE),           
            med_val_perm = sum(median_val > unlist(null_dist)),
            perm_pval = med_val_perm/list_len, # Permutation p-value for the median
           perm_pval_lower = sum(unlist(conf_interval_lower[1][1]) > unlist(null_dist))/list_len, # Permutation p-value for the lower CI bound
           perm_pval_upper = sum(unlist(conf_interval_upper[1][1]) > unlist(null_dist))/list_len, # Permutation p-value for the upper CI bound
           null_dist_mean = mean(unlist(null_dist), na.rm = TRUE),
           null_dist_sd = sd(unlist(null_dist), na.rm = TRUE),
           pnorm_pval = map(null_dist, ~calc_pvals(median_val, null_dist_mean, null_dist_sd)),
           z_score = (median_val - null_dist_mean) / null_dist_sd) %>% # calculate Z score 
           zscore_pval = pnorm(z_score, lower.tail = TRUE) %>% # calculate p-value based on Z score 
        mutate(pnorm_pval = map_dbl(pnorm_pval, ~ .x[[1]])) %>%
        ungroup()
    return(results)
} 


# look_group <- c("scACR_51864", "scACR_18990")#, "scACR_43190", "scACR_53777")

# test_ACR_group <- merged_bootstraps_nulls %>% 
#     #ungroup() %>% 
# #    dplyr::filter(ACR_values %in% look_group)
# #test_stats <- calculate_pvalues_per_ACR(test_ACR_group)

message("Calculating Pvalues....")
calculated_Pvals_all <- calculate_pvalues_per_ACR(merged_bootstraps_nulls)
## Save the RDS object here after all the hard work is done... 

message("Saving Intermediate Data")
save_pval_null_dists <- paste0(prefix, ".combined_data.rds")
saveRDS(calculated_Pvals_all, file = save_pval_null_dists)

message("Filtering Pvalues...")
calculated_Pvals_all_quantify <- calculated_Pvals_all %>% 
    dplyr::select(ACR_values, cell_type, perm_pval, perm_pval_lower, perm_pval_upper, pnorm_pval, z_score) 

save_pval_null_dists <- paste0(prefix, ".all_pvalues.csv")
write_delim(calculated_Pvals_all_quantify, file = save_pval_null_dists, delim=",")


## Clean up happen later after script is finalized. Stops from StackOverflow
rm(calculated_Pvals_all)
rm(merged_bootstraps_nulls)
gc()

#setwd("/scratch/jpm73279/comparative_single_cell/dev_location/entropy_final")

quantify_cell_type_specific_acrs <- function(signifigant_matrix, pval_slot, pval_filter, prefix) {
        
  # Save histogram of p-values
    hist_save <- paste0(prefix, ".pvalue_dist.png")
    png(hist_save)
    hist(as.numeric(as.matrix(signifigant_matrix[,pval_slot])), breaks = 100)
    dev.off()
    
    
    signifigant_matrix <- signifigant_matrix %>% 
        rename(pval = !!sym(pval_slot))
    
    print(head(signifigant_matrix))
    ct_specific_acrs_filtered_cell_types <- signifigant_matrix %>% 
        dplyr::filter(z_score < -1) %>%
        dplyr::filter(pval < pval_filter)
    print(head(ct_specific_acrs_filtered_cell_types))

    
    counts_ct <- ct_specific_acrs_filtered_cell_types  %>% 
        group_by(ACR_values) %>% 
        summarise(count_n_cell_types = n())  %>% 
        arrange(desc(count_n_cell_types))  %>% 
        mutate(class_acr = case_when(count_n_cell_types == 1 ~ "cts_acr",
                                count_n_cell_types > 1 & count_n_cell_types <=3 ~ "ctr_acr", 
                                count_n_cell_types > 3 ~ "broadly_accessible_acr"))

    
    
      # Check 1
      if (any(duplicated(counts_ct$ACR_values))) {
        stop("Error: Some ACRs are assigned to multiple categories.")
      }


    ## Different filtering to assign and filter different classes of ACRs
    cell_type_specific_restricted_ACRs <- counts_ct  %>% 
        dplyr::filter(class_acr == "cts_acr" | class_acr == "ctr_acr" )

    cell_type_specific_ACRs <- counts_ct  %>% 
        dplyr::filter(class_acr == "cts_acr")

    cell_type_restricted_ACRs <- counts_ct  %>% 
        dplyr::filter(class_acr == "ctr_acr" )
    

    ct_specific_acrs_filtered_cell_types.filtered <- ct_specific_acrs_filtered_cell_types  %>% 
        dplyr::filter(ACR_values %in% cell_type_specific_restricted_ACRs$ACR_values)
    

    write_delim(ct_specific_acrs_filtered_cell_types.filtered, file = paste0(prefix, ".pvalues.csv", collapse = "."), delim=",")
    

    ##########################################
    ###Isoalte write cell-type specific ACRs
    ##########################################
    ct_specific_acrs_filtered_cell_types_joined_names <- ct_specific_acrs_filtered_cell_types  %>% 
        ungroup() %>% 
        dplyr::filter(ACR_values %in% cell_type_specific_ACRs$ACR_values)  %>% 
        mutate(sc_acr_name = str_c(ACR_values, cell_type, sep = ";"))
    
    
     # Check 2
      if (nrow(unique(ct_specific_acrs_filtered_cell_types_joined_names["ACR_values"])) != nrow(cell_type_specific_ACRs)) {
        stop("Error: Mismatch in number of ACRs assigned to cell-type specific category.")
      }
    

    combined <- left_join(ct_specific_acrs_filtered_cell_types_joined_names, bed_file_read, by = c("ACR_values" = "acr_number"))

    # Write Cell Type Specific Peaks ALl
    cell_type_specific_acrs <- combined  %>% 
        dplyr::ungroup() %>% 
        #mutate_at(vars(!!sym(pval_slot)), character)  %>% 
        dplyr::select(chrom, start, stop, sc_acr_name, pval) 

    write_tsv(cell_type_specific_acrs, file = paste0(prefix, ".all_cts.ACRs.bed", collapse = "."), col_names = FALSE)
    
    # Write Cell Type Specific Peaks by Cell Type 
    combined  %>% 
        group_by(cell_type)  %>% 
        dplyr::select(chrom, start, stop, sc_acr_name, pval)  %>% 
        group_walk(~ write_delim(.x, paste0(prefix,".", .y$cell_type, ".cts.ACRs.bed", collapse = "."), delim = "\t", col_names = FALSE))

    
    ##########################################
    #Select cell-type restricted ACRs and Write
    ##########################################
    ctr_acrs_filtered_cell_types_joined_names <- ct_specific_acrs_filtered_cell_types  %>% 
        dplyr::filter(ACR_values %in% cell_type_restricted_ACRs$ACR_values) %>%
        group_by(ACR_values) %>% 
        summarize(combined_cell_type = paste(cell_type, collapse = ','),
              pval = paste(pval, collapse = ","))  %>% 
        ungroup()  %>% 
        mutate(sc_acr_name = str_c(ACR_values, combined_cell_type, sep = ";"))

    combined_ctr_acrs <- left_join(ctr_acrs_filtered_cell_types_joined_names, bed_file_read, by = c("ACR_values" = "acr_number"))

    # Write Cell Type restricted Peaks ALl
    cell_type_restricted_acrs <- combined_ctr_acrs  %>% 
        dplyr::select(chrom, start, stop, sc_acr_name, pval)

    write_tsv(cell_type_restricted_acrs, file = paste0(prefix, ".all_ctr.ACRs.bed", collapse = "."), col_names = FALSE)

 
    ##########################################   
    #Select Broadly Accessible ACRs and write 
    ##########################################
    '%ni%' <- Negate("%in%")
    broadly_accessible_acrs <- bed_file_read  %>% 
        dplyr::filter(acr_number %ni% cell_type_specific_restricted_ACRs$ACR_values)  %>% 
        dplyr::mutate(sc_acr_name = str_c(acr_number, "broadly_accessible", sep = ";"))  %>% 
        dplyr::select(-acr_number) %>%
        dplyr::select(-accessability) %>%
        mutate(pval = NA) %>% 
        dplyr::select(chrom, start, stop, sc_acr_name, pval)

    #write_tsv(broadly_accessible_acrs, file = paste0(prefix, "broadly_accessible.ACRs.bed", collapse = "."), col_names = FALSE)
    
    #Combine All ACRs and Write
    
    cell_type_specific_acrs <- cell_type_specific_acrs %>% 
            mutate(pval = as.character(pval))
    
    cell_type_restricted_acrs <- cell_type_restricted_acrs %>% 
        mutate(pval = as.character(pval))
    
    all_acrs_combined <- bind_rows(broadly_accessible_acrs, cell_type_specific_acrs, cell_type_restricted_acrs)
    write_tsv(all_acrs_combined, file = paste0(prefix, ".all_ACRs.classified.bed", collapse = "."), col_names = FALSE)
    
    return(all_acrs_combined)
}



if (stat_test == "perm") {
  quantify_cell_type_specific_acrs(calculated_Pvals_all_quantify, "perm_pval", threshold, prefix)
} else if (stat_test == "pnorm") {
  quantify_cell_type_specific_acrs(calculated_Pvals_all_quantify, "pnorm_pval", threshold, prefix)
} else {
  stop("Invalid stat_test value")
}

message("Done! Check your outputs Dork! :D ")


