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


working_dir <- "/scratch/jpm73279/comparative_single_cell/dev_location/entropy_calc.CTs/R_implemenation_test"

# load arguments
args <- commandArgs(T)
input_data <- as.character(args[1])
meta <- as.character(args[2])
peak_file <- as.character(args[3])
meta_slot <- as.character(args[4])
prefix <- as.character(args[5])

## Read Inputs 
input <- input_data
bed_file_read <- read_delim(peak_file, col_names = c("chrom", "start", "stop", "acr_number", "accessability"))
meta_data <- read.delim(meta)

## Use the column for meta_data for cell type ACR calling 
meta_slot_var <- c(meta_slot)
##!!sym(meta_slot_var)


raw_cpm_counts_all_genes <- read_delim(input, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
    dplyr::mutate(cellID = barcode)  %>%
    dplyr::mutate(geneID = gene_name)


merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"))  %>%
    #mutate(safe_cluster_name = str_c("Louvain_C", LouvainClusters, sep ="_"))  %>%
    #dplyr::select(-LouvainClusters)  %>%
    group_by(!!sym(meta_slot_var), geneID)  %>%
    summarise(counts = sum(accessability, na.rm = TRUE))

### Alt CPM Calc
merged_meta_cpm_information_copied <- merged_meta_cpm_information
catch <- merged_meta_cpm_information_copied  %>%
    group_by(!!sym(meta_slot_var)) %>%
    group_map(~(cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>%
    unlist()



caught_values <- as_tibble(catch)
see <- ungroup(merged_meta_cpm_information_copied)
merged_meta_cpm_information_copied<- bind_cols(merged_meta_cpm_information_copied,caught_values)  %>% 
    rename(grouped_CPM = value)

specificity <- function(x, threads=30){
    
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
    
    # find least variable peaks
    var.peak <- apply(x, 1, var)
    var.peak <- var.peak[order(var.peak, decreasing=F)]
    num.peak <- floor(length(var.peak)*0.5)
    least.var.peaks <- names(var.peak[1:num.peak])
    null.sd <- apply(sp[least.var.peaks,], 1, sd)
    null.mean <- apply(sp[least.var.peaks,], 1, mean)
    fit <- lm(null.sd~null.mean)
    df <- apply(sp, 1, mean)
    fitted.sd <- predict(fit, newdata=data.frame(null.mean=df), type="response")
    
    # iterate over each peak, building a null distribution using the mean
    pvals <- mclapply(seq(1:nrow(sp)), function(z){
        if((z %% 1000) == 0){message(" - iterated over ",z," ACRs ...")}
        qp <- sp[z,]
        qp1 <- qp[is.finite(qp)]
        ave <- mean(qp1, na.rm=T)
        edist <- rnorm(500000, mean=ave, sd=fitted.sd[z])
        #n.fit <- fitdist(edist, "norm")
        #pvals <- pnorm(qp, mean=ave, sd=sd(edist))
        edist <- edist[order(edist, decreasing=F)]
        vals <- lapply(1:length(qp), function(y){
            obs <- qp[y]
            length(edist[edist < obs])/length(edist)
        })
        return(unlist(vals))
    }, mc.cores=threads)
    pvals <- do.call(rbind, pvals)
    colnames(pvals) <- colnames(x)
    rownames(pvals) <- rownames(x)
    return(pvals)
}

head(merged_meta_cpm_information_copied)

three_col_prep <- merged_meta_cpm_information_copied  %>% 
    dplyr::rename("V2" = "geneID")  %>% 
    dplyr::rename("V1" = !!sym(meta_slot_var))  %>% 
    select(-counts)



    # make sure bins/cells are factors
    three_col_prep$V1 <- factor(three_col_prep$V1)
    three_col_prep$V2 <- factor(three_col_prep$V2)


    # convert to sparseMatrix format
    sparse_count_matrix <- Matrix::sparseMatrix(i=as.numeric(three_col_prep$V1),
                              j=as.numeric(three_col_prep$V2),
                              x=as.numeric(three_col_prep$grouped_CPM),
                             dimnames=list(levels(three_col_prep$V1),levels(three_col_prep$V2)))





transposed_ACRs_by_ct <- as.matrix(t(sparse_count_matrix))

transposed_ACRs_by_ct


library("parallel")
final <- specificity(transposed_ACRs_by_ct, 5)

#final$ACR < rownames(final)
#x <- as_tibble(final)  %>% 
#    pivot_longer(names_to = "cell_type", values_to="pval")

ct_specific_acrs_filtered_cell_types <- as_tibble(as.data.frame(as.table(final)))  %>% 
    dplyr::filter(Freq < .01)


counts_ct <- ct_specific_acrs_filtered_cell_types  %>% 
    group_by(Var1) %>% 
    summarise(count_n_cell_types = n()) 

counts_ct <- ct_specific_acrs_filtered_cell_types  %>% 
    group_by(Var1) %>% 
    summarise(count_n_cell_types = n())  %>% 
    arrange(desc(count_n_cell_types))  %>% 
    dplyr::filter(count_n_cell_types == 1)


ct_specific_acrs_filtered_cell_types  %>% 
    dplyr::filter(Var1 %in% counts_ct$Var1)


write_delim(ct_specific_acrs_filtered_cell_types, file = paste0(prefix, ".pvalues.csv"), delim=",")


ct_specific_acrs_filtered_cell_types_joined_names <- ct_specific_acrs_filtered_cell_types  %>% 
    mutate(sc_acr_name = str_c(Var1, Var2, sep = ";"))

ct_specific_acrs_filtered_cell_types_joined_names

combined <- left_join(ct_specific_acrs_filtered_cell_types_joined_names, bed_file_read, by = c("Var1" = "acr_number"))

# Write Cell Type Specific Peaks ALl
#prefix="Zm.V4_annot."
cell_type_acrs_combined <- combined  %>% 
    dplyr::select(chrom, start, stop, sc_acr_name, Freq)

write_tsv(cell_type_acrs_combined, file = paste0(prefix, "all_cts.ACRs.bed"))


# Write Cell Type Specific Peaks by Cell Type 
combined  %>% 
    group_by(Var2)  %>% 
    dplyr::select(chrom, start, stop, sc_acr_name, Freq)  %>% 
    group_walk(~ write_delim(.x, paste0(prefix, .y$Var2, ".cts.ACRs.bed"), delim = "\t"))


head(ct_specific_acrs_filtered_cell_types_joined_names)

#Select Broadly Accessible ACRs and write 
'%ni%' <- Negate("%in%")
broadly_accessible_acrs <- bed_file_read  %>% 
    dplyr::filter(acr_number %ni% ct_specific_acrs_filtered_cell_types_joined_names$Var2)  %>% 
        mutate(sc_acr_name = str_c(acr_number, "broadly_accessible", sep = ";"))

write_tsv(broadly_accessible_acrs, file = paste0(prefix, "broadly_accessible.ACRs.bed"))

#Combine All ACRs and Write
all_acrs_combined <- bind_rows(broadly_accessible_acrs, cell_type_acrs_combined)
write_tsv(all_acrs_combined, file = paste0(prefix, "all_ACRs.classified.bed"))
