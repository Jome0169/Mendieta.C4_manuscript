######################
## Cicero trajectories
##updating 051321 set the argument of different s and window size and distance_constraint
##updation 120120 this version will help us to correct FDR based the FDR function
##we need to add the negative score case

library(cicero)
library(Matrix)
library(parallel)
library(doSNOW)
library(methods)
library(tcltk)
library(iterators)
library(itertools)

# default variables
threads <- 1
dims <- 50

#commandline arguments
args = commandArgs(TRUE)

#read in commandline arguments
threads  <- as.numeric(args[1])
output   <- as.character(args[2])
input    <- as.character(args[3])
metafile <- as.character(args[4])
genome_fl <- as.character(args[5])
gene_bed_fl <- as.character(args[6])
output_dir <- as.character(args[7])

##other parameter setting
s <- as.numeric(args[8])
window <- as.numeric(args[9])
distance_constraint <- as.numeric(args[10])

##create dir to store the temp file
raw_and_shuffled_dir <- paste0(output_dir,'/raw_and_shuffled_dir')
dir.create(raw_and_shuffled_dir)
gene_acc_dir <- paste0(output_dir,'/gene_acc_dir')
dir.create(gene_acc_dir)
store_temp_gene_acc_dir <- paste0(output_dir,'/store_temp_gene_acc_dir')
dir.create(store_temp_gene_acc_dir)

#######################
## Load data, functions
# add annotation
loadMeta <- function(cds, jac){
  
    # svd and raw
    ids <- colnames(cds)
    print(head(ids, n=10))
    jac <- jac[ids,]
    print(head(jac, n=10))
    jac$UMAP1 <- jac$umap1
    jac$UMAP2 <- jac$umap2
    ids.2 <- rownames(jac)
    print(head(ids.2, n=10))
    cds <- cds[,colnames(cds) %in% ids.2]
    umap.d <- t(jac[,c("UMAP1","UMAP2")])
    
    # UMAP output
    cds@reducedDimA <- umap.d
    cds@reducedDimS <- umap.d
    cds@dim_reduce_type <- "UMAP"
    ##pData is kind of project that we can add the items on it
    pData(cds)$tissue_cluster <- jac$LouvainClusters
    pData(cds)$Cluster <- as.factor(jac$LouvainClusters)
    pData(cds)$preclust <- as.factor(pData(cds)$LouvainClusters)
    
    # return data
    return(cds)
}


# load data
message(" ... Loading data")
a <- read.table(input)
zm <- read.table(genome_fl)
genes <- read.table(gene_bed_fl)
##change the genes column name
colnames(genes) <- c('chromosome','start','end','strand','feature','gene','transcript','symbol')
##updating 020421 debug chagne the chromosome name of genes
genes$chromosome <- genes$chromosome

#print(genes)

meta <- read.table(metafile)

################
## Cluster cells 
# create cicero object
message(" ... Creating CDS")
##V1 is the peak name and V2 is the cell name 
##have an intersection for these two files
a <- a[as.character(a$V2) %in% rownames(meta),]
#print(a)
shuf <- a
##shuffle the peaks and the cell id
shuf$V1 <- shuf$V1[sample(length(shuf$V1))]
shuf$V2 <- shuf$V2[sample(length(shuf$V2))]
##the make_atac_cds is from the cicero
cds <- make_atac_cds(a, binarize=T)
#print(cds)
##generate shuffed cds
shufcds <- make_atac_cds(shuf, binarize=T)

##calculate the col and row sum of cds
c.colSums <- Matrix::colSums(exprs(cds))
c.rowSums <- Matrix::rowSums(exprs(cds))
s.colSums <- Matrix::colSums(exprs(shufcds))
s.rowSums <- Matrix::rowSums(exprs(shufcds))

message("... Adding Meta data")
##add the metadata infor to the cds object
##aim of pData is store the pheno data such as the meta information
cds <- cds[,colnames(cds) %in% rownames(meta)]
shufcds <- shufcds[,colnames(cds) %in% rownames(meta)]

print("Here are the number of colnames in CDS which correspond to ROW names in meta")
print(dim(cds))

pData(cds) <- meta[colnames(exprs(cds)),]
pData(shufcds) <- meta[colnames(exprs(shufcds)),]
cds <- cds[Matrix::rowSums(exprs(cds))>0,]
cds <- cds[,Matrix::colSums(exprs(cds))>0]
shufcds <- shufcds[Matrix::rowSums(exprs(shufcds))>0,]
shufcds <- shufcds[,Matrix::colSums(exprs(shufcds))>0]

##now the cds is an object
cds <- detectGenes(cds) ##Counts how many cells each feature in a CellDataSet object that are detectably expressed above a minimum threshold.
shufcds <- detectGenes(shufcds) 
cds <- estimateSizeFactors(cds)
shufcds <- estimateSizeFactors(shufcds)

# load results from jaccard
message(" ... Loading Jaccard-based clustering results and reduced dimensions")
cds <- loadMeta(cds, meta)
shufcds <- loadMeta(shufcds, meta)

message(" ... Loaded Meta")

####################################################
## Estimate connections, modules and gene activities			       

##add a function of run_cicero
run_cicero_new <- function(cds,
                       genomic_coords,
                       window = 500000,
                       silent=FALSE,
                       sample_num = 100,
                       distance_constraint = 250000,
                       s=0.75) {
  # Check input
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(is.logical(silent))
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(sample_num))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  
  if (!silent) print("Starting Cicero")
  if (!silent) print("Calculating distance_parameter value")
  distance_parameters <- estimate_distance_parameter(cds, window=window,
                                                     maxit=100, sample_num = sample_num,
                                                     distance_constraint = distance_constraint,
                                                     s=s,
                                                     distance_parameter_convergence = 1e-22,
                                                     genomic_coords = genomic_coords)
  
  mean_distance_parameter <- mean(unlist(distance_parameters))
  
  if (!silent) print("Running models")
  cicero_out <-
    generate_cicero_models(cds,
                           distance_parameter = mean_distance_parameter,
                           window = window,
                           genomic_coords = genomic_coords)
  
  if (!silent) print("Assembling connections")
  all_cons <- assemble_connections(cicero_out, silent=silent)
  
  if (!silent) print("Done")
  all_cons
}



message(" ... Initializing per cluster cicero run")

##run for each cluster
# run cicero to get co-accessible sites and modules BY CLUSTER
meta2 <- pData(cds)
meta2$Cluster <- as.character(meta2$Cluster)
print(table(meta2$Cluster))
clusts <- unique(meta2$Cluster)
cell_ids <- c()

# iterate
its <- 0

# foreach parameters
##this function aims for the parallel analysis
cl <- makeSOCKcluster(threads)
registerDoSNOW(cl)
tasks <- length(clusts)
pb <- txtProgressBar(max = tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
package.labs <- c("cicero", "Matrix")
message(" ... Initializing per cluster cicero run")

##store the temp files
saveRDS(meta2,paste0(store_temp_gene_acc_dir,'/meta2.rds'))
saveRDS(cds,paste0(store_temp_gene_acc_dir,'/cds.rds'))
saveRDS(shufcds,paste0(store_temp_gene_acc_dir,'/shufcds.rds'))
saveRDS(clusts,paste0(store_temp_gene_acc_dir,'/clusts.rds'))

# run in parallel
gact <- list()
gact <- foreach(i=clusts, .combine='c', .packages=package.labs, .options.snow=opts) %dopar% {

    # get umap coordinates and make cicero CDS
    message("###--- Creating cicero object, cluster",i, " ---###")
    ids <- rownames(meta2[meta2$Cluster==i,])
    index.keep <- colnames(exprs(cds)) %in% ids ##this is to extract the ids information
    s.cds <- cds[,index.keep]

    # only consider sites accessible in at least 1% of cells in cluster
    ##the 0.00 should be 0.01?
    s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>(ncol(exprs(s.cds))*0.00),]
    s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
    ##exprs indicate the matrix of cds
    print(head(exprs(s.cds)[,1:5]))
    message(" - number of sites for cluster ", i, " = ", nrow(s.cds))

    # get UMAP coordinates
    umap_coords <- t(reducedDimA(s.cds))
    umap_coords <- umap_coords[colnames(s.cds),]
    message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
    rownames(umap_coords) <- colnames(exprs(s.cds))
    ##we will generate cicero cds based on the previous developed umap coordinates
    ##Next, we access the tSNE coordinates from the input CDS object where they are stored by Monocle and run
    cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=20)

    # run cicero (connections)
    message(" ... Running cicero")
    ##try window 300000 pre is 100000
    conns <- run_cicero_new(cicero_cds, zm, window=window, sample_num=100,distance_constraint=distance_constraint, s=s)
    
    #window <- as.numeric(args[10])
    #distance_constraint <- as.numeric(args[11])
    
    #?run_cicero
    
    # write results to disk
    write.table(conns, file=paste(raw_and_shuffled_dir,"/bycluster",i,".",output,".cicero.loops.txt",sep=""),
                sep="\t",quote=F, col.names=F, row.names=F)
}
close(pb)
stopCluster(cl)


###################################################################################################
## SHUFFLE ----------------------------------------------------------------------------------------
###################################################################################################
##The idea is to find the cut-off in the shuffled/random data that removes 95% of coACRs. Then use this cut-off on the real data
##the shuffled/random coACRs are “true FP”

# for each cluster
meta2 <- pData(shufcds)
meta2$Cluster <- as.character(meta2$Cluster)
#print(table(meta2$Cluster))
clusts <- unique(meta2$Cluster)
cell_ids <- c()

# iterate
its <- 0

# foreach parameters
cl <- makeSOCKcluster(threads)
registerDoSNOW(cl)
tasks <- length(clusts)
pb <- txtProgressBar(max = tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
package.labs <- c("cicero", "Matrix")
message(" ... Initializing shuffled per cluster cicero run")

## shuffled ##
gact <- list()
gact <- foreach(i=clusts, .combine='c', .packages=package.labs, .options.snow=opts) %dopar% {

    # get umap coordinates and make cicero CDS
    message("###--- Creating cicero object, cluster",i, " ---###")
    ids <- rownames(meta2[meta2$Cluster==i,])
    index.keep <- colnames(exprs(shufcds)) %in% ids
    s.cds <- shufcds[,index.keep]

    # only consider sites accessible in at least 1% of cells in cluster
    s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>(ncol(exprs(s.cds))*0.00),]
    s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
    #print(head(exprs(s.cds)[,1:5]))
    message(" - number of sites for cluster ", i, " = ", nrow(s.cds))

    # get UMAP coordinates
    umap_coords <- t(reducedDimA(s.cds))
    umap_coords <- umap_coords[colnames(s.cds),]
    message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
    rownames(umap_coords) <- colnames(exprs(s.cds))
    cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=20)

    # run cicero (connections)
    message(" ... Running cicero")
    ##try window 300000 pre is 100000
    conns <- run_cicero_new(cicero_cds, zm, window=window, sample_num=100,distance_constraint=distance_constraint, s=s)
   
    # write results to disk
    write.table(conns, file=paste(raw_and_shuffled_dir,"/shuffled_bycluster",i,".",output,".cicero.loops.txt",sep=""),
                sep="\t",quote=F, col.names=F, row.names=F)

}
close(pb)
stopCluster(cl)


########################  
## COMPUTE GENE ACTIVITY

##add the positive and negative conns for the FDR
##updation 120120 
getFDR <- function(obs,
                   exp,
                   fdr=0.05,
                   grid=100,
                   verbose=F){
  
  ##updation 120120 we change the x to V3
  # split into +/-
  n.exp <- subset(exp, exp$V3 < 0)
  p.exp <- subset(exp, exp$V3 > 0)
  
  # get counts
  pos.nexp <- nrow(p.exp)
  neg.nexp <- nrow(n.exp)
  pos.nobs <- nrow(subset(obs, obs$V3 > 0))
  neg.nobs <- nrow(subset(obs, obs$V3 < 0))
  if(verbose){message(" - number expected links = (+) ",pos.nexp, " | (-) ",neg.nexp)}
  if(verbose){message(" - number observed links = (+) ",pos.nobs, " | (-) ",neg.nobs)}
  
  # generate range of thresholds
  p.vals <- seq(from=0, to=1, length.out=grid)
  n.vals <- seq(from=0, to= -1, length.out=grid)
  
  # iterate over grid
  if(verbose){message(" - scanning positive thresholds ...")}
  p.thresh <- c()
  for(i in p.vals){
    num.exp <- sum(p.exp$V3 > as.numeric(i))
    c.fdr <- num.exp/(pos.nexp)
    if(is.na(c.fdr)){
      c.fdr <- 0
    }
    p.thresh <- c(p.thresh, c.fdr)
    if(verbose){message(" - (+) correlation threshold = ", i, " | FDR = ", c.fdr)}
    if(c.fdr < fdr){
      break
    }
  }
  if(verbose){message(" - scanning negative thresholds ...")}
  n.thresh <- c()
  for(i in n.vals){
    num.exp <- sum(n.exp$V3 < as.numeric(i))
    c.fdr <- num.exp/(neg.nexp)
    if(is.na(c.fdr)){
      c.fdr <- 0
    }
    n.thresh <- c(n.thresh, c.fdr)
    if(verbose){message(" - (-) correlation threshold = ", i, " | FDR = ", c.fdr)}
    if(c.fdr < fdr){
      break
    }
  }
  
  # select cut-offs
  p.threshold <- min(p.vals[which(p.thresh <= fdr)])
  n.threshold <- max(n.vals[which(n.thresh <= fdr)])
  
  # filter
  obs <- subset(obs, obs$V3 > p.threshold | obs$V3 < n.threshold)
  
  # verbose number of +/- linkages
  pos.links <- nrow(subset(obs, obs$V3 > 0))
  neg.links <- nrow(subset(obs, obs$V3 < 0))
  if(verbose){message(" - found ",pos.links, " + and ", neg.links," - ACR-ACR links ...")}
  
  # return
  return(obs)
}


# for each cluster
meta2 <- pData(cds)
meta2$Cluster <- as.character(meta2$Cluster)
#print(table(meta2$Cluster))
clusts <- unique(meta2$Cluster)
cell_ids <- c()

# iterate
its <- 0

# foreach parameters
message(" ... Initializing per cluster cicero run - GENE ACTIVITY")
gascores <- mclapply(clusts, function(x){

    # load connections
    id.true <- paste(raw_and_shuffled_dir,"/bycluster",x,".",output,".cicero.loops.txt",sep="")
    id.false <- paste(raw_and_shuffled_dir,"/shuffled_bycluster",x,".",output,".cicero.loops.txt",sep="")
    t.conns <- read.table(id.true)
    s.conns <- read.table(id.false)
    
    # filter loops --------------------------------------------------------------------------------
    
    # empty vector
    #b.sub <- c()
    #lims <- seq(from=0, to=0.99, length.out=100)
    
    # find cut-off
    ##j is from the 0 to 0.99
    ##this is the co-accessibility score
    ##we want to detect the minmum positive co-acc score
    ##such as threshold we obtain is 0.03
    ##so the connections with lower than 0.03 will be filtered out
    ##how to idenitfy the 0.03,
    ##we allow the fdr < 0.05 which means the conn acc is 0.03
    #for(j in lims){
    #    b.sub <- c(b.sub, nrow(subset(s.conns, s.conns$V3 >= j)))
    #}
    ##b.sub is generated from different threshold level
    ##if we get higher threshold level, the fdr value is getting lower 
    ##fdr
    #[1] 4.776554e-01 1.053479e-01 5.821455e-02 4.090530e-02 3.171600e-02
    #[6] 2.548796e-02 2.118709e-02 1.759128e-02 1.479453e-02 1.245608e-02
    #[11] 1.057592e-02 8.636999e-03 7.320885e-03 6.263293e-03 5.240955e-03
    
  
    ##the main idea is that the s.conns is the FP, we will remove 95% of them, and the left 5% is the T,
    ##it looks like we have 500 FP and these 500 FP has threshold like 0.01 10, 0.02, 20, 0.03 25
    ##5% of the 500 FP is the 25, which means the 25 corresponding 0.03 will be the threshold 
    
    #fdr <- b.sub/nrow(s.conns)
    #threshold <- min(lims[which(fdr < 0.05)])
    #message(" - threshold = ", threshold, " | ", id.true)
    
    ##updation 120120 use the getFDR function
    # filter loops
    a.sub <- getFDR(t.conns, s.conns, fdr=0.05, verbose=F)
    #a.sub <- subset(t.conns, t.conns$V3 >= threshold)
    id <- gsub("bycluster", "filtered", id.true)
    write.table(a.sub, file=paste0(id), quote=F, row.names=F, col.names=F, sep="\t")
    colnames(a.sub) <- c("Peak1", "Peak2", "coaccess")
    
    # get gene activity scores --------------------------------------------------------------------
    message("--- estimating gene activity scores for cluster ",x)
    ##cell ids
    ids <- rownames(meta2[meta2$Cluster==x,])
    index.keep <- colnames(exprs(cds)) %in% ids
    s.cds <- cds[,index.keep]
    
    # only consider sites accessible in at least 1% of cells in cluster
    s.cds <- s.cds[Matrix::rowSums(exprs(s.cds))>0,]
    s.cds <- s.cds[,Matrix::colSums(exprs(s.cds))>0]
    print(head(exprs(s.cds)[,1:2]))
    message(" - number of sites for cluster ", x, " = ", nrow(s.cds))
    
    # get UMAP coordinates
    ##we may not use this information since the cicero_cds is not used in the downstream analysis
    umap_coords <- t(reducedDimA(s.cds))
    umap_coords <- umap_coords[colnames(s.cds),]
    message("# UMAP coords = ", nrow(umap_coords), " | # cells = ", ncol(s.cds))
    rownames(umap_coords) <- colnames(exprs(s.cds))
    cicero_cds <- make_cicero_cds(s.cds, reduced_coordinates=umap_coords, k=20)
    
    ##estimate the gene activity
    ##it will generate the promoter region 1kb to the real start
    ##the real end is the start of the transcript
    # estimate gene activity
    message(" ... Estimating gene activity scores")
    pos <- subset(genes, strand == "+")
    pos <- pos[order(pos$start),]
    pos <- pos[!duplicated(pos$gene),]
    pos$end <- pos$start

    pos$start <- pos$start - 1000
    pos$start[pos$start<0] <- 0
    
    neg <- subset(genes, strand == "-")
    neg <- neg[order(neg$start, decreasing = T),] 
    neg <- neg[!duplicated(neg$gene),] 
    neg$start <- neg$end
    neg$end <- neg$end + 1000
    
    # merge
    gene_ann2 <- rbind(pos, neg)
    gene_ann2 <- gene_ann2[,c(1:3, 8)]
    gene_ann2 <- gene_ann2[order(gene_ann2$start, decreasing=F),]
    colnames(gene_ann2)[4] <- "gene"

    message("- Here are the annotations for the annotate_by_CDS ...")
    #print(head(gene_ann3, n = 10)
    #annotate genes
    message("- annotate genes by peaks ...")
    s.cds <- annotate_cds_by_site(s.cds, gene_ann2, all=F, verbose = TRUE)

    # estimate un-normalized activity
    message("- build gene activity matrix ... ")
    ##s.cds is the object
    ##a.sub is the conn file
    unnorm_ga <- build_gene_activity_matrix(s.cds, a.sub, dist_thresh = 100000, coaccess_cutoff = 0)
    unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]
  
    # gene activity per cluster
    num_genes <- pData(s.cds)$num_genes_expressed
    names(num_genes) <- row.names(pData(s.cds))
    
    # normalize
    cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
    geneact <- as.data.frame(summary(cicero_gene_activities))
    geneact$i <- rownames(cicero_gene_activities)[geneact$i]
    geneact$j <- colnames(cicero_gene_activities)[geneact$j]
    
    # output
    write.table(geneact, file=paste(gene_acc_dir,"/filtered",x,".",output,".cicero.geneActivity.txt",sep=""),
                sep="\t",quote=F, col.names=F, row.names=F)
    
    # return
    return(unnorm_ga)
}, mc.cores=threads)




message(" ... Gather All Gene Activity Scores together. This should be quick!")
gen_name <- function(name_1, gene_acc_dir, output_loc){
    final <- paste0(gene_acc_dir,"/filtered",name_1,".",output_loc,".cicero.geneActivity.txt")
}


file_names <- lapply(clusts, gen_name, gene_acc_dir, output)
merged_data_list <- lapply(file_names, read.table)
final_merged<- do.call(rbind,merged_data_list)

write.table(final_merged, file=paste(output,".all.cicero.geneActivity.txt",sep=""),
               sep="\t",quote=F, col.names=F, row.names=F)



