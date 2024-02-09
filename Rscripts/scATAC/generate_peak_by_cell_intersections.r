.libPaths(c("/home/jpm73279/R/x86_64-conda-linux-gnu-library/4.1","/home/jpm73279/.conda/envs/R_final_install/lib/R/library", .libPaths()))
library(devtools)
library(reshape2)
library(tidyverse)
load_all('/home/jpm73279/Socrates')


# arguments
#args <- commandArgs(TRUE)
#tn5_bed <- as.character(args[1])
#peak_bed <- as.character(args[1])
#metadata <- as.character(args[3])
#ann <- as.character(args[4])
#chr <- as.character(args[4])
#prefix <- as.character(args[5])

# arguments
args <- commandArgs(TRUE)
tn5_bed  <- as.character(args[1])
peak_bed  <- as.character(args[2])
metadata  <- as.character(args[3])
ann <-  as.character(args[4])
chr <-  as.character(args[5])
prefix <- as.character(args[6])


#tn5_bed <- "/scratch/jpm73279/comparative_single_cell/01.alignments_annotations/zea_mays/Zm_merged_repliates.unique.mpq10.tn5.sorted.bed"
#peak_bed <- "/scratch/jpm73279/comparative_single_cell/07.call.ACRs/Zm_peak_calls/Zm.v4.final/Zm.peaks.annot_v4.500bp_peaks.bed"
#metadata <- "/home/jpm73279/Mendieta_et_al_comparative_single_cell/metrics/annotations/Zm_annot_v4/Zm.leaf_annot.V4.meta.final.txt"
#ann <- "/home/jpm73279/genome_downloads/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
#chr <- "/home/jpm73279/genome_downloads/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.chrom.size"
#prefix <- "test"

#Load Object
Zm.subcluster <- loadBEDandGenomeData(tn5_bed, ann, chr)

Zm.subcluster$meta <- read.table(metadata)
Zm.subcluster$acr <- read.table(peak_bed)

###################################################################################################
###################################################################################################
###################################################################################################
#' generateMatrix
#'
#' This function generates the sparse matrix from equally sized genomic bins or ACRs.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges tileGenome
#' @importFrom IRanges subsetByOverlaps
#'
#' @param obj Object output from findCells or buildMetaData. Required.
#' @param filtered Logical. Whether or not to use the filtered set of cells. Defaults to TRUE.
#' @param windows Integer. Window size to build bins. If the 'peaks' parameter is set to TRUE,
#' this argument is over-ridden.
#' @param peaks Logical. If TRUE, use ACRs to build sparse matrix instead of genomic tiles.
#' Default is set to FALSE.
#' @param blacklist in bed format. If given removes black list regions from either peaks or 
#' generated bins. Default is set to null.
#' Default is set to FALSE.
#' @param verbose Logical. Whether or not to print progress.
#'
#' @rdname generateMatrix
#' @export
#'
generateMatrix_sparse <- function(obj,
                           filtered=T,
                           windows=1000,
                           peaks=FALSE,
                           blacklist=NULL,
                           organelle_scaffolds = NULL,
                           verbose=T){


    # convert tn5 bed to Granges
    tn5.gr <- GRanges(seqnames=as.character(obj$bed$V1),
                      ranges=IRanges(start=as.numeric(obj$bed$V2),
                                     end=as.numeric(obj$bed$V3)),
                      strand=as.character(obj$bed$V5),
                      names=as.character(obj$bed$V4))
    
    
    # Remove Organell Scaffolds if given
    if(is.null(organelle_scaffolds) == FALSE) {
        
        tn5.gr <- dropSeqlevels(tn5.gr, organelle_scaffolds)


    } else {
        tn5.gr <- tn5.gr
    }


    
    # Read in baclist if given
    if(is.null(blacklist) == FALSE) {
        blacklist_r <- read.table(as.character(blacklist))


        blacklist.gr <- GRanges(seqnames=as.character(blacklist_r$V1),
          
                ranges=IRanges(start=as.numeric(blacklist_r$V2),
                                         end=as.numeric(blacklist_r$V3)),
                          names=as.character(blacklist_r$V4))
    } else {
        blacklist.gr <- NULL
    }


    # use filtered barcodes?
    if(filtered){
        use <- obj$meta.v3
    }else{
        use <- obj$meta
    }




    # generate intervals
    if(!peaks){
        # build bins from specified tile length
        chr.seq.lengths <- as.numeric(obj$chr$V2)
        names(chr.seq.lengths) <- obj$chr$V1
        intervals <- tileGenome(chr.seq.lengths, tilewidth=windows, cut.last.tile.in.chrom=TRUE)


        #Remove if black list included
        #Remove procedure learned from: https://www.biostars.org/p/263214/
       if (is.null(blacklist.gr) == FALSE){


            intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),] 
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")




        }else{
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
        }




    }else{


        # generate intervals from ACRs
        intervals <- GRanges(seqnames=as.character(obj$acr$V1),
                             ranges=IRanges(start=as.numeric(obj$acr$V2),
                                            end=as.numeric(obj$acr$V3)))


        if (is.null(blacklist.gr) == FALSE){
            intervals <- intervals[-queryHits(findOverlaps(intervals, blacklist.gr, type="any")),] 
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
        }else{
            regions <- as.data.frame(intervals)
            regions <- paste(regions$seqnames, regions$start, regions$end, sep="_")
        }
    }




    # get intervals overlapping Tn5 sites by barcode
    hits <- as.data.frame(findOverlaps(tn5.gr, intervals))
    df <- data.frame(regions=regions[hits$subjectHits], barcodes=as.character(obj$bed$V4)[hits$queryHits])
    df <- df[!duplicated(df),]
    df$binary <- 1
    colnames(df) <- c("V1","V2","V3")
    
    
    #9/26/2022 include for sake of calculation of isCell 
    # make sure nSites is calculated
    #Integration Sites
    a <- df
    
    #Meta data to interset
    b <- use
    a$V1 <- factor(a$V1)
    a$V2 <- factor(a$V2)


    #Generate sparse matrix
    a <- Matrix::sparseMatrix(i=as.numeric(a$V1),
                              j=as.numeric(a$V2),
                              x=as.numeric(a$V3),
                              dimnames=list(levels(a$V1),levels(a$V2)))


    # align barcodes
    both <- intersect(rownames(b), colnames(a))
    a <- a[,both]
    b <- b[both,]


    # make sure nSites is calculated
    b$nSites   <- Matrix::colSums(a)
    b$log10nSites <- log10(b$nSites)


    # return
    obj$counts <- df
    obj$meta <- b 


    return(obj) 
    
    }



message("Generating Sparse Matrix")

Zm.subcluster.vasculature <- generateMatrix_sparse(Zm.subcluster, 
                                              filtered=FALSE,
                                              peaks = TRUE,
                                              verbose=FALSE)

#Zm.subcluster.vasculature <- convertSparseData(Zm.subcluster.vasculature, verbose = TRUE)

message("Sparse Matrix Done! Preparing to write to ouptput")
    sparse_count_matrix <- Zm.subcluster.vasculature$counts
    # make sure bins/cells are factors
    sparse_count_matrix$V1 <- factor(sparse_count_matrix$V1)
    sparse_count_matrix$V2 <- factor(sparse_count_matrix$V2)




    # convert to sparseMatrix format
    sparse_count_matrix <- Matrix::sparseMatrix(i=as.numeric(sparse_count_matrix$V1),
                              j=as.numeric(sparse_count_matrix$V2),
                              x=as.numeric(sparse_count_matrix$V3),
                             dimnames=list(levels(sparse_count_matrix$V1),levels(sparse_count_matrix$V2)))






message("Saving Output")
saveRDS(sparse_count_matrix, file = paste0(prefix,".peaks_by_intersections.rds"))
