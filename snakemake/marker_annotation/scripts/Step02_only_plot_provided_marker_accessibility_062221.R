###################################################################################################
###################################################################################################
##                                                                                               ##
##                  functions for plotting accessibility of markers from cicero                  ##
##                                                                                               ##
###################################################################################################
###################################################################################################

##only plot the defined markers

# load libraries
library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
library(varistran)
library(edgeR)
library(parallel)
library(png)


###################################################################################################
###################################################################################################
###################################################################################################

args <- commandArgs(trailingOnly=T)

all_gene_impute_acc <- as.character(args[1])
plot_each_CT <- as.character(args[2])
GAobj_rds <- as.character(args[3]) ##GAobj.rds
markers <- as.character(args[4])
output_dir <- as.character(args[5])

lim <- as.numeric(args[6])




##create eachCT marker dir
marker_opt_dir <- paste0(output_dir,'/opt_store_eachCT_marker_dir')
if (!dir.exists(marker_opt_dir)){
  dir.create(marker_opt_dir)
} else {
  print("Dir already exists!")
}

##prepare the input
dat <- readRDS(GAobj_rds)
all.b <- dat$b
impute.activity <- readRDS(all_gene_impute_acc)
marker.info <- read.table(markers, header=T)
rownames(marker.info) <- marker.info$geneID



plot.act.scores    <- function(df,output_dir,
                               acts=acts, 
                               info=NULL, 
                               top=NULL,
                               logT=F,
                               marker.dist=NULL,
                               outname="markerActivityScores.pdf", 
                               lim=0.95){
  
  # prep data
  df <- df[rownames(df) %in% colnames(acts),]
  acts <- acts[,which(rownames(df) %in% colnames(acts))]
  
  # reorder rows
  rownames(info) <- info$geneID
  info <- info[order(info$type),]
  info.genes <- rownames(info)
  act.genes <- rownames(acts)
  rd.cells <- rownames(df)
  
  # common genes
  common <- intersect(info.genes, act.genes)
  info <- info[which(rownames(info) %in% common),]
  info.ordered <- rownames(info)
  sub.scores <- acts[info.ordered,]
  gids <- info.ordered
  
  # setup plot size
  nrows <- ceiling(length(gids)/6)
  totals <- nrows*6
  ratio <- nrows/6
  
  # params
  png(file=paste0(output_dir,'/',outname), width=12, height=ratio*12, units="in", res=500, type="cairo")
  layout(matrix(c(1:totals), ncol=6, byrow=T))
  par(mar=c(2,2,1,1))
  
  # adjust cluster IDs
  message("begin plotting pre-defined markers...")
  for (i in 1:length(gids)){
    
    # copy meta data
    gene.index <- which(rownames(sub.scores) == gids[i])
    acv <- sub.scores[gene.index,]
    
    # set up plot cols/sizes
    orderRow <- order(acv, decreasing=F)
    #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
    #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
    #cols <- inferno(100)
    #cols <- plasma(100)
    cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
    acv <- as.numeric(acv[orderRow])
    if(logT==T){
      acv <- log2(acv+1)
    }
    df2 <- df[orderRow,]
    acv[is.na(acv)] <- 0
    acv[is.infinite(acv)] <- 0
    upper.lim <- quantile(acv, lim)
    acv[acv > upper.lim] <- upper.lim
    if(!is.null(marker.dist)){
      message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
      colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
    }else{
      min.acv <- min(acv) - (1e-6*min(acv))
      max.acv <- max(acv) + (1e-6*max(acv))
      message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
      if(min.acv == max.acv){
        next
      }
      colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
    }
    colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
    colvec[is.na(colvec) & acv == 0] <- cols[1]
    sizes <- rescale(acv, c(0.25, 0.3))
    
    # plot
    plot(df2$umap1, df2$umap2, col=colvec,
         main=paste(info$name[i],info$type[i],sep="-"),
         #main=info$name[i],
         cex.main=0.6,
         xlab="", ylab="", bty="n",
         xaxt="n", yaxt="n", pch=16, cex=0.25)
    
  }
  
  # turn device off
  dev.off()
  
}

plot.act.scores_eachCT    <- function(marker_opt_dir,df, 
                                      acts=acts, 
                                      info=NULL, 
                                      top=NULL,
                                      logT=F,
                                      marker.dist=NULL,
                                      outname="markerActivityScores.pdf", 
                                      lim=0.95){
  
  ##updating 062221
  ##save the each cell type to a folder
  unique_type_list <- unique(info$type)
  
  for (i in 1:length(unique_type_list)) {
    
    target_celltype <- unique_type_list[i]
    
    message(paste0(" - analyze target celltype ",target_celltype))
    
    # prep data
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,which(rownames(df) %in% colnames(acts))]
    
    # reorder rows
    rownames(info) <- info$geneID
    info <- info[order(info$type),]
    
    ##updating 062221
    ##select the target cells
    info <- info[info$type == target_celltype,]
    
    info.genes <- rownames(info)
    act.genes <- rownames(acts)
    rd.cells <- rownames(df)
    
    # common genes
    common <- intersect(info.genes, act.genes)
    info <- info[which(rownames(info) %in% common),]
    info.ordered <- rownames(info)
    sub.scores <- acts[info.ordered,]
    gids <- info.ordered
    
    # setup plot size
    nrows <- ceiling(length(gids)/6)
    totals <- nrows*6
    ratio <- nrows/6
    
    # params
    #png(file=outname, width=12, height=ratio*12, units="in", res=500, type="cairo")
    #png(file=paste0(marker_opt_dir,'/',target_celltype,'.png'), width=12, height=ratio*12, units="in", res=500, type="cairo")
    
    outname = paste0(marker_opt_dir,'/',target_celltype,'.pdf')
    png(file=outname, width=12, height=ratio*12, units="in", res=500, type="cairo")
    
    layout(matrix(c(1:totals), ncol=6, byrow=T))
    par(mar=c(2,2,1,1))
    
    # adjust cluster IDs
    message("begin plotting pre-defined markers...")
    for (i in 1:length(gids)){
      
      # copy meta data
      gene.index <- which(rownames(sub.scores) == gids[i])
      acv <- sub.scores[gene.index,]
      
      # set up plot cols/sizes
      orderRow <- order(acv, decreasing=F)
      #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
      #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
      #cols <- inferno(100)
      #cols <- plasma(100)
      cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.5)(100)
      acv <- as.numeric(acv[orderRow])
      if(logT==T){
        acv <- log2(acv+1)
      }
      df2 <- df[orderRow,]
      acv[is.na(acv)] <- 0
      acv[is.infinite(acv)] <- 0
      upper.lim <- quantile(acv, lim)
      acv[acv > upper.lim] <- upper.lim
      if(!is.null(marker.dist)){
        message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
        colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
      }else{
        min.acv <- min(acv) - (1e-6*min(acv))
        max.acv <- max(acv) + (1e-6*max(acv))
        message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
        if(min.acv == max.acv){
          next
        }
        colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
      }
      colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
      colvec[is.na(colvec) & acv == 0] <- cols[1]
      sizes <- rescale(acv, c(0.25, 0.3))
      
      # plot
      plot(df2$umap1, df2$umap2, col=colvec,
           main=paste(info$name[i],info$type[i],sep="-"),
           #main=info$name[i],
           cex.main=0.6,
           xlab="", ylab="", bty="n",
           xaxt="n", yaxt="n", pch=16, cex=0.25)
      
    }
    
    # turn device off
    dev.off()
    
  }
}



if (plot_each_CT == 'yes'){
  ##plot the markers in each folder
  ##and the outname will not work
  plot.act.scores_eachCT(marker_opt_dir,
                         all.b,
                         acts=impute.activity,
                         info=marker.info,
                         logT=F,
                         lim=lim,
                         marker.dist=NULL,
                         outname=paste0("combined.impute.known.Markers.png"))
  
}else{
  plot.act.scores(all.b,output_dir,
                  acts=impute.activity,
                  info=marker.info,
                  logT=F,
                  lim=lim,
                  marker.dist=NULL,
                  outname=paste0("combined.impute.known.Markers.png"))

}















