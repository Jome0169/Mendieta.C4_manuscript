###############################################
##plotting accessibility of markers from cicero
##updation 120420 change the loading function to the original one and do not gzip
##updation 120420 change the Louv cluster since we conduct a step to remove the genotype effects

##this script will generate two types of normalized gene acc
##1) normalized acc: 1, scale based on the Mclust threshold 2, normalized by the Geometric mean
##2) impute normalized acc: 1,2, 3, use the smooth function to have an normalization

# libs
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
library(phytools)

# load data
loadData           <- function(meta, pcs, geneact, markers){
    
    # verbose
    message(" - loading data ...")
    
    # load
    message("Loading meta data")
    b <- read.table(meta)

    message("Validating Clusters")
    b <- validateCluster(b)

    message("Reading Marker info")
    marker.info <- read.table(markers, header=T)
    rownames(marker.info) <- marker.info$geneID

    message("Reading PCs")
    h.pcs <- read.table(pcs)
    #activity <- read.table(gzfile(geneact))

    message("Reading Activitiy from cicero")
    activity <- read.table(geneact,stringsAsFactors = T)
    b <- b[rownames(b) %in% rownames(h.pcs),]
    h.pcs <- h.pcs[rownames(h.pcs) %in% rownames(b),]
    
    # reformat sparse
    activity <- sparseMatrix(i=as.numeric(activity$V1),
                             j=as.numeric(activity$V2),
                             x=as.numeric(activity$V3),
                             dimnames=list(levels(activity$V1), levels(activity$V2)))
    #activity <- readRDS(geneact)
    ##geneact is 'cicero_gene_activities_sparse_mtx.rds' for TEs
    
    activity <- activity[,colnames(activity) %in% rownames(b)]
    activity <- activity[,Matrix::colSums(activity>0) > 0]
    activity <- activity[Matrix::rowSums(activity>0) > 0,]
    activity <- activity[,Matrix::colSums(activity>0) > 0]
    b <- b[colnames(activity),]
    h.pcs <- h.pcs[colnames(activity),]
    
    # sanity check
    sanity <- 0
    if(sanity > 0){
	print(head(b))
        print(head(marker.info))
        print(head(h.pcs[,1:5]))
        print(head(activity[,1:5]))
    }

    # output
    return(list(h.pcs=h.pcs, b=b, activity=activity, marker.info=marker.info))
    
}
validateCluster    <- function(df){
    if(! "Cluster" %in% colnames(df)){
        if("tissue_cluster" %in% colnames(df)){
            df$Cluster <- paste0("cluster", as.numeric(as.factor(df$tissue_cluster)))
            return(df)
        }else if("LouvainClusters" %in% colnames(df)){
            df$Cluster <- paste0("cluster",df$LouvainClusters)
            return(df)
        }
    }else{
        if(length(grepl("cluster",df$Cluster)==T)==0){
            df$Cluster<- paste0("cluster",df$Cluster)
            return(df)
        }else{
            return(df)
        }
    }
}

# tfidf
safe_tfidf         <- function(tf, idf,  block_size=2000e6){
    result = tryCatch({
        result = tf * idf
        result
    }, error = function(e) {
        options(DelayedArray.block.size=block_size)
        DelayedArray:::set_verbose_block_processing(TRUE)
        
        tf = DelayedArray(tf)
        idf = as.matrix(idf)
        
        result = tf * idf
        result
    })
    return(result)
}
tfidf              <- function(bmat, frequencies=F, log_scale_tf=T, scale_factor=100000){
    
    # Use either raw counts or divide by total counts in each cell
    if (frequencies) {
        # "term frequency" method
        tf = t(t(bmat) / Matrix::colSums(bmat))
    } else {
        # "raw count" method
        tf = bmat
    }
    
    # Either TF method can optionally be log scaled
    if (log_scale_tf) {
        if (frequencies) {
            tf@x = log1p(tf@x * scale_factor)
        } else {
            tf@x = log1p(tf@x * 1)
        }
    }
    
    # IDF w/ "inverse document frequency smooth" method
    idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
    
    # TF-IDF
    tf_idf_counts = safe_tfidf(tf, idf)
    rownames(tf_idf_counts) = rownames(bmat)
    colnames(tf_idf_counts) = colnames(bmat)
    return(Matrix(tf_idf_counts, sparse=T))
}

# processing
sparse_apply       <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
    if (convert_to_dense){
        if (MARGIN == 1){
            Sp_X <- Matrix::t(Sp_X)
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(as.matrix(Sp_X[,i]), ...)
            }, FUN, ...)
        }else{
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(as.matrix(Sp_X[,i]), ...)
            }, FUN, ...)
        }
    }else{
        if (MARGIN == 1){
            Sp_X <- Matrix::t(Sp_X)
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(Sp_X[,i], ...)
            }, FUN, ...)
        }else{
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(Sp_X[,i], ...)
            }, FUN, ...)
        }
    }
    
    if(MARGIN == 1){
        new <- Matrix(t(as.matrix(do.call(cbind, res))), sparse=T)
        return(new)
    }else{
        new <- Matrix(as.matrix(do.call(cbind, res)), sparse=T)
        return(new)
    }
    
}
normalizeSF        <- function(adj.act, verbose=F){
    
    # verbose
    if(verbose==T){
        message(" - unnormalized activities:")
        print(head(adj.act[,1:5]))
    }
    
    # get per sample norm factors
    norm.factor <- estimate_sf_sparse(adj.act)
    if(verbose==T){
        message(" - normalization factors:")
        print(head(norm.factor))
    }
    
    # normalized counts
    ##this is a typical way to use the Diagonal to have a normalization for the data
    norm.act <- adj.act %*% Diagonal(x=1/norm.factor)
    colnames(norm.act) <- colnames(adj.act)
    rownames(norm.act) <- rownames(adj.act)
    if(verbose==T){
        message(" - normalized activities:")
        print(head(norm.act[,1:5]))
    }
    
    # output
    return(list(norm.act=norm.act, norm.factor=norm.factor))
}
estimate_sf_sparse <- function(counts, round_exprs=T, method="mean-geometric-mean-total"){
    if (round_exprs)
        counts <- round(counts)
    
    ##calculate the total for each cell since the column name is the cell name
    ##why we use the sf?
    if(method == 'mean-geometric-mean-total') {
        cell_total <- Matrix::colSums(counts)
        sfs <- cell_total / exp(mean(log(cell_total)))
    }else if(method == 'mean-geometric-mean-log-total') {
        cell_total <- Matrix::colSums(counts)
        sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
    }
    
    sfs[is.na(sfs)] <- 1
    sfs
}
shannon.entropy    <- function(p){
    if (min(p) < 0 || sum(p) <= 0) return(NA)
    p.norm <- p[p>0]/sum(p)
    -sum(log2(p.norm)*p.norm)
}
gene.score         <- function(data=data, markers=markers, rd=rd){
    
    # reformat data into sparseMatrix
    if(is.data.frame(data)){
        message("reformating data into sparseMatrix ...")
        data <- sparseMatrix(i=as.numeric(data$V1),
                             j=as.numeric(data$V2),
                             x=as.numeric(data$V3),
                             dimnames=list(levels(data$V1), levels(data$V2)))
    }
    
    # reorder matrices
    data <- data[,rownames(rd)]
    cellave <- Matrix::colSums(data)
    
    # collect unique genes
    genes <- unique(markers$V4)
    mat <- matrix(ncol=length(genes), nrow=nrow(rd), dimnames=list(rownames(rd), genes))
    cnts <- 0
    
    # iterate over unique genes
    for(i in genes){
        cnts <- cnts+1
        message("estimating accessibility for gene: ",i, " ...")
        peak.sub <- subset(markers, markers$V4==i)
        gene.peak <- as.character(droplevels(peak.sub$V8))
        if(length(gene.peak) > 1){
            gene.scores <- data[rownames(data) %in% gene.peak,]
            summed.score <- log1p((Matrix::colSums(gene.scores)/cellave)*10000)
            summed.score <- summed.score - mean(summed.score, na.rm=T)
            summed.score <- rescale(summed.score, c(0,1))
        }else{
            summed.score <- log1p((data[rownames(data) %in% gene.peak,]/cellave)*10000)
            summed.score <- summed.score - mean(summed.score, na.rm=T)
            summed.score <- rescale(summed.score, c(0,1))
        }
        ordered.score <- summed.score
        mat[,cnts] <- ordered.score
    }
    
    # return new reduced dimension matrix
    df <- cbind(rd, mat)
    return(df)
}
normalize.activity <- function(df, acts, n.random=NULL, output="output", logTransform=F, scaleP=F){
    
    # verbose
    message(" - normalizing libraries ...")
    
    # scale by ploidy
    ##the scale ploidy is off
    if(scaleP==T){
        message(" - scaling cell activities by relative ploidy ...")
        nonzero <- as.data.frame(summary(acts))
        nonzero$i <- rownames(acts)[nonzero$i]
        nonzero$j <- colnames(acts)[nonzero$j]
        ploidy <- aggregate(nonzero$x~nonzero$j, FUN=function(x) as.numeric(quantile(x, c(0.5))))
        colnames(ploidy) <- c("cellID", "aveActivity")
        rownames(ploidy) <- ploidy$cellID
        ploidy <- ploidy[colnames(acts),]
        acts <- acts %*% Diagonal(x=1/ploidy$aveActivity)
        colnames(acts) <- rownames(df)
        message(" - num genes: ",nrow(acts), " | num cells: ", ncol(acts))
    }
    
    # quick clean
    acts <- acts[Matrix::rowSums(acts>0)>0,]
    acts <- acts[,Matrix::colSums(acts>0)>0]
    
    # check data frames/matrices
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,colnames(acts) %in% rownames(df)]
    df <- df[colnames(acts),]
    acts.o <- acts
    
    # if select random
    if(!is.null(n.random)){
        acts <- acts[sample(nrow(acts), n.random),]
    }
    
    # find cluster means
    its <- 0
    message(" - extracting average gene activity across cells per cluster ...")
    ##the apply usually gives us information about analyzing the acts for the cluster from 1 to 18
    ##the 1 in the apply suggests the margin 1 for row, it means we set a function for each TE
    ##since the acts row name is TE and colname is cluster
    ##x indicates the rows of acts
    clust.means <- apply(acts, 1, function(x, clust){
        its <<- its + 1 ##it seems no use
        message(paste0('analyzed row is ',its))
        gene.cell <- as.numeric(x)
        ##split the acts by diff clust
        gene.cell.clust <- split(gene.cell, factor(clust, levels=unique(clust)))
        ##gene.cell.clust is for each name Mutator_TIR_transposon_4635$clustercluster_14 and below is the value for the cells under its cluster
        ##lapply() function is useful for performing operations on list objects and returns a list object of same length of original set.
        clust.aves <- as.numeric(lapply(gene.cell.clust, function(y){mean(y, na.rm=T)}))
        ##if we add the as.numeric and the format will be change to matrix
        #clust.aves <- lapply(gene.cell.clust, function(y){mean(y, na.rm=T)})
        ##add the row name
        names(clust.aves) <- names(gene.cell.clust)
        clust.aves
    }, clust=df$Cluster) ##we provide the clust information here
    dim(clust.means) ##18 1078 we have 18 cluster calculate the average acc for each cluster for all genes
    head(clust.means[1:3,1:3])
    
    # plot raw cluster means
    clust.means <- as.matrix(t(clust.means))
    c.means <- clust.means[,mixedorder(colnames(clust.means))]
    row.o <- apply(c.means, 1, which.max)
    class(row.o) ##integer
    c.means <- c.means[order(row.o, decreasing=F),]
    class(c.means) ##change to a matrix
    row.o <- rownames(c.means)
    ##generate a raw heatmap without normalization
    ##we do not generate it since it will cost a lot of time
    #pdf(paste0(output,".raw_heatmap.pdf"), width=5, height=6)
    #heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
    #          useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
    #dev.off()
    ##this file stores the row cluster means for each gene
    write.table(c.means, file=paste0(output,".RawClusterMeans.txt"), 
                quote=F, row.names=T, col.names=T, sep="\t")
    
    
    
    #############
    ##have a test fo the Mclust
    #mod1 <- Mclust(iris[,1:6])
    #summary(mod1)
    #mod1$classification
    #top <- names(mod1$parameters$mean)[which.max(mod1$parameters$mean)]
    #top
    
    #mod1$parameters$mean[which.max(mod1$parameters$mean)]
    #quantile(cars$mpg)
    #mod1$bic
    
    
    ##use mixture model to set a threshold
    # fit 2-guassian mixture model, return higher threshold
    message(" - fitting mixture model ...")
    ##analyze for each column instead of row
    ##row is the cluster name information
    ##what's the meaning of the mod$parameters$mean
    ##since relative gene accessiblity scores exhibited a bimodal distrbution 
    thresholds <- apply(clust.means, 2, function(x){
        mod <- Mclust(x, G=2, verbose=F) ##Mclust is the Gaussian finit mixture models fitted via EM algorithm for model-based clustering classification, and density estimation, including Bayesian regularization, dimension reduction for visualisation, and resampling-based inference.
        ?Mclust
        #specify G:2, An integer vector specifying the numbers of mixture components (clusters) for which the BIC is to be calculated. The default is G=1:9.
        ##BIC values are an approximation to intergrated (not maximum) likelihood, and you want the model with the greatest integrated likelihood (Bayes factor) so you choose the model with the largest BIC.
        top <- names(mod$parameters$mean)[which.max(mod$parameters$mean)]
        upper.cells <- x[which(mod$classification==top)]
        val <- quantile(upper.cells[upper.cells>0], c(0.05))
        #val <- mean(upper.cells, na.rm=T)
        if(val == 0){
            val <- min(upper.cells[upper.cells>0])
        }
        val
        #max(mod$parameters$mean)
    })
    
    # scale cell activity by cluster-specific threshold for expression
    message(" - scaling by cluster averages ...")
    ##the thresholds are values for each cluster: eg. 
    ##clustercluster_1  clustercluster_5 clustercluster_18  clustercluster_3 clustercluster_14 clustercluster_16  clustercluster_9 
    ##0.001389722       0.001471488    0.001372896       0.001923661       0.001318943
    ##each cell has its own cluster information
    ##assign the each cell with the cluster threshold information
    adj.thresh <- thresholds[df$Cluster]
    print(head(adj.thresh, n=10))
    ##dim(Diagonal(x=1/adj.thresh)) is 4655 4655 this means we generate a diagnoal based on the adj.thresh to multiply the orignal act
    adj.act <- acts.o %*% Diagonal(x=1/adj.thresh)
    adj.act@x[is.infinite(adj.act@x)] <- 0
    adj.act@x[is.na(adj.act@x)] <- 0
    adj.act@x <- round(adj.act@x, digits=0) ##take round for the adj.act
    adj.act <- adj.act[Matrix::rowSums(adj.act>0)>0,]
    adj.act <- adj.act[,Matrix::colSums(adj.act>0)>0]
    
    
    ##write out the un-normalized counts
    ##so the above scripts are to scale the whole dataset based on the calculated adj.thresh
    ##but why???????
    # un-normalized counts
    ua.out <- as.data.frame(summary(adj.act))
    ua.out$i <- rownames(adj.act)[ua.out$i]
    ua.out$j <- colnames(adj.act)[ua.out$j]
    ua.out <- ua.out[ua.out$x>0,]
    write.table(ua.out, file=paste0(output,".countsActivity.sparse"), 
                quote=F, row.names=F, col.names=F, sep="\t")
    
    # normalize by size factors
    if(logTransform==F){
        message(" - lib size before size factors = ")
        print(head(Matrix::colSums(adj.act)))
    }
    
    ##this step is to do the normalization
    ## normalizeSF using the 
    # do normalization mean-geometric-mean-total
    message(" - estimating normalization factors ...")
    results <- normalizeSF(adj.act, verbose=T)
    norm.act <- results$norm.act
    norm.factor <- results$norm.factor ##this norm.factor is the output of the sf normalization
    
    ##log transform is FALSE for the default
    # log transform?
    if(logTransform==T){
        message(" - square-root transforming counts activity ...")
        norm.act <- sqrt(norm.act)
        norm.act@x[is.na(norm.act@x)] <- 0
        norm.act@x[is.infinite(norm.act@x)] <- 0
        message(" - lib size afer square-root transformation = ")
        print(head(Matrix::colSums(norm.act)))
    }
    
    # write size factors to meta data file
    ##write out the norm.factor
    df$size_factors <- norm.factor[rownames(df)]
    write.table(df, file=paste0(output,".jaccard.cluster.reducedDims.SF.txt"), 
                quote=F, row.names=T, col.names=T, sep="\t")
    
    # verbose
    message(" - lib size after size factors = ")
    print(head(Matrix::colSums(norm.act)))
    
    # remove empty cells/genes
    norm.act <- norm.act[Matrix::rowSums(norm.act>0)>0,]
    norm.act <- norm.act[,Matrix::colSums(norm.act>0)>0]
    message(" - cells = ",ncol(norm.act), " | genes = ", nrow(norm.act))
    print(head(norm.act[,1:10]))
    
    # print output to disk
    ##write out results to sparse
    ia.out <- as.data.frame(summary(norm.act))
    ia.out$i <- rownames(adj.act)[ia.out$i]
    ia.out$j <- colnames(adj.act)[ia.out$j]
    ia.out <- ia.out[ia.out$x>0,]
    write.table(ia.out, file=paste0(output,".normalizedActivity.sparse"), quote=F, row.names=F,
                col.names=F, sep="\t")
    
    # return
    message(" - returning normalized matrix ...")
    return(list(norm.act=norm.act, norm.factor=norm.factor, adj.act=adj.act, row.o=row.o))
}
smooth.data        <- function(x, k=25, step=3, npcs=30, cleanExp=F, df=NULL, rds=NULL, output="output"){
    
    # verbose
    message(" - imputing gene activity ...")
    
    # input
    data.use <- x
    
    # verbose
    if(!is.null(rds)){
        
        if(!is.null(df)){
            message("   * using UMAP manifold for smoothing ...")
            pcs <- df[,c("umap1","umap2")]
        }else{
            message("   * using prior PC space as manifold ...")
            pcs <- rds[colnames(x),c(1:npcs)]
        }
    }else{
        
        # LSI
        message("   * PC manifold set to NULL, running LSI (TFIDF)...")
        x[x>0] <- 1
        tf.idf <- tfidf(x)
        
        # get PCS
        message("   * PC manifold set to NULL, running LSI ...")
        pc <- irlba(t(tf.idf), npcs)
        pcs <- pc$u 
        rownames(pcs) <- colnames(x)
        colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
        
        # do l2-norm
        pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
    }
    
    # get KNN
    message("   * finding knn graph ...")
    knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
    j <- as.numeric(x = t(x = knn.graph))
    i <- ((1:length(x = j)) - 1) %/% k + 1
    edgeList = data.frame(i, j, 1)
    A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
    
    # Smooth graph
    message("   * smoothing graph ...")
    A = A + t(A)
    A = A / Matrix::rowSums(A)
    step.size = step
    if(step.size > 1){
        for(i in 1:step.size){
            message("     ~ step ",i)
            A = A %*% A
        }
    }
    
    # smooth data
    message("   * smoothing activity ...")
    out <- t(A %*% t(data.use))
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    print(head(out[,1:10]))
    
    # find non-expressed genes
    if(cleanExp==T){
        message("   * finding activity thresholds ...")
        its <- 0
        new <- sparse_apply(out, 2, function(z){
            its <<- its + 1
            mod <- Mclust(z, G=2)
            top.parameter <- which.max(mod$parameters$mean)
            vals <- split(z, mod$classification)[[top.parameter]]
            threshold <- quantile(vals[vals>0], c(0.05))[1]
            if(its %% 100 == 0){message("   * iterated over ",its, " cells (min thresh = ",threshold,") ...")}
            z[z < threshold] <- 0
            return(z)
        }, convert_to_dense=F)
        impute.activity <- Matrix(new, sparse=T)
    }else{
        impute.activity <- out
    }
    
    # clean empty rows (if present) and round to two decimal places
    print(head(impute.activity[,1:10]))
    message("   * clean after imputation ...")
    impute.activity <- impute.activity[Matrix::rowSums(impute.activity)>0,]
    impute.activity <- impute.activity[,Matrix::colSums(impute.activity)>0]
    impute.activity <- impute.activity %*% Diagonal(x=1e6/Matrix::colSums(impute.activity))
    impute.activity@x <- round(impute.activity@x, digits=2)
    
    # write to disk
    # ia.out <- as.data.frame(summary(impute.activity))
    # ia.out$i <- rownames(impute.activity)[ia.out$i]
    # ia.out$j <- colnames(impute.activity)[ia.out$j]
    # write.table(ia.out, file=paste0(output,".imputedActivity.sparse"), quote=F, row.names=F,
    #             col.names=F, sep="\t")
    
    # return sparse Matrix
    return(impute.activity)
}
clusterAves        <- function(df, acts, output, markers){
    
    # verbose start-up
    if(is.null(df) | is.null(acts)){
        stop("ERROR: must provide metadata and activity matrix")
    }
    if(! "Cluster" %in% colnames(df)){
        df$Cluster <- df$umapclust
    }
    rownames(markers) <- markers$geneID
    
    # estimate significant differences
    df <- df[colnames(acts),]
    clust.o <- mixedsort(unique(df$Cluster))
    df$Cluster <- factor(df$Cluster, levels=clust.o)
    
    # iterate over genes
    message(" - transversing gene activity ...")
    it <- 0
    df.acts <- as.data.frame(t(as.matrix(acts)))
    df.splits <- split(df.acts, df$Cluster)
    difs.mean <- lapply(df.splits, function(z){colMeans(z)})
    difs.mean <- as.matrix(do.call(cbind, difs.mean))
    
    # plot
    c.means <- difs.mean[,mixedorder(colnames(difs.mean))]
    c.means <- c.means[rownames(c.means) %in% rownames(markers),]
    markers <- markers[rownames(c.means),]
    row.o <- apply(c.means, 1, which.max)
    markers <-  markers[order(row.o, decreasing=F),]
    c.means <- c.means[order(row.o, decreasing=F),]
    
    # rowSide cols
    r.cols <- colorRampPalette(brewer.pal(8,"Set2"))(length(unique(markers$type)))
    
    # plot
    pdf(paste0(output,"_markers_heatmap.pdf"), width=10, height=12)
    heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
              useRaster=F, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), 
              margins=c(10,10), cexRow = 0.7, RowSideColors = r.cols[factor(markers$type)],
              labRow=paste(markers$name,markers$type,sep="-"))
    dev.off()
    
    # estimate shannon
    sh <- apply(c.means, 1, shannon.entropy)
    total <- sum(sh)
    print(sh)
}

# plotting
plot.scores        <- function(df=df, genes=NULL, info=info){
    
    # plot all genes
    if(is.null(genes)){
        gids <- colnames(df)[grepl("Zm00001d",colnames(df))]
        info.ordered <- info[gids,]
        nrows <- ceiling(length(gids)/4)
        totals <- nrows*4
        ratio <- nrows/4
        
        # params
        pdf(file="markerAccessibility.pdf", width=16, height=ratio*16)
        layout(matrix(c(1:totals), ncol=4, byrow=T))
        par(mar=c(2,2,1,1))
        
        # iterate over each gene
        for(i in 1:length(gids)){
            geneID <- gids[i]
            acv <- df[,geneID]
            cols <- colorRampPalette(c("grey75", magma(10)))(11)
            colvec <- cols[cut(acv, breaks=12)]
            sizes <- rescale(acv, c(0.2, 0.5))
            
            # plot
            plot(df$umap1, df$umap2, col=colvec, 
                 main=paste(gids[i],info.ordered$type[i],sep="-"),
                 xlab="", ylab="", bty="n",
                 xaxt="n", yaxt="n", pch=16, cex=sizes)
        }
        
        # turn device off
        dev.off()
        
    }else{
        # plot subset of genes
        stop("gene subsets not yet supported")
    }
}
plot.act.scores    <- function(df, acts=acts, info=NULL, top=NULL, logT=F, marker.dist=NULL,
                               outname="markerActivityScores.pdf", lim=0.95){
    
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
    pdf(file=outname, width=16, height=ratio*16)
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
        cols <- colorRampPalette(c("grey80","grey80","grey75",brewer.pal(9, "BuPu")[2:9]), bias=1)(100)
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
            min.acv <- -0.1
            max.acv <- max(acv) + (0.05*max(acv))
            message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
            colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
        }
        colvec[is.na(colvec)] <- cols[1]
        sizes <- rescale(acv, c(0.25, 0.3))
        
        # plot
        plot(df2$umap1, df2$umap2, col=colvec, 
             main=paste(info$name[i],info$type[i],sep="-"),
             xlab="", ylab="", bty="n",
             xaxt="n", yaxt="n", pch=16, cex=0.15)
      
    }
    
    # turn device off
    dev.off()
    
}

#############
##have a test to decipher the DEGs in each cell
#plot.new.markers(df=b, 
#                 acts=activity, 
#                 top=5, 
#                 normT='mean.dif',
#                 logT=T,
#                 lim=0.99,
#                 outname=paste0(output,".normalized"))
#function(b, activity, impute.activity, output="all", marker.info)
#runMajorDeNovo(out$b, out$activity, out$impute.activity, marker.info.dat, output="all")

#out <- readRDS('out.rds')
#dat <- readRDS('dat.rds')
#df <-  out$b
#acts <- out$activity
#marker.info.dat <- marker.info.dat <- dat$marker.info



plot.new.markers   <- function(df=NULL, acts=NULL, outname="ActivityScores.pdf", row.o=NULL,
                               top=5, normT='prop.dif', logT=F, threshold=0.01, lim=0.95){
    
    # verbose start-up
    if(is.null(df) | is.null(acts)){
        stop("ERROR: must provide metadata and activity matrix")
    }
    if(! "Cluster" %in% colnames(df)){
        df$Cluster <- df$umapclust
    }
    
    # estimate significant differences
    df <- df[colnames(acts),]
    df$Cluster <- factor(df$Cluster)
    
    # iterate over genes
    message(" - transversing gene activity ...")
    ##rowname is the gene name and colname is the cell name
    it <- 0
    total.mean <- Matrix::rowMeans(acts) ##generate the mean of acts for each te from all the cells
    props <- acts
    props@x <- rep(1, length(props@x)) ##generate props with 1 to replace the acts
    total.prop <- Matrix::rowMeans(props) ##generate mean of '1' act for each te from all the cells
    df.splits <- split(colnames(acts), df$Cluster) ##split the name of cells to different cluster
    
    # estimate mean, proportions
    ##we will choose what types of analyses we want to use, the proportion which means that we will use the 1 instead of exact act number to indicate
    message(" - estimating difference in proportions and means ...")
    ##use the lapply to have a function to analyze each splitted cells for each cluster
    ##calculate the average of acts from cells of one cluster for each TE
    ##if this is for the proportion that means for each cluster, what the proportion of cells are expressed
    difs.prop <- lapply(df.splits, function(z){
        sub.act <- acts[,z] ##sub.act indicates the act information for each cluster
        Matrix::rowMeans(sub.act >0)
    })
    length(difs.prop)
    
    difs.mean <- lapply(df.splits, function(z){
        sub.act <- acts[,z]
        rowMeans(sub.act)
    })
    
    length(difs.mean)
    
    difs.prop <- do.call(cbind, difs.prop) ##combine the each cluster information to one matrix
    difs.mean <- do.call(cbind, difs.mean)
    colnames(difs.prop) <- names(df.splits)
    colnames(difs.mean) <- names(df.splits)
    
    # estimate pairwise log2 fold change
    message(" - estimating mean log2 fold changes")
    ##this apply function means we will analyze on each row or say each TE
    ave.difs.mean <- t(apply(difs.mean, 1, function(x){
        sapply(x, function(z){mean(log2((z+1)/(x+1)), na.rm=T)})
    }))
    ##sapply(X, FUN)
    #Arguments:
    #  -X: A vector or an object
    #-FUN: Function applied to each element of x
    #it means it will compare value from one to another one cluster for all the combinations and next take the mean of the value of the log2 fold changes.
    #z means each element which is the exact value
    #x means a vector of vaules that contains values from each cluster
    #the positive value is higher expressed 
    
    
    # save output
    ##all the average log2 fold change
    write.table(ave.difs.mean, file=paste0(outname,".NormalizedAveLog2FC.txt"), 
                quote=F, row.names=T, col.names=T, sep="\t")
    write.table(difs.mean, file=paste0(outname,".NormalizedClusterMeans.txt"), 
                quote=F, row.names=T, col.names=T, sep="\t")
    write.table(difs.prop, file=paste0(outname,".NormalizedClusterProportions.txt"), 
                quote=F, row.names=T, col.names=T, sep="\t")
    
    # output options
    if(normT=="mean.dif"){
        difs <- ave.difs.mean
    }else if(normT=="prop.dif"){
        difs <- difs.prop
    }else if(normT=="adj.dif"){
        difs <- ave.difs.mean * difs.prop
    }
    ##since we set the normT as the difs
    ##we will use the prop.dif
    ##the mixedorder means we order the colnames of the difs
    difs <- difs[,mixedorder(colnames(difs))]
    
    # # plot normalized heatmap
    # c.means <- difs.mean[,mixedorder(colnames(difs.mean))]
    # pdf(paste0(output,".normalized_heatmap.pdf"), width=5, height=6)
    # heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
    #           useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
    # dev.off()
    
    # filter out genes by proportion in at least one cluster
    message(" - selecting genes passing minimum proportion threshold ...")
    ##the threshold we set 0.01
    ##for the each row 
    #threshold <- 0.01
    ##we will analyze on each row
    ##if the each value of each row is higher than the 0.01, we will count 1, we count how many value are over the threshold. the total number will of 18 since 18 clusters
    ##be careful we use the difs.prop to check the pass.genes
    ##the passed genes will be used on the difs
    pass.genes <- apply(difs.prop, 1, function(x){length(x[x>threshold])})
    difs <- difs[pass.genes > 0,]
    
    # melt
    reduced <- melt(difs)
    reduced <- reduced[order(reduced$value, decreasing=T),]
    ##select the top 5 tes for each cluster
    top.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=top))
    
    ##we need to write out all the top
    all_top.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=1000000000))
    ##write out the all_top.genes information for each cluster
    write.csv(all_top.genes,paste0(outname,"opt_rank_difs.csv"))
    
    ##only left the value is over 1 which means this te has value in all cells of that cluster
    ##default is prop but actually we use the mean so the value will be over 1
    ##if we set the prop the max value will be 1
    ##the reason why we set > 1 since we log the proportion log2(1.5/1) if 1.5 > 1 it means the value will be over 1
    top.genes.out <- subset(reduced, reduced$value >= 1)
    write.table(top.genes.out, file=paste0(outname,".log2FC1genes.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    topN <- top.genes$Var1
    ##topN means we select the top 5 for each cluster since we have 18 clusters we have 90 TEs
    message(" - # genes = ", nrow(top.genes), " | # per cluster: ", top)
    print(head(top.genes))
    
    ##so the final genes pass the threshold is in the all.normalized.log2FC1genes.txt file
    
    
    # plot top 50 -----------------------------------------------------------------------------
    nrows <- ceiling(length(topN)/top)
    totals <- nrows*top
    ratio <- nrows/top
    
    # params
    pdf(file=paste0(outname,".denovo.pdf"), width=20, height=ratio*20)
    layout(matrix(c(1:totals), ncol=top, byrow=T))
    par(mar=c(2,2,1,1))
    
    # adjust cluster IDs
    message("begin plotting ...")
    #i <- 1
    for (i in 1:nrow(top.genes)){
        
        ##we plot based on the values frome the acts form
        ##we do not have log
      
        # copy meta data
        dff <- df
        geneID <- top.genes$Var1[i]
        ##extract acv for the specific gene
        acv <- as.numeric(acts[rownames(acts) %in% geneID,])
        if(logT==T){
            acv <- log2(acv+1)
        }
        message(" - gene: ", geneID)
        
        # set up plot cols/sizes
        orderRow <- order(acv, decreasing=F)
        #cols <- colorRampPalette(c("grey85",brewer.pal(9,"YlGnBu")), bias=1)(100)
        #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
        #cols <- inferno(100)
        #cols <- plasma(100)
        ##generate 100 colors 
        cols <- colorRampPalette(viridis(100), bias=0.5)(100)
        #cols <- colorRampPalette(c("grey75","darkorange","firebrick3"), bias=0.75)(100)
        acv <- acv[orderRow]
        df2 <- df[orderRow,]
        
        # set up upper limits
        ##generate 95% upper limits 15
        upper.lim <- quantile(acv, lim)
        ##if the acv in one cell > 15, we allow the acv to be 15
        acv[acv > upper.lim] <- upper.lim
        
        min.acv <- -0.1
        max.acv <- max(acv) + (0.05*max(acv))
        ##we generate colvec which means we divide the gap between min.acv and max.acv to have 101 colors 
        colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv,length.out=101))]
        colvec[is.na(colvec)] <- cols[1] ##set na is to cols[1]
        sizes <- rescale(acv, c(0.35, 0.4)) ##rescale the sizes to 0.35 to 4
        
        # plot
        plot(df2$umap1, df2$umap2, 
             col=colvec, 
             main=paste(top.genes$Var1[i],top.genes$Var2[i],sep="-"),
             xlab="", ylab="", bty="n", xaxt="n", yaxt="n", pch=16, cex=0.25)
        add.color.bar(6, colvec, title=NULL, lims=c(0,1),prompt=F,x=min(df2$umap1),y=min(df2$umap2)+1)
        
    }
    
    # turn device off
    dev.off()
    
}
plot.dists         <- function(impute.activity, out="normalized_geneactivity_distribution.pdf"){
    
    # output
    pdf(out, width=6, height=5)
    
    # first plot
    den1 <- density(impute.activity[,1], from=0,to=max(impute.activity))
    den1$y <- rescale(den1$y, c(0,1))
    plot(den1, col=alpha("grey50", 0.5), lwd=0.75)
    
    # all other plots
    for(i in 2:ncol(impute.activity)){
        den <- density(impute.activity[,i], from=0,to=max(impute.activity))
        den$y <- rescale(den$y, c(0,1))
        lines(den$x, den$y, col=alpha("grey50", 0.5), lwd=0.75)
    }
    
    # write to disk
    dev.off()
    
}

########
##no use
findTopMakers      <- function(acts, df, sub.list=NULL, group=NULL){
    # verbose start-up
    if(is.null(df) | is.null(acts)){
        stop("ERROR: must provide metadata and activity matrix")
    }
    if(! "Cluster" %in% colnames(df)){
        df$Cluster <- df$umapclust
    }
    if(!is.null(sub.list)){
        acts <- acts[rownames(acts) %in% sub.list,]
    }
    
    # estimate significant differences
    acts <- acts[,rownames(df)]
    clust.o <- mixedsort(unique(df$Cluster))
    df$Cluster <- factor(df$Cluster, levels=clust.o)
    
    # iterate over genes
    message(" - transversing gene activity ...")
    it <- 0
    difs <- t(apply(acts, 1, function(x, y){
        it <<- it + 1
        clust.gene <- split(x, y)
        iit <- 0
        m.dif <- sapply(clust.gene, FUN=function(z, all){
            iit <- iit + 1
            clust.exp <- as.numeric(z)
            nonclust.exp <- as.numeric(all[!names(all) %in% names(clust.exp)])
            p1 <- length(clust.exp[clust.exp > 0]) / length(clust.exp)
            p2 <- length(nonclust.exp[nonclust.exp > 0]) / length(nonclust.exp)
            prop.dif <- log2((p1+0.0001)/(p2+0.0001))
            mean.dif <- log2( (mean(clust.exp, na.rm=T)+0.01) / 
                                  (mean(nonclust.exp, na.rm=T)+0.01) )
            adj.dif <- prop.dif * abs(mean.dif)
            
            # return data type
            if(normT=='adj.dif'){
                return(adj.dif)
            }else if(normT=='mean.dif'){
                return(mean.dif)
            }else if(normT=='prop.dif'){
                return(prop.dif)
            }else{
                return(adj.dif)
            }
        }, all=x)
        if(it %% 1000 == 0){message("   * iterated over ", it, " genes ...")}
        #if(it %% 1000 == 0){message("     ~ ", head(m.dif)[1]," ",head(m.dif)[2],
        #                            " ",head(m.dif)[3]," ",head(m.dif)[4]," ",head(m.dif)[5])}
        return(m.dif)
        
    }, y=df$Cluster))
    colnames(difs) <- levels(df$Cluster)
    
    # melt
    reduced <- melt(difs)
    reduced <- reduced[order(reduced$value, decreasing=T),]
    top.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=top))
    topN <- top.genes$Var1
    message(" - # genes = ", nrow(top.genes), " | # per cluster: ", top)
    print(head(top.genes))
}

#################
##updation 112020
##have a test for the script
#args    
#meta <- 'AtRoot.meta.v5.celltypes.txt'
#geneact <- 'cicero_gene_activities_sparse_mtx.rds'
#pcs <- 'Arabidopsis.root.SVD.txt'
#mark <- 'TE.bed'
#threads <- '1'

#dat <- loadData(meta, pcs, geneact, mark)    
#b.meta <- dat$b
#activity.all <- dat$activity
#h.pcs1 <- dat$h.pcs
#marker.info.dat <- dat$marker.info

#all.b <- b.meta
#all.activity <- activity.all
#all.hpcs <- h.pcs1
#marker.info <- marker.info.dat

#df <- all.b
#acts <- all.activity


# parallel implementation over major clusters
runMajorPriori     <- function(all.b, all.activity, all.hpcs, marker.info, threads=1, output="all",
                               smooth.markers=T, i = 3, kn = 20){
    
    ##intersect with the meta data
    # preprocess
    ids <- intersect(rownames(all.b), colnames(all.activity))
    ids <- intersect(ids, rownames(all.hpcs))
    all.b <- all.b[ids,]
    all.activity <- all.activity[,ids]
    all.hpcs <- all.hpcs[ids,]
    
    # subset by major clusters
    all.b$LouvainClusters <- as.character(all.b$LouvainClusters)
    all.activity <- all.activity[,rownames(all.b)]
    ##do filtration that allows the rowsum of activity to be over 0
    ##and make the final all.hpcs and all.b
    all.activity <- all.activity[Matrix::rowSums(all.activity)>0,]
    all.activity <- all.activity[,Matrix::colSums(all.activity)>0]
    all.hpcs <- all.hpcs[colnames(all.activity),]
    all.b <- all.b[colnames(all.activity),]

    # normalize per cell activity by cluster average and size factors
    results <- normalize.activity(all.b, all.activity, output=output, logTransform=F, scaleP=F)
    activity <- results$norm.act
    row.o <- results$row.o
    rm(results)
    gc()
    
    ##smooth the gene 
    # if smooth genes only
    if(smooth.markers){
        activity <- activity[rownames(activity) %in% as.character(marker.info$geneID),]
    }
    
    # impute gene accessibility scores
    impute.activity <- smooth.data(activity, k=kn, step=i, npcs=ncol(all.hpcs), df=NULL,
                                       rds=all.hpcs, cleanExp=F, output=output)
    ia <- as.data.frame(summary(impute.activity))
    ia$i <- rownames(impute.activity)[as.numeric(ia$i)]
    ia$j <- colnames(impute.activity)[as.numeric(ia$j)]
    write.table(ia, file=paste0(output,".smoothed.sparse"), quote=F, row.names=F, col.names=F, sep="\t")
    
    # ranges
    marker.impact <- impute.activity[rownames(impute.activity) %in% as.character(marker.info$geneID),]
    marker.ranges <- lapply(rownames(marker.impact), function(x){
        x[is.na(x)] <- 0
        min.x <- min(marker.impact[x,], na.rm=T)
        max.x <- max(marker.impact[x,], na.rm=T)
        rng <- c(min.x,max.x)
        rng[is.na(rng)] <- 0
        rng[1] <- rng[1] - 0.1
        rng[2] <- rng[2] + (0.05*rng[2])
        return(rng)
    })
    names(marker.ranges) <- rownames(marker.impact)
    
    # plot all
    ##plot all the marker information
    ##so we have 1400 TEs we will plot all the 1400 TEs together
    plot.act.scores(all.b, acts=activity, 
                    info=marker.info, 
                    logT=T,
                    lim=0.99,
                    marker.dist=NULL,
                    outname=paste0(output,"combined.",".normalized.known.Markers.pdf"))
    
    plot.act.scores(all.b, acts=impute.activity, 
                    info=marker.info, 
                    logT=F,
                    lim=0.995,
                    marker.dist=marker.ranges,
                    outname=paste0(output,"combined.",".impute.known.Markers.pdf"))
    
    write.csv(all.b,'afterrunMajorPriori_meta.csv')
    #write.csv(impute.activity,'impute.activity.csv')
    # return
    return(list(b=all.b, activity=activity, impute.activity=impute.activity))
    
}
runMajorDeNovo     <- function(b, activity, impute.activity, output="all", marker.info){
    
    # verbose
    message(" - finding de novo enriched markers ...")
    
    # find de novo markers
    plot.new.markers(df=b, 
                     acts=activity, 
                     top=5, 
                     normT='mean.dif',
                     logT=T,
                     lim=0.99,
                     outname=paste0(output,".normalized"))
    
    plot.new.markers(df=b, 
                     acts=impute.activity, 
                     top=5, 
                     normT='mean.dif', ##we use the mean.dif
                     logT=F,
                     lim=0.99,
                     outname=paste0(output,".imputed"))
    
    ##modify from hd do not generate these two
    # cluster aves
    #clusterAves(b, activity, paste0(output,".normalized"), marker.info)
    #clusterAves(b, impute.activity, paste0(output,".imputed"), marker.info)

}


