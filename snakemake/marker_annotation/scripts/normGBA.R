# normalize raw gene body accessibility ##

# load arguments
args <- commandArgs(T)
if(length(args)!=5){stop("Rscript normGBA.R <gene.sparse> <meta> <Zea_mays.AGPv4.36.Allgene.nuclear.bed> <prefix> <F>")}
input <- as.character(args[1])
meta <- as.character(args[2])
gene <- as.character(args[3])
prefix <- as.character(args[4])
earlyExit <- as.logical(args[5])

# load libraries
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(sctransform)
library(ggplot2)

# load functions
rowStdSparse <- function(x){
    r.var <- sqrt(Matrix::rowSums((Matrix::rowMeans(x)-x)^2) / ncol(x))
    Diagonal(x=1/r.var) %*% (x-Matrix::rowMeans(x))
}

# load data
message(" - loading input data ...")
a <- read.table(input)

message(" - loading meta data ...")
b <- read.table(meta)
rownames(b) <- b$cellID

message(" - loading Annotation data ...")
g <- read.table(gene)

# process
a <- sparseMatrix(i=as.numeric(factor(a$V1)),
                  j=as.numeric(factor(a$V2)),
                  x=as.numeric(a$V3),
                  dimnames=list(levels(factor(a$V1)),levels(factor(a$V2))))
a <- a[,colnames(a) %in% rownames(b)]
b <- b[colnames(a),]
rownames(g) <- g$V4
g <- g[rownames(a),]
g$gene.len <- g$V3 - g$V2
g <- subset(g, g$gene.len > 100)
a <- a[rownames(g),]

print("This is g")
print(head(g, n = 10))


print("This is a")
print(head(a), n = 10)



# gene attribute
gene_attr <- data.frame(mean = Matrix::rowMeans(a), 
                        detection_rate = Matrix::rowMeans(a > 0),
                        var = apply(a, 1, var))

print("Gene attr")
print(head(gene_attr, n = 10))

gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(a)
cell_attr <- data.frame(n_umi = Matrix::colSums(a),
                        n_gene = Matrix::colSums(a > 0),
                        log_umi = log10(Matrix::colSums(a)))

print("cell attr")
print(head(cell_attr, n = 10))
rownames(cell_attr) <- colnames(a)

# normalize accessibility by gene length
message(" - normalizing gene accessibility ...")
norm <- sctransform::vst(a, cell_attr = cell_attr,
                         latent_var = c('log_umi'), 
                         return_gene_attr = TRUE, 
                         return_cell_attr = TRUE, 
                         verbosity = 2)
norm <- correct(norm, do_round = T)
norm <- Matrix(norm, sparse=T)

# filter rows/columns - 100 genes accessible per cell & at least 10 cells accessible per gene
nors  <- norm[,Matrix::colSums(norm>0)>0]
norm <- norm[Matrix::rowSums(norm>0)>0,]

# write as counts
cnts <- as.data.frame(summary(norm))
cnts$i <- rownames(norm)[cnts$i]
cnts$j <- colnames(norm)[cnts$j]
write.table(cnts, file=paste0(prefix,".sctGBAcounts.sparse"),
            quote=F, row.names=F, col.names=F, sep="\t")

# normalize to sum gene accessibility to 1
ids <- colnames(norm)
norm <- norm %*% Diagonal(x=1/Matrix::colSums(norm))
colnames(norm) <- ids

# write output
out <- as.data.frame(summary(norm))
out$i <- rownames(norm)[out$i]
out$j <- colnames(norm)[out$j]
write.table(out, file=paste0(prefix,".GBaccessibility.sparse"),
            quote=F, row.names=F, col.names=F, sep="\t")

# exit early
if(earlyExit){
    stop(" - finished - ")
}

# format as z-score matrix
message(" - estimating z-scores ...")
z.sp <- rowStdSparse(norm)
colnames(z.sp) <- colnames(norm)
rownames(z.sp) <- rownames(norm)

# pcs
message(" - running SVD ...")
svd.out <- irlba(t(z.sp), 6)
svd.var <- diag(svd.out$d)
svd.var[1,1] <- 0
LSI_out <-  t(t(svd.var %*% t(svd.out$u)) %*% t(svd.out$v))
colnames(LSI_out) <- colnames(z.sp)
rownames(LSI_out) <- rownames(z.sp)
LSI_out <- t(scale(t(LSI_out)))
quants <- quantile(LSI_out, c(0.01, 0.99))
max.z <- quants[2]
min.z <- quants[1]
message(" - LSI caps: ", min.z, " | ", max.z)
LSI_out[LSI_out > max.z] <- max.z
LSI_out[LSI_out < min.z] <- min.z
print(head(LSI_out[,1:5]))

## cluster cells and peaks
#message(" - clustering cells and peaks ...")
#cell.pcs <- svd.out$u %*% svd.var
#gene.pcs <- svd.out$v %*% svd.var
#hclust_cells <- hclust(proxy::dist(cell.pcs, method="correlation"), method="ward.D2")
#hclust_genes <- hclust(proxy::dist(gene.pcs, method="correlation"), method="ward.D2")
#
## output
#message(" - plotting heatmap ...")
#
## colors
#cols <- colorRampPalette(c("cyan","grey25","gold"))(30)
#pdf(paste0(prefix,".geneAccessiblity.heatmap.pdf"), 
#     width=6, height=7)
#
## plot
#heatmap.2(LSI_out, 
#          trace='none', dendrogram='both', scale='none',
#          col=cols, useRaster=T,
#          Rowv=as.dendrogram(hclust_genes),
#          Colv=as.dendrogram(hclust_cells),
#          labRow=F, labCol=F)
#
## dev.off
#dev.off()
