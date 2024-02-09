###################################
##merge gene body and gene activity

##load libraries
library(Matrix)

args <- commandArgs(T)
geneact <- as.character(args[1]) ##geneActivity
genebod <- as.character(args[2]) ##genebody
output  <- as.character(args[3]) ##set prefix

##load data
act <- read.table(geneact)
bod <- read.table(genebod)

##keep the same cells
shared.cells <- intersect(unique(as.character(act$V2)), unique(as.character(bod$V2)))
bod <- bod[as.character(bod$V2) %in% shared.cells,]
act <- act[as.character(act$V2) %in% shared.cells,]

##merge data
colnames(act) <- c("gene","cell","activity")
colnames(bod) <- c("gene","cell","geneaccess")

message("Merging Cells and Genes")
##set the weight to merge
cboth <- merge(act, bod, by=c("gene","cell"), all=T)
cboth[is.na(cboth)] <- 0
cboth <- as.data.frame(cboth)
cboth$combined <- cboth$activity + (3*cboth$geneaccess)
cboth$activity <- NULL
cboth$geneaccess <- NULL
cboth$gene <- as.factor(cboth$gene)
cboth$cell <- as.factor(cboth$cell)

message("Normalizing Matrix")
##normalize the columns to be 1
combine_mat <- sparseMatrix(i=as.numeric(cboth$gene),
			 j=as.numeric(cboth$cell),
			 x=as.numeric(cboth$combined),
			 dimnames=list(levels(cboth$gene),levels(cboth$cell)))

message("Col Names Command")
mat_id <- colnames(combine_mat)

message("Combined Matrix")
combine_mat <- combine_mat %*% Diagonal(x=1/Matrix::colSums(combine_mat))
colnames(combine_mat) <- mat_id

##transfer to the sparse
cboth <- as.data.frame(summary(combine_mat))
cboth$i <- rownames(combine_mat)[cboth$i]
cboth$j <- colnames(combine_mat)[cboth$j]
write.table(cboth, file=paste0(output,".GAadjusted.sparse"),
            row.names=F, col.names=F, quote=F, sep="\t")

