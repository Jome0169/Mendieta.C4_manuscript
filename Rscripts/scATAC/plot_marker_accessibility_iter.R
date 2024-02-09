############################################
## pipeline plot marker accessibility scores

# load arguments
args <- commandArgs(trailingOnly=T)

#args    
meta <- as.character(args[1])
geneact <- as.character(args[2])
pcs <- as.character(args[3])
mark <- as.character(args[4])
threads <- as.numeric(args[5])
knear <- as.numeric(args[6])
steps <- as.numeric(args[7])
output_dir <- as.numeric(args[8])

# load functions
source("functions.plot_marker_accessibility_iter.R")

# load data
if(file.exists("GAobj.rds")){
    dat <- readRDS("GAobj.rds")
}else{
    dat <- loadData(meta, pcs, geneact, mark)    
}
b.meta <- dat$b
activity.all <- dat$activity
h.pcs1 <- dat$h.pcs
marker.info.dat <- dat$marker.info

# iterate over each major cluster
out <- runMajorPriori(b.meta, activity.all, h.pcs1, marker.info.dat, 
                      threads=threads, smooth.markers=F, 
                      i = steps, kn = knear, output = output_dir )
saveRDS(out, paste0(output_dir,'out.rds'))
#runMajorDeNovo(out$b, out$activity, out$impute.activity, marker.info.dat, output="all")
