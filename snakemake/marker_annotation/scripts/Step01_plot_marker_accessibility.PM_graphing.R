###################################################################################################
## plot marker accessibility scores
###################################################################################################
##updating 062221 add the target cluster argument
##updating 062221 create a fold to store each type of markers in a seperate pdf that is suitable for the many markers within one cell type

# load arguments
args <- commandArgs(trailingOnly=T)
#if(length(args) != 5){stop("Rscript plot_marker_accessibility.R [meta] [gene_activity] [pcs.txt] [markers.bed] [threads]")}

 .libPaths("/home/jpm73279/.conda/envs/R_final_install/lib/R/library")
#args    
meta <- as.character(args[1])
geneact <- as.character(args[2])
pcs <- as.character(args[3])
mark <- as.character(args[4])
threads <- as.numeric(args[5])

##updating 062221
target_cluster <- as.character(args[6])
plot_each_CT <- as.character(args[7]) ##use yes to initiate this argument
output_dir <- as.character(args[8])
function_script <- as.character(args[9])

knn <- as.numeric(args[10])
smooth_val <- as.numeric(args[11])
output_base_name <- as.character(args[12])

# load functions
#source("functions.plot_marker_accessibility.R")
source(function_script)


# load data
#if(file.exists(paste0(output_dir,"/GAobj.rds"))){
#    dat <- readRDS(paste0(output_dir,"/GAobj.rds"))
#}else{
#    ##updating 062221 save to the GAobj
#    saveRDS(dat,paste0(output_dir,"/GAobj.rds"))
#}

print(target_cluster)
dat <- loadData(meta, pcs, geneact, mark, target_cluster)
b.meta <- dat$b
activity.all <- dat$activity
h.pcs1 <- dat$h.pcs
marker.info.dat <- dat$marker.info

##udpating 062221 create a marker opt dir
##updation add the checking of dir
marker_opt_dir <- paste0(output_dir,'/opt_store_eachCT_marker_dir')
if (!dir.exists(marker_opt_dir)){
  dir.create(marker_opt_dir)
} else {
  print("Dir already exists!")
}


# iterate over each major cluster
out <- runMajorPriori(b.meta, activity.all, h.pcs1, marker.info.dat, plot_each_CT, 
                      marker_opt_dir,output_dir,threads=threads, 
                      smooth.markers=F, knn_val = knn, smooth = smooth_val, output = output_base_name)

gen_output_name <- paste0(output_base_name, ".out.rds")
saveRDS(out,paste0(output_dir,gen_output_name))

#de_novo <- runMajorDeNovo(output_dir,out$b, out$activity, out$impute.activity, marker.info.dat,output=output_base_name)
#gen_output_name_novo <- paste0(output_base_name, ".out.de_novo.rds")
#saveRDS(de_novo,paste0(output_dir,gen_output_name_novo))


