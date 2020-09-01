library(graphics)
source("scripts/outbred.pruned.R")
source("scripts/emma.R")

#get list of samples
phenotypes = read.table("/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt" ,header=T, sep="\t")
samples = phenotypes$"Animal_ID"

#make kinship.matrices object
#only need to do once; afterward can load Rdata file 
load("/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/grms/kinship_matrices.Rdata")
#print("making kinship matrix rdata files...")
#K.chr <- vector("list", 19)
#names(K.chr) <- 1:19
k_dir="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/grms/"
g_dir="/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/"
for( i in 1:19){
  print(i)
  load(paste0(k_dir, "combined.chr", toString(i), ".oxinfo90.kinship.RData", sep=""))
  a = read.table(paste0(g_dir, "combined.chr", toString(i), ".oxinfo90.annot", sep=""), header=FALSE)
  numSnps = nrow(a)
  colnames(k) <- samples
  rownames(k) <- samples
  K.chr[[i]]$K = k
  K.chr[[i]]$n = numSnps
}
kinship.matrices <- vector("list", 1)
kinship.matrices[["K.chr"]] = K.chr

save(kinship.matrices, file="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/grms/kinship_matrices.Rdata")


