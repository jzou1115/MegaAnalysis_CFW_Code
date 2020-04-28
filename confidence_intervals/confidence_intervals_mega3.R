library(parallel)
#library(emma)
library(graphics)
source("/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_Code/confidence_intervals/outbred.pruned.R")
#library("RevoUtilsMath") #JZ: error installing?
source("/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_Code/confidence_intervals/emma.R")

#get list of samples
phenotypes = read.table("/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt" ,header=T, sep="\t")
samples = phenotypes$"Animal_ID"

#make kinship.matrices object
#only need to do once; afterward can load Rdata file 
#load("/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/grms/oxinfo90/kinship_matrices.Rdata")

#K.chr <- vector("list", 19)
#names(K.chr) <- 1:19
k_dir="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/grms/oxinfo90/"
g_dir="/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/"
#n <- 3125920
#for( i in 1:19){
#  k = read.table(paste0(k_dir, "combined.chr", toString(i), ".oxinfo90.cXX.txt", sep=""), header=FALSE)
#  a = read.table(paste0(g_dir, "combined.chr", toString(i), ".oxinfo90.annot", sep=""), header=FALSE)
#  numSnps = n - nrow(a)
#  print(numSnps)
#  colnames(k) <- samples
#  rownames(k) <- samples
#  K.chr[[i]]$K = k
#  K.chr[[i]]$n = numSnps
#}
#
#kinship.matrices <- vector("list", 1)
#kinship.matrices[["K.chr"]] = K.chr
#
#save(kinship.matrices, file="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/grms/oxinfo90/kinship_matrices.Rdata")

#create mixed model objects
#phenotypeF <- "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt"
outdir <- "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/"
#mixed.model.heritabilities(kinship.matrices, phenotypes=NULL, out.dir=outdir,  normalise=TRUE, out.file="heritabilities_alleles.txt",phenotype.file=phenotypeF,mc.cores = 10)

#make dosage files
#phenotypes <- read.table("/u/home/j/jzou1115/project-zarlab/CFW/phenotypes/pheno.onlygenoed.combined.meta.txt", header=FALSE)
#samples <- phenotypes$V1
#
#print("making dosage rdata files...")
#for(chrm in 1:19){
#	print(chrm)
#	dosageF <- read.table(paste0(g_dir, "/combined.chr", toString(chrm), ".oxinfo90.dosages.gz"), header=F)
#	annot <- read.table(paste0(g_dir, "/combined.chr", toString(chrm),".oxinfo90.annot"), header=F)
#	colnames(annot) <- c("rsid", "bp", "chr")
#	#remove columns with snps and a1/a0
#	dosages_only <- dosageF[,4:3237] 
#
#	dosages <- list(ids = samples, dosages = dosages_only, map = annot)
#	save(dosages, file= paste0(g_dir, "/chr", toString(chrm), ".oxinfo90.Rdata"))
#}
#
#
##genome scans
#print("genome scans...")
outscans="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scans"
#mixed.model.pruned.dosage.scan.wrapper( mm.dir=outdir, 
#                                        out.dir=outscans, 
#                                        mc.cores=1, 
#                                        allele.dir=g_dir,
#                                        suffix=".oxinfo90.Rdata",
#                                        permute=100, 
#                                        decreasing=FALSE, 
#                                        per.chromosome=TRUE, 
#                                        type="allele" )
#

#CI
ci_out = "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/confidence_intervals"
#setwd(ci_out)

#simulate.confidence.intervals( qtl.file="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_perm.txt", 
#                               n.sim=10000, 
#                               window=5.0e6, 
#                               out.dir=ci_out, 
#                               per.chromosome=TRUE, 
#                               mixed.model=TRUE, 
#                               mm.dir=outdir, 
#                               dir=g_dir, 
#                               suffix=".oxinfo90.Rdata" , mc.cores1=10, mc.cores2=10, reverse=FALSE)

evaluate.confidence.intervals.logP( qtl.dir = ci_out, 
                                    scan.dir =outscans, pdf.file="CI.pdf", 
                                    quantiles=c(0.25,0.5,0.75,0.9,0.95), 
                                    dist.thresh=5.0e6, 
                                    CI.file =paste0(ci_out,"/CI_pruned.txt"),  #output file
									fine.dir = NULL,
                                    mc.cores=1)
