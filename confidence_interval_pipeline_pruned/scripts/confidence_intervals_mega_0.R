library(graphics)
source("scripts/outbred.pruned.R")
source("scripts/emma.R")

#print("making ixed model objects...")
#phenotypeF <- "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt"
#outdir <- "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity30"

k <- snakemake@input[["kinship"]]
phenotypeF <-snakemake@input[["pheno"]]
outdir <- snakemake@params[["outdir"]]


load(k)
mixed.model.heritabilities(kinship.matrices, phenotypes=NULL, out.dir=outdir,  normalise=TRUE, out.file="heritabilities_alleles.txt",phenotype.file=phenotypeF,mc.cores = 10)

