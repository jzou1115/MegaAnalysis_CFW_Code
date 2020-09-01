#input/output from snakemake
dosageF <- snakemake@input[["dosage"]]
annotF <- snakemake@input[["annot"]]
pruneF <- snakemake@input[["prune"]]
phenoF <- snakemake@input[["pheno"]]

outF <- snakemake@output[["out"]]
outdir <- snakemake@params[["outdir"]]
if(!dir.exists(outdir)){
	mkdir(outdir)
}

#read in data
dosage <- read.table(dosageF, header=F)
annot <- read.table(annotF, header=F)
prune <- read.table(pruneF, header=F)
prune_snps <- unique(prune$V1)

#get list of samples
phenotypes = read.table(phenoF ,header=T, sep="\t")
samples = phenotypes$"Animal_ID"

#filter snps
keep <- match(prune_snps, annot$V1)
dosage <- dosage[keep,]
annot <- annot[keep,]

#make RData output file
colnames(annot) <- c("rsid", "bp", "chr")
dosages <- list(ids = samples, dosages = dosage[,4:ncol(dosage)], map = annot)
save(dosages, file= outF)


