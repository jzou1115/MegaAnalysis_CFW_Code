
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/combined.chr2.oxinfo90.dosages.gz', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/combined.chr2.oxinfo90.annot', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/combined.chr2.oxinfo90.txt', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt', "dosage" = '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/combined.chr2.oxinfo90.dosages.gz', "annot" = '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/combined.chr2.oxinfo90.annot', "prune" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/combined.chr2.oxinfo90.txt', "pheno" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt'),
    output = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr2.oxinfo90.Rdata', "out" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr2.oxinfo90.Rdata'),
    params = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes', "outdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes'),
    wildcards = list('2', "chr" = '2'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("outdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned', "chromosomes" = 'data/chromosomes.txt', "phenotypes" = 'data/phenotypes.txt', "qtls" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold3.txt', "kinship" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/grms/kinship_matrices.Rdata', "scans" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scans', "mm" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2', "qtl_ids" = 'data/qtls_decorrelated_info90_all.txt'),
    rule = 'subsetPrunedSNPs',
    bench_iteration = as.numeric(NA),
    scriptdir = '/u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/scripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## Original script #########
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


