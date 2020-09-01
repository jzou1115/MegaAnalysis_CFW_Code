
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
    input = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/kinship_matrices.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt', "kinship" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/kinship_matrices.RData', "pheno" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/phenotypes/pheno.combined.txt'),
    output = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.activity30.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.activity5.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.activity.begin.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.activity.end.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.activity.middle.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.activity.total.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.basal.activity.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.bmd.a.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.bmd.n.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.bmd.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.decay.activity.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.EDL.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.baseline.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.context.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.context.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.cue.baseline.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.cue.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.cue.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.uncond.freeze.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.fc.uncond.freeze.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.gastroc.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.glucose.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.habit.diff.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.habit.ratio.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.plantaris.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.pp12.ppi.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.pp6.ppi.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.soleus.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.startle.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.tail.length.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.TA.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.tibia.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.weight.bmi.tibia.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm/emma.weight.RData'),
    params = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm', "outdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/mm'),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("outdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned', "chromosomes" = 'data/chromosomes.txt', "phenotypes" = 'data/phenotypes.txt', "qtls" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt', "kinship" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/grms/kinship_matrices.Rdata', "scans" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scans', "mm" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2', "qtl_ids" = 'data/qtls_decorrelated_info90.txt'),
    rule = 'mixedModelHeritabilities',
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
mixed.model.heritabilities(kinship.matrices, phenotypes=NULL, out.dir=outdir,  normalise=TRUE, out.file="heritabilities_alleles.txt",phenotype.file=phenotypeF,mc.cores = 1)

