
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
    input = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr1.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr2.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr3.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr4.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr5.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr6.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr7.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr8.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr9.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr10.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr11.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr12.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr13.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr14.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr15.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr16.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr17.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr18.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr19.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity30.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity5.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.begin.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.end.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.middle.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.total.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.basal.activity.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.a.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.n.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.decay.activity.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.EDL.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.baseline.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.baseline.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.gastroc.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.glucose.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.diff.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.ratio.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.plantaris.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp12.ppi.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp6.ppi.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.soleus.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.startle.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tail.length.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.TA.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tibia.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.bmi.tibia.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.RData', "qtls" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt', "genoFiles" = c('/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr1.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr2.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr3.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr4.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr5.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr6.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr7.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr8.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr9.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr10.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr11.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr12.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr13.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr14.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr15.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr16.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr17.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr18.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr19.oxinfo90.Rdata'), "mmFiles" = c('/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity30.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity5.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.begin.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.end.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.middle.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.total.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.basal.activity.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.a.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.n.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.decay.activity.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.EDL.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.baseline.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.baseline.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.corr.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.gastroc.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.glucose.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.diff.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.ratio.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.plantaris.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp12.ppi.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp6.ppi.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.soleus.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.startle.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tail.length.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.TA.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tibia.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.bmi.tibia.RData', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.RData')),
    output = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals/tibia.12.91175287.CI.RData'),
    params = list('tibia.12.91175287', '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals', "qtl" = 'tibia.12.91175287', "g_dir" = '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90', "mm_dir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2', "ci_out" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals'),
    wildcards = list('tibia.12.91175287', "q" = 'tibia.12.91175287'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("outdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out', "mmdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2', "gdir" = '/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90', "kdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/grms', "cidir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals', "scandir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scans', "chromosomes" = 'data/chromosomes.txt', "phenotypes" = 'data/phenotypes.txt', "qtls" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt', "qtl_ids" = 'data/qtls_decorrelated_info90.txt'),
    rule = 'simulateConfidenceIntervals',
    bench_iteration = as.numeric(NA),
    scriptdir = '/u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/scripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## Original script #########
library(parallel)
library(graphics)
source("scripts/outbred.pruned.R")
source("scripts/emma.R")

qtl <- snakemake@params[["qtl"]]
qtl.file <- snakemake@input[["qtls"]]
g_dir <- snakemake@params[["g_dir"]]
mm_dir <- snakemake@params[["mm_dir"]]
outdir <- snakemake@params[["ci_out"]]
if(!dir.exists(outdir)){
    mkdir(outdir)
}

setwd(outdir)

n.sim=10000
window=5.0e6
out.dir=outdir
per.chromosome=TRUE
mixed.model=TRUE
mm.dir=mm_dir
dir=g_dir
suffix=".oxinfo90.Rdata"
mc.cores1=1
mc.cores2=1
reverse=FALSE

#simulate.confidence.intervals
qtls = read.table(qtl.file, header=T)


if ( ! file.exists( out.dir ) ) dir.create( out.dir)

file = paste(out.dir, "/", qtl, ".CI.RData", sep="")
if ( ! file.exists(file)) {

	tokens = strsplit(qtl,  ".", fixed = TRUE)
	p=tokens[[1]][1]
	chr = as.numeric(tokens[[1]][2])
	bp = as.numeric(tokens[[1]][3])
	i = which((qtls$phenotype==p) & (qtls$chr==chr) & (qtls$bp==bp) )
	if(length(i)>1){
		print("Warning: qtl not unique")
		i = i[1]
	}	
	
	dosages = load.pruned.allele.dosages.for.chr2( chr=chr, dir=dir, suffix=suffix, chicago=chicago  )
	
	mm.file = paste(mm.dir, "/emma.", p, ".RData", sep="")
	cat(mm.file, "\n")
	load(mm.file)
	mm$per.chromosome = per.chromosome
	mm$mixed.model = mixed.model
	cat("attempting ", file, "\n")
	sim = simulate.qtl.confidence.interval( chr, bp, mm, dosages, n.sim=n.sim, window=window, out.dir=out.dir, mc.cores=mc.cores2 )
	sim.data = list( qtl=qtls.chr[i,], n.sim=n.sim, window=window, sim=sim )
	save(sim.data,file=file)
}
else{
	cat( file, " exists\n")
}		


#simulate.confidence.intervals <- function( qtl.file="/Net/dense/data/nicod/QTLs_areas_1Dec14.txt", n.sim=10000, window=5.0e6, out.dir="./", per.chromosome=TRUE, mixed.model=TRUE, mm.dir="/Net/dense/data/nicod/Ultimate_Analysis_R5.5/norm_scaled/", dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData" , mc.cores1=10, mc.cores2=10, reverse=FALSE, chicago=chicago) {
#  #JZ: replace w/ read table so that input is read in correctly and not shifted over by 1 col  
#  #qtls = read.delim(qtl.file, stringsAsFactors=FALSE)
#    qtls = read.table(qtl.file, header=T)
#    print(qtls)
#    
#    if ( ! file.exists( out.dir ) ) dir.create( out.dir)
#    
#    chrs = sort(unique(qtls$chr) )
#    if ( reverse ) chrs = rev(chrs)
#    for( chr in chrs ) {
#        
#        qtls.chr = qtls[qtls$chr==chr,]
#        files=paste(out.dir, "/", qtls.chr$phenotype, ".", chr, ".", qtls.chr$bp, ".CI.RData", sep="")
#        file.existance = sapply( files, file.exists )
#        qtls.chr = qtls.chr[!file.existance,]
#
#        if ( nrow(qtls.chr) > 0 ) {
#             dosages = load.pruned.allele.dosages.for.chr2( chr=chr, dir=dir, suffix=suffix, chicago=chicago  )
#            mclapply( 1:nrow(qtls.chr), function( i, qtls.chr, out.dir, mm.dir, chr, dosages, n.sim, window, mc.cores2=mc.cores2) {
#                p = qtls.chr$phenotype[i]
#    
#                file=paste(out.dir, "/", p, ".", chr, ".", qtls.chr$bp[i], ".CI.RData", sep="")
#                if ( ! file.exists(file)) {
#                    mm.file = paste(mm.dir, "/emma.", p, ".RData", sep="")
#                    cat(mm.file, "\n")
#                    load(mm.file)
#                    mm$per.chromosome = per.chromosome
#                    mm$mixed.model = mixed.model
#                    
#                    cat("attempting ", file, "\n")
#                    sim = simulate.qtl.confidence.interval( chr, qtls.chr$bp[i], mm, dosages, n.sim=n.sim, window=window, out.dir=out.dir, mc.cores=mc.cores2 )
#                    sim.data = list( qtl=qtls.chr[i,], n.sim=n.sim, window=window, sim=sim )
#                    save(sim.data,file=file)
#                }
#                else {
#                    cat( file, " exists\n")
#                }
#            },  qtls.chr, out.dir, mm.dir, chr, dosages, n.sim, window, mc.cores2=mc.cores2, mc.cores=mc.cores1 )
#        }
#    }
#}
#        
