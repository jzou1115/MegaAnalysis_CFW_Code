library(parallel)
library(graphics)
source("scripts/outbred.pruned.R")
source("scripts/emma.R")

scandir <- snakemake@params[["scan_dir"]]
ci_out <- snakemake@params[["ci_out"]]

evaluate.confidence.intervals.logP( qtl.dir = ci_out, 
                                    scan.dir =scandir, pdf.file="CI.pdf", 
                                    quantiles=c(0.25,0.5,0.75,0.9,0.95), 
                                    dist.thresh=5.0e6, 
                                    CI.file =paste0(ci_out,"/CI_pruned.txt"),  #output file
									fine.dir = NULL,
                                    mc.cores=1)
