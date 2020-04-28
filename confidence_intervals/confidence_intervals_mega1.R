library(parallel)
#library(emma)
library(graphics)
source("/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_Code/confidence_intervals/outbred.pruned.R")
#library("RevoUtilsMath") #JZ: error installing?
source("/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_Code/confidence_intervals/emma.R")
args = commandArgs(trailingOnly=TRUE)

k_dir="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/grms/"
g_dir="/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/"

 mixed.model.data.file <- args[1]
#genome scans
print("genome scans...")
outscans="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scans"

#mixed.model.pruned.dosage.scan.wrapper( mm.dir=outdir, 
#                                        out.dir=outscans, 
#                                        mc.cores=5, 
#                                        allele.dir=g_dir,
#                                        suffix=".oxinfo90.Rdata",
#                                        permute=100, 
#                                        decreasing=FALSE, 
#                                        per.chromosome=TRUE, 
#                                        type="allele" )

#mixed.model.pruned.dosage.scan.genome( mixed.model.data.file, 
#                                                   haplotype.dir="/Net/sparse/data/scratch/MOTT/OUTBRED/haplotype.dosages.0.02", 
#                                                   suffix=".oxinfo90.Rdata", 
#                                                   chrs=c(1:19), 
#                                                   out.dir=outscans,  
#                                                   permute=100, 
#                                                   mixed.model=TRUE, 
#                                                   per.chromosome=TRUE, 
#                                                   type="allele" )
#

mm = NULL
#cat(sum(fe), length(fe), "loading ",  mixed.model.data.file, "\n")
        
load( mixed.model.data.file) # loads mixed.model.data into object called mm. This is produced by a previous call to mixed.model.heritabilities()
        
cat("loaded OK ",  mixed.model.data.file, "\n")

mm$mixed.model = TRUE
mm$per.chromosome = TRUE
mm$perms = 100

for (chr in 1:19){
		print(chr)
		mixed.model.pruned.allele.dosage.scan.per.chr( chr, 
                                                           mm, 
                                                           out.dir=outscans, 
                                                           dir=g_dir, 
                                                           suffix=".oxinfo90.Rdata", 
                                                           chicago=FALSE, 
                                                           mc.cores=1 )

}
