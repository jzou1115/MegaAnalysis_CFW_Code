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
	i_chr = length(tokens[[1]])-1
	chr = as.numeric(tokens[[1]][i_chr])
	bp = as.numeric(tokens[[1]][length(tokens[[1]])])
	p = paste(tokens[[1]][1:i_chr-1], sep=".", collapse=".")

#	p=tokens[[1]][1]
#	chr = as.numeric(tokens[[1]][2])
#	bp = as.numeric(tokens[[1]][3])
	i = which((qtls$phenotype==p) & (qtls$chr==chr) & (qtls$bp==bp) )
	print(i)
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
	sim.data = list( qtl=qtls[i,], n.sim=n.sim, window=window, sim=sim )
	save(sim.data,file=file)
} else{
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
