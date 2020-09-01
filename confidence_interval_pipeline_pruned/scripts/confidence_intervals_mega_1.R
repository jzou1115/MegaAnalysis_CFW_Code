library(graphics)
source("scripts/outbred.pruned.R")
source("scripts/emma.R")

#snakemake input
mixed.model.data.file  <-  snakemake@input[["mmF"]]
outscans <- snakemake@params[["scan_dir"]]
g_dir <- snakemake@params[["g_dir"]]

if(!dir.exists(outscans)){
    mkdir(outdir)
}


mm = NULL

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
