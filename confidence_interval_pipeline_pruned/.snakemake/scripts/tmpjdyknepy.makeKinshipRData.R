
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
    input = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr1.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr2.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr3.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr4.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr5.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr6.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr7.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr8.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr9.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr10.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr11.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr12.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr13.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr14.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr15.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr16.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr17.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr18.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr19.oxinfo90.Rdata', "genotypes" = c('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr1.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr2.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr3.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr4.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr5.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr6.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr7.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr8.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr9.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr10.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr11.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr12.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr13.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr14.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr15.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr16.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr17.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr18.oxinfo90.Rdata', '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr19.oxinfo90.Rdata')),
    output = list('/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/kinship_matrices.RData', "output" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/kinship_matrices.RData'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("outdir" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned', "chromosomes" = 'data/chromosomes.txt', "phenotypes" = 'data/phenotypes.txt', "qtls" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt', "kinship" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/grms/kinship_matrices.Rdata', "scans" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scans', "mm" = '/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2', "qtl_ids" = 'data/qtls_decorrelated_info90.txt'),
    rule = 'makeKinship',
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
genotypes <- snakemake@input[["genotypes"]]
out <- snakemake@params[["output"]]
print(out)

kinship.matrices <- list()
kinship.matrices[["K.chr"]] <- list()

for(genoF in genotypes){
	load(genoF)
	print(genoF)

	base <- basename(genoF)
	tokens <- strsplit(base, split=".", fixed=T)
	chrm <- tokens[[1]][1]
	chrm <- substr(chrm, 4, nchar(chrm))

	dos = dosages$dosages
	K = cor(dos, use="complete.obs")
	n = nrow(dos)

	kinship.matrices$K.chr[[chrm]] = list("K"=K, "n"=n)	
	break

}

save(kinship.matrices, file=out)

#> load("/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/chr1.oxinfo90.Rdata")
#> ls()
#[1] "dosages"
#> names(dosages)
#[1] "ids"     "dosages" "map"  
#> names(dosages$map)
#[1] "rsid" "bp"   "chr" 
#> length(dosages$ids)
#[1] 3234
#> dim(dosages$dosages)


#> ls()
#[1] "kinship.matrices"
#> names(kinship.matrices)
#[1] ""      "K.chr"
#> names(kinship.matrices$K.chr)
# [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
#[16] "16" "17" "18" "19"
#> names(kinship.matrices$K.chr$1)
#Error: unexpected numeric constant in "names(kinship.matrices$K.chr$1"
#> names(kinship.matrices$K.chr$"1")
#[1] "K" "n"
#> names(kinship.matrices$K.chr$"1"$K)
#> names(kinship.matrices$K.chr$"1"$n)
#NULL
#> dim(kinship.matrices$K.chr$"1"$K)
#[1] 3234 3234
#> kinship.matrices$K.chr$"1"$n
#[1] 2874436
