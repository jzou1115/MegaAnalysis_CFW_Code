genotypes <- snakemake@input[["genotypes"]]
out <- snakemake@output[["output"]]
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

	ids = dosages$ids
	dos = dosages$dosages
	K = cov(dos, use="complete.obs")
	row.names(K) = ids
	colnames(K) = ids
	n = nrow(dos)

	kinship.matrices$K.chr[[chrm]] = list("K"=K, "n"=n)	

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
