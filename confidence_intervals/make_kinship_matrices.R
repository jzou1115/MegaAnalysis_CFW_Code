args = commandArgs(trailingOnly=TRUE)

dosageF = args[1]
outF = args[2]

dosages = read.table(dosageF, header=F)
dosages = dosages[,4:ncol(dosages)] #remove first 3 columns with rsid/a1/a0
print(dim(dosages))

n = nrow(dosages)
downsample=0.05
size = round(downsample*n)
use = sample(n, size=size, replace=FALSE )
dosages = dosages[use,]

dosages = t(dosages)
dosages = scale(dosages)
dosages = t(dosages)

k = cor(dosages, use="complete.obs")
print(dim(k))

dosages = dosages 

save(k, file = outF)
