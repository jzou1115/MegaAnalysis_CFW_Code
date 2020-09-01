args <- commandArgs(trailing=T)

f = args[1]
print(f)
z = read.table(args[2], header=T)

prefix = args[3]

load(f)

pip <- fitted_rss$pip
write.table(pip, file=paste0(prefix, "_pip.txt"), row.names=F, col.names=F)

sets <- fitted_rss$sets$cs

ind <- c()
for(set in names(sets)){
	temp <-array(as.numeric(unlist( fitted_rss$sets$cs[set])))
	ind <- c(ind, temp)
}
ind <- unique(ind)
keep <- z[ind, "rs"]


write.table(keep, file=paste0(prefix, "_set.txt"), row.names=F, col.names=F, quote=F)
