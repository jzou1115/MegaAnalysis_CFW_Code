library(susieR)
options(echo=TRUE)
#args <- commandArgs(trailingOnly=FALSE)

zf <-  snakemake@input[["z"]]
ldf <- snakemake@input[["ld"]]
prefix <- snakemake@params[["prefix"]]

print(args)

z = read.table(zf, header=T, sep="\t")
print(dim(z))
rs = z$rs
print(length(rs))

ld = read.table(ldf, header=T, sep="\t")
print(dim(ld))
ld <- data.matrix(ld)

# Need to compute R (ld between snps)
fitted_rss <- susie_rss(z$z, ld, L = 1,
                        estimate_residual_variance = TRUE, 
                        estimate_prior_variance = TRUE)
print(summary(fitted_rss)$cs)


save(fitted_rss, file=paste0(prefix, ".RData"))


variants = c()
for( set in names(fitted_rss$sets$cs)){
	vars = fitted_rss$sets$cs[set] 
	variants <- c(variants, vars)
}
variants = unlist(variants)

print(variants)
rs = rs[variants]
write.table(rs, file=paste0(prefix, "_snps.txt"), row.names=F, col.names=F, quote=F)

plot_f <- paste0(prefix, ".png")
png(plot_f)
susie_plot(fitted_rss, y="PIP")
dev.off()
