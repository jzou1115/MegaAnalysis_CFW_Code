library(data.table)

args = commandArgs(trailingOnly=TRUE)

pheno = args[1]

data <- list()

for (n in 1:19) {
	f <- paste0("/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association/gwas_out_cov/ox1_pheno_", toString(pheno),"_chr",toString(n), ".assoc.txt") 
	data[[n]] = read.table(f, header=T)
}

data_all <- rbindlist(data)

observed <- sort(data_all$p_lrt)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))


outfile=paste0("/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association/gwas_out_cov/ox1_pheno_", toString(pheno), "_qqplot.png")
png(outfile, width=6, height=6, units="in", res = 300)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

#data_all$s1 = data_all$beta / data_all$se

