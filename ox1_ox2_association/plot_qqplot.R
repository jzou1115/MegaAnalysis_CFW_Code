library(data.table)

args = commandArgs(trailingOnly=TRUE)

inputF = args[1]
outF = args[2]

data_all <- read.table(inputF, header=T)

observed <- sort(data_all$p_lrt)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))


png(outF, width=6, height=6, units="in", res = 300)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

#data_all$s1 = data_all$beta / data_all$se

