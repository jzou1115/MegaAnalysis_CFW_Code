#JZ: 
#Input: 
#1. chr- chromosome of qtl performing simulation for (lead snp?) 
#2. bp - position of qtl performing simulation for (lead snp?) 
#3. mm - ?
#4. dosages - only ever use  dosages$dosages, which is probably a nxm matrix, where n is number of samples and m is number of snps
#5. n.sim - number of simulations to perform for qtl
#6. window - largest distance to consider for locus
#7. out.dir - output directory
#8. mc.cores - number of cores to use?  Doesn't actually seem to be used in this function

simulate.qtl.confidence.interval <- function( chr, bp, mm, dosages, n.sim=1000, window=5.0e6, out.dir="./",  mc.cores=48 ) {
  
  did = dosages$ids
  map = dosages$map
  
  #JZ: What are the adjustments for?   
  if ( mm$per.chromosome ) { # use chromosome-specific adjustments
    mm.use = mm$mm.chr[[chr]]
    cat("Using per.chromosome\n")
  }
  else {
    mm.use = mm$mm.all # use a common adjustment
  }
  #JZ: subset SNPs for some reason.  perhaps to intersect snps between mm.use and dosages$ids.
  use = match( mm.use$ids, did, nomatch=0)
  D = dosages$dosages
  D = D[,use]
  
  #JZ: find SNPs within the window (default = 5MB) and subset the dosages and map files
  snps.subset2 = which( abs(map$bp -bp) < window ) # outer subset defining the search window
  D.qtl = D[snps.subset2,]
  map.qtl = map[snps.subset2,]
  bp.idx = which( map.qtl$bp == bp )
  if ( length(bp.idx) == 0 ) {
    warning("ERROR  coordinate" ,chr, bp , "not found\n")
    return(NULL);
  }
  
  #JZ: n is number of SNPs within the window
  n = nrow(D.qtl)
  
  #JZ: what is the point of this?  compute the association statistic?  What is the multiplier?
  if ( mm$mixed.model == TRUE ) {
    D.qtl = t(mm.use$multiplier) %*% t(D.qtl) # dimension nsub (rows) * nsnps (cols)
    r = cor( D.qtl,  mm.use$y )
    lm.snp = lm( mm.use$y ~ D.qtl[,bp.idx],y=TRUE)
  }
  else {
    D.qtl = t(D.qtl)
    r = cor( D.qtl,  mm.use$y.original )
    lm.snp = lm( mm.use$y.original ~ D.qtl[,bp.idx], y=TRUE)
  }
  
  
  r2 = r*r
  r2 = ifelse(is.na(r2), 0, r2) #JZ: why would this be NA?
  
  phen.r2 = r2[bp.idx]
  phen = lm.snp$y
  phen.resid = resid(lm.snp) #JZ:empirical residuals of association statistic
  
  #JZ: find all snps that are within window/2 and correlated with phenotype.  Why window/2 when line 31 uses window?
  snps.subset1 = which( abs(map.qtl$bp -bp) < window/2 & r2>0 )
  
  #JZ: Select 1000 random snps to be causal and perform simulations with. 
  snps.sample = sort(sample(snps.subset1, n.sim, replace=TRUE)) # the snp sites to simulate from
  
  #JZ: compute summary statistics for 1000 snps
  var.snp = apply( D.qtl, 2, var, na.rm=TRUE ) # variability of the dosages at each snp
  snps.beta = r[bp.idx]*sqrt(var(phen))/sqrt(var.snp)
  snps.sample = c( bp.idx, snps.sample )
  
  #JZ: apply simulation framework to 1000 snps.  Assume that snp is causal; 
  sims = mclapply(snps.sample, function( snp, snps.beta, var.snp, phen.resid, D.qtl, map.qtl, target.r2 ) {
    #            y = sample(phen.resid) + snps.beta[snp]*D.qtl[snp,]
    #JZ: simulate phenotype assuming that SNP is causal.  randomly assign empirical residuals to different samples.  (doesn't seem equivalent to bootstrapping samples?)
    y = sample(phen.resid) + snps.beta[snp]*D.qtl[,snp] 
    #JZ: compute association statistics using simulated phenotype
    nsub = length(y)
    r = cor( D.qtl,  y )
    r2 = r*r
    r2 = ifelse(is.na(r2), 0, r2)
    F = (r2*(nsub-2))/(1-r2+1.0e-20)
    logP = -pf(F, 1, nsub-2, lower.tail = FALSE, log.p = TRUE)/log(10)
    logP[snp] = -logP[snp] # exclude the causal snp
    #JZ: find SNP (not including causal snp) with the most significant association statistic
    w.max = which.max( logP )
    #JZ: find bp location of most significant association statistic in simulations
    bp.max = map.qtl$bp[w.max]
    #JZ: compute difference distance between causal snp and the  most significant association statistic in simulations
    bp.delta = bp.max-map.qtl$bp[snp]
    return( c(snp, map.qtl$bp[snp], logP[snp], map.qtl$bp[w.max], logP[w.max],  bp.delta, snps.beta[snp], var.snp[snp] ))
  }, snps.beta, var.snp, phen.resid, D.qtl, map.qtl, phen.r2, mc.cores=mc.cores )
  
  sims = do.call( "rbind", sims )
  colnames(sims) = c("snp", "bp", "snp.logP", "max.bp", "max.logP", "bp.delta", "snp.beta", "snp.var")
  return(t(sims))
  
}