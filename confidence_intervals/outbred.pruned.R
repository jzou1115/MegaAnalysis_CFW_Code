library(parallel)
#library(emma)
library(graphics)
source("outbred.rwd.R")
source("emma.R")
#JZ:troubling installing
#library("RevoUtilsMath")
#setMKLthreads(1)

load.pruned.allele.dosages.for.chr <- function( chr, dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData", chicago=FALSE ) {

    file = paste(dir, "/chr", chr, suffix, sep="")
    cat("loading", file, "\n")
    load(file)
    map=pruned_pos
    names(map) = c( "chr", "bp", "ref", "alt")
    nameList = unlist(nameList)
    if ( chicago ) {
        ids = sub(".bam$", "", nameList, perl=TRUE)
        ids = sub("_recal$", "", nameList, perl=TRUE)
    } else {
        ids = sub("_recal.bam$", "", nameList, perl=TRUE)
        ids = sub("___", "/", ids, perl=TRUE)
    }
    colnames(pruned_dosages) = ids
    rownames(pruned_dosages) = map$bp
    data=list( chr=chr, file=file, nameList=nameList, ids=ids, dosages=pruned_dosages, map=map, type="allele" )
    return(data)
}

#JZ: works for my dosage RData file
load.pruned.allele.dosages.for.chr2 <- function( chr, dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData", chicago=FALSE ) {
  
  file = paste(dir, "/chr", chr, suffix, sep="")
  cat("loading", file, "\n")
  load(file)
 
  colnames(dosages$dosages) = dosages$ids
  rownames(dosages$dosages) = dosages$map$bp
  data=list( chr=dosages$chr, file=dosages$file, ids=dosages$ids, dosages=dosages$dosages, map=dosages$map, type="allele" )
  return(data)
}

load.complete.allele.dosages.for.locus <- function( ff.data, chr, from.bp, to.bp ) {

    if ( ff.data$chr == chr & ff.data$type == "allele" ) {
        which.snps = which( ff.data$map$bp >= from.bp & ff.data$map$bp <= to.bp )
        if ( length(which.snps) > 0 ) {
            dosages.locus = ff.data$dosages[which.snps,]
            data = list( chr=chr, ids=ff.data$ids, Nsnps=length(which.snps), map=ff.data$map[which.snps,], dosages=dosages.locus, type="allele", range=c(from.bp, to.bp))
            return(data)
        }
        else {
            warning("No SNPs in interval", chr, from.bp, to.bp, "\n")
        }
    }
    else {
        warning( "Incorrect chromosome or type",  ff.data$chr, chr , ff.data$type, "allele" , "\n")
    }
    return(NULL)
}

make.haplotype.dosages <- function( chrs=c(1:19,"X"), version="5.5", hap.subset=0.05, outdir="/Net/sparse/data/scratch/MOTT/OUTBRED/haplotype.dosages.0.05", suffix=".haplotypes.RData" ) {

    if ( ! file.exists( outdir ) ) dir.create( outdir)
    mclapply(  chrs , function(chr, outdir, suffix, hap.subset ) {
        outfile = paste(outdir, "/chr", chr, suffix, sep="")
        if ( ! file.exists(outfile)) {
            cat( "chr ", chr, "\n")
            dosages = rwd.load.dosages.for.chr( chr, type="haplotype", version="5.5", hap.subset=hap.subset, min.maf=0.001 )
            d.dim = dim(dosages)
            nhap = d.dim[3]-1
            dosages$dosages = dosages$dosages[,,1:nhap]
            save(dosages,file=outfile)
            cat( "wrote", outfile, "\n"  ) 
        }
    }, outdir, suffix, hap.subset, mc.cores=2 )
}

load.haplotype.dosages.for.chr <- function( chr, dir="/Net/sparse/data/scratch/MOTT/OUTBRED/haplotype.dosages.0.02", suffix=".haplotypes.RData"  ) {
    file = paste(dir, "/chr", chr, suffix, sep="")
    load(file)
    return(dosages)
                                        #    dosages = list( chr=chr, ids=ids, intervals=intervals, map=map, dosages=dd, type=type, hap.subset=hap.subset))
}

haplotype.dosage.stats <- function( chrs=c(1:19,"X"), dir="/Net/sparse/data/scratch/MOTT/OUTBRED/haplotype.dosages.0.02", pdf.file="haplotype.plots.pdf" ) {

    pdf( pdf.file )
    par(mfrow=c(3,1))
    zz = NULL
    y = c( 0.5, 0.6, 0.7, 0.8, 0.9 )
    
    for( chr in chrs ) {
        dosages = load.haplotype.dosages.for.chr(chr , dir=dir );
        d = dosages$dosages
        d.mean = apply( d, c(1,3), mean )
        dd = matrix( nrow=nrow(d.mean), ncol=1+ncol(d.mean))
        dd[,1:ncol(d.mean)] = d.mean
        dd[,ncol(dd)] = 1-apply(d.mean,1,sum)
        d.sorted = t(apply( dd, 1, sort, decreasing=TRUE ) )
        mx = d.sorted[,1]+d.sorted[,2]
        ent = apply( dd, 1, function(x) { -sum( log2(x+1.0e-10)*x )} )
        data=data.frame(chr=chr, dosages$map, max.2=mx, entropy=ent)
        save(data, file=paste( "haplotype.stats.", chr, ".RData", sep=""))
        plot(dosages$map$bp/1.0e6, mx, t="l", ylim=c(0,1), col="black", xlab=paste(chr, "Mb"), ylab="entropy, max2probs" )
#        plot(dosages$map$bp/1.0e6, mx, t="l", ylim=c(0,2), col="black", xlab=paste(chr, "Mb"), ylab="entropy, max2probs" )
#        lines(dosages$map$bp/1.0e6, ent, col="red")
        
        z = sapply( y, function( x, mx) { mean(mx>x) }, mx )
        zz = rbind(zz, z )
    }
    zz = data.frame( chr=chrs, zz )
    names(zz) = c( "chr", y )
    print(zz)
    dev.off()
    return(zz)
}


make.pruned.kinship.matrix <- function (chrs=c(1:19,"X"), type="allele", file.prefix=paste("pruned.kinship",  type,  sep="."), dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData", scale.dosages=TRUE, mc.cores=20, chicago=FALSE) {

    K.chr = list()
    K.chr = mclapply( chrs, function( chr, type, dir, suffix, chicago ) {
                                        #                                         K = lapply( chrs, function( chr, type ) {
        pruned.dosages = load.pruned.allele.dosages.for.chr( chr=chr, dir=dir, suffix=suffix, chicago=chicago )
        D = pruned.dosages$dosages
        if ( scale.dosages == TRUE ) {
            tD = t(D)
            tD = scale(tD)
            tD = ifelse( is.na(tD), 0, tD)
            D = t(tD)
            cat( "rescaled kinship matrix\n")
        }
        n = nrow(D)
#        D=D-rowSums(D)/ncol(D)
        n = nrow(D)
        cat("chr ", chr, n, "\n")
        k = cor(D)
        rownames(k) = as.character(pruned.dosages$ids)
        colnames(k) = as.character(pruned.dosages$ids)
        print(k[1:3,1:3])
        return( list( K=k, n=n))
    }, type, dir, suffix, chicago, mc.cores=mc.cores)
    names(K.chr) = chrs
    K.all = NULL
    K.n = 0
    for ( chr in names(K.chr)) {
        cat("merging chr ", chr, "\n")
        if ( is.null(K.all)) {
            K.all = K.chr[[chr]]$K
            n.all = K.chr[[chr]]$n
        }
        else {
            K.all = K.all + K.chr[[chr]]$K*K.chr[[chr]]$n
            K.n = K.n + K.chr[[chr]]$n
        }
    }
    K.all = K.all/K.n
    kinship2gcta(K.all, file.prefix)
    file=paste(file.prefix, ".RData", sep="")
    kinship.matrices = list( K.all = K.all, K.n=K.n, K.chr=K.chr, type=type, filter=filter, scale.dosages=scale.dosages ) # Note output has been changed to give per chromosome files as well as total kinship
    save( kinship.matrices, file=file)
    cat("wrote kinship matrix to", file, "\n")
}

                                        # mixed.model.pruned.dosage.scan.wrapper( dir="./kinship/", out.dir="./haplotype.permuted.per.chrom/", mc.cores=40, permute=100, decreasing=FALSE, mixed.model=TRUE, per.chromosome=TRUE, type="haplotype" )

mixed.model.pruned.dosage.scan.wrapper <- function( mm.dir, out.dir, mc.cores=20, permute=0, decreasing=FALSE, mixed.model=TRUE, per.chromosome=TRUE, type="allele", haplotype.dir="/Net/sparse/data/scratch/MOTT/OUTBRED/STATShaplotype.dosages.0.02", allele.dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData" ) {
    files=sort(list.files( path=mm.dir, pattern=".RData", full=TRUE), decreasing=decreasing )
                                        #    files = sample(files)
    mclapply ( files , function ( f, out.dir, permute, mixed.model, per.chromosome, type, haplotype.dir, dir, suffix ) {
        if ( FALSE ) {
            x = grep("ALP", f )
            cat( f, x,"\n")
            if ( !is.na(x[1])) {
                cat( "matches\n")
                
                mixed.model.pruned.dosage.scan.genome( f, out.dir=out.dir, permute=permute, mixed.model=mixed.model, per.chromosome=per.chromosome, type=type , haplotype.dir=haplotype.dir, dir=dir, suffix=suffix  )
            }
        }
        else {
            mixed.model.pruned.dosage.scan.genome( f, out.dir=out.dir, permute=permute, mixed.model=mixed.model, per.chromosome=per.chromosome, type=type, haplotype.dir=haplotype.dir,  dir=dir, suffix=suffix  )
        }            
    },out.dir, permute, mixed.model, per.chromosome, type=type, haplotype.dir=haplotype.dir ,  dir=allele.dir, suffix=suffix, mc.cores=mc.cores )
}



mixed.model.pruned.dosage.scan.genome <- function( mixed.model.data.file, haplotype.dir="/Net/sparse/data/scratch/MOTT/OUTBRED/haplotype.dosages.0.02", dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData", chrs=c(1:19,"X"), out.dir=NULL,  permute=0, mixed.model=TRUE, per.chromosome=TRUE, type="allele" ) {

    if ( !is.null ( out.dir ) ) {
        if (! file.exists(out.dir))
            dir.create( out.dir )
    }
    else
        out.dir="./"

    txt = sub("^.+emma.", "", mixed.model.data.file)
    txt = sub(".RData$", "", txt)
    
    files=paste(out.dir, "/", txt, ".", chrs, ".RData", sep="")
    
    if ( permute>0) {
        files=paste(out.dir, "/", txt, ".", chrs, ".permute.RData", sep="")
    }
    
    fe = file.exists(files)
    
    if ( sum(fe) < length(fe) ) {
        
        mm = NULL
        cat(sum(fe), length(fe), "loading ",  mixed.model.data.file, "\n")
        
        load( mixed.model.data.file) # loads mixed.model.data into object called mm. This is produced by a previous call to mixed.model.heritabilities()
        
        cat("loaded OK ",  mixed.model.data.file, "\n")
        
        mm$mixed.model = mixed.model
        mm$per.chromosome = per.chromosome
        
        mm$perms = permute
        
        cat("ATTEMPTING", mm$p, "\n")
        if ( type == "allele" ) {
            lapply( chrs, mixed.model.pruned.allele.dosage.scan.per.chr, mm, dir=dir, suffix=suffix,out.dir=out.dir )
        }
        else if ( type == "haplotype" ) {
            lapply( chrs, mixed.model.pruned.haplotype.dosage.scan.per.chr, mm, out.dir=out.dir, haplotype.dir=haplotype.dir )
        }
    }
}

mixed.model.pruned.haplotype.dosage.scan.per.chr <- function( chr, mm, out.dir=NULL,  haplotype.dir="/Net/sparse/data/scratch/MOTT/OUTBRED/haplotype.dosages.0.02", mc.cores=20 ) {

    file=paste(out.dir, "/", mm$p, ".", chr, ".RData", sep="")
    if ( mm$perms>0 ) {
        file=paste(out.dir, "/", mm$p, ".", chr, ".permute.RData", sep="")
    }
    
    if ( ! file.exists(file)) {
        cat("attempting ", file, "\n")
        dosages = load.haplotype.dosages.for.chr( chr=chr, dir=haplotype.dir )
        did = dosages$ids
        map = dosages$map
        if ( mm$per.chromosome ) { # use chromosome-specific adjustments
            mm.use = mm$mm.chr[[chr]]
            cat("Using per.chromosome\n")           
        }
        else {
            mm.use = mm$mm.all # use a common adjustment
        }

        use = match( mm.use$ids, did, nomatch=0)
        D = dosages$dosages
        D = D[,use,]
        D.dim = dim(D)
        n = D.dim[2]
        
        if ( mm$mixed.model == TRUE ) {
            if ( mm$perms>0 ) {
                yy = sapply( 1:mm$perms, function( i,y) { sample(y, replace=FALSE) }, mm.use$y )  # matrix of permuted phenotypes dimension  nsubjects (rows) * npermutations (cols)
            } else {
                yy = NULL
            }
            lm.f=t(apply( D, 1, function( X, y, multiplier, yy ) {
                multiplier %*% X
                f=lm.fit( X,y )
                res = c( f$df, sum(f$residuals*f$residuals))
                if ( !is.null(yy) ) {
                    res.p = apply( yy, 2, function( y.perm, f ) {
                        r = qr.resid( f$qr, y.perm )
                        return(sum(r*r))
                    }, f )
                    res = c( res, res.p )  
                } else {
                    return(res)
                }
            }, mm.use$y, mm.use$multiplier, yy ))
            mu = mean(mm.use$y)
            
        } else {
            if ( mm$perms>0 ) {
                yy = sapply( 1:mm$perms, function( i,y) { sample(y, replace=FALSE) }, mm.use$y.original )  # matrix of permuted phenotypes dimension  nsubjects (rows) * npermutations (cols)
            } else {
                yy = NULL
            }
            lm.f=t(apply( D, 1, function( X, y, yy ) {
                res = c( f$df, sum(f$residuals*f$residuals))
                if ( !is.null(yy) ) {
                    res.p = apply( yy, 2, function( y.perm, f ) {
                        r = qr.resid( f$qr, y.perm )
                        return(sum(r*r))
                    }, f )
                    res = c( res, res.p )  
                } else {
                    return(res)
                }
            }, mm.use$y.original, yy ))
            mu = mean(mm.use$y.original)
        }
        
        tss = sum(mm.use$y*mm.use$y)-mu*mu*n
        N = ncol(lm.f)
        if ( N == 2 ) { # no permutations
            rss = lm.f[,2]
            fss = tss - rss
            rdf = lm.f[,1]
            fdf = n - rdf
            F = (fss/fdf)/(rss/rdf)
            logP = -pf(F,fdf,rdf, lower.tail=FALSE, log.p=TRUE)/log(10)
            df = data.frame( bp=map$bp, logP=logP )
            data = list( phenotype=mm$p, chr=chr, df=df)
        }
        else {
            rss = lm.f[,2:N]
            perm.logP = t(apply( lm.f, 1, function( row, tss, n ) {
                rss = row[2:length(row)]
                fss = tss - rss
                rdf = row[1]
                fdf = n - rdf
                F = (fss/fdf)/(rss/rdf)
                logP = -pf(F,fdf,rdf, lower.tail=FALSE, log.p=TRUE)/log(10)
            }, tss, n ))
            logP = perm.logP[,1]
            perm.logP = perm.logP[,2:ncol(perm.logP)]
            df = data.frame( bp=map$bp, logP=logP  )
            data = list( phenotype=mm$p, chr=chr, df=df, perm.logP=perm.logP, type="haplotype")
        }
        cat("writing", file, "\n")
        save(data, file=file)
    }
}

mixed.model.pruned.allele.dosage.scan.per.chr <- function( chr, mm, out.dir=NULL, dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData", chicago=FALSE, mc.cores=1 ) {

    file=paste(out.dir, "/", mm$p, ".", chr, ".RData", sep="")
    if ( !is.null(mm$perms) ) {
        file=paste(out.dir, "/", mm$p, ".", chr, ".permute.RData", sep="")
    }
    
    if ( ! file.exists(file)) {
        cat("attempting ", file, "\n")
        dosages = load.pruned.allele.dosages.for.chr2( chr=chr, dir=dir, suffix=suffix, chicago=chicago  )
        did = dosages$ids
        map = dosages$map
        if ( mm$per.chromosome ) { # use chromosome-specific adjustments
            mm.use = mm$mm.chr[[chr]]
            cat("Using per.chromosome\n")
        }
        else {
            mm.use = mm$mm.all # use a common adjustment
        }
        use = match( mm.use$ids, did, nomatch=0)
        D = dosages$dosages
        D = D[,use]
        n = nrow(D)
        if ( mm$mixed.model == TRUE ) {
            D = t(mm.use$multiplier) %*% t(D) # dimension nsub (rows) * nsnps (cols)
            r = cor( D,  mm.use$y )
            nsub = length(mm.use$y)
        }
        else {
            D = t(D)
            r = cor( D,  mm.use$y.original )
            nsub = length(mm.use$y.original)
        }
        
        r2 = r*r
        r2 = ifelse(is.na(r2), 0, r2)
        F = (r2*(nsub-2))/(1-r2+1.0e-20)
        logP = -pf(F, 1, nsub-2, lower.tail = FALSE, log.p = TRUE)/log(10)
        df = data.frame( bp=map$bp, logP=logP )
        data = list( phenotype=mm$p, chr=chr, df=df, type="allele")

        if ( mm$perms>0 ) {
#        outfile    
            yy = sapply( 1:mm$perms, function( i,y) { sample(y, replace=FALSE) }, mm.use$y )  # matrix of permuted phenotypes dimension  nsubjects (rows) * npermutations (cols)
            
            r = cor( D, yy ) # dimension nsnps (rows) * nperms (cols) # NEEDS FIXING
            r2 = r*r
            r2 = ifelse(is.na(r2), 0, r2)
            F = (r2*(nsub-2))/(1-r2+1.0e-20)
            data$perm.logP = apply( F, 2,
                function( FF, nsub ) {
                    -pf(FF, 1, nsub-2, lower.tail = FALSE, log.p = TRUE)/log(10) # record  logP of each permutation
                }, nsub )
        }
        cat("writing", file, "\n")
        save(data, file=file)
    }
}

find.qtls <- function( logP, low.thresh ) {

    d = diff( logP ) 
    d1 = c( 0, d )
    d2 = c( d, 0 )
    w = which( d1 >= 0 & d2 <= 0 & logP>low.thresh ) # find local maxima exceeding a low initial threshold
    nw = length(logP)
    if ( length(w) ) {
        mx = t(sapply( w, function( i, nw ) {
           
            lp = logP[i]
            ok =  which(logP > lp )
            ok.right = ok[ok>i]
            wright = length(logP)
            right.edge = TRUE
            if ( length(ok.right) > 0 ) {
                wright = min(ok.right)  # nearest coordinate to the right with logP exceeding the current max
                right.edge = FALSE
            }
            ok.left = ok[ok<i]
            wleft = 1
            left.edge = TRUE
            if ( length(ok.left) > 0 ) {
                wleft = max(ok.left) # nearest coord to the left exceeding the current max
                left.edge = FALSE
            }
                                        #            cat(  "i", i, "l", wleft, "r", wright, "ok.r", ok.right,"\n")
            return( c(i, wleft,  wright, logP[i], left.edge, right.edge ) )   # left.edge, right.edge say whether neareast coord is at the edge
        }, nw ))

        return(mx)
    }
    return(NULL)
}

find.qtls.in.chr<- function( df, phenotype, chr, low.thresh, dist.thresh ) {

    mx = find.qtls( df$logP, low.thresh)

    if ( !is.null(mx)) {
        nw = nrow(mx)
        dfw = df[mx[,1],]
        mxbp = max(dfw$bp)
        dist = apply( data.frame(dfw$bp-df$bp[mx[,2]], df$bp[mx[,3]]-dfw$bp, mx[,5], mx[,6] ), 1,
            function(d) {
                s = d[3]+d[4]   # number of edge events
                if ( s == 0 ) { return(min(d[1],d[2])) }
                else if ( s == 2 ) { return(mxbp) } #   fake value to ensure QTL is called
                else if ( d[3] ) { return( d[2] ) }
                else { return(d[1]) }
            })
                                                                                                                                                   
        local.maxima = data.frame(phenotype=rep(phenotype, nw), chr=rep(chr, nw), dfw, left.bp=df$bp[mx[,2]], right.bp=df$bp[mx[,3]], dist=dist )
        local.maxima = local.maxima[local.maxima$dist >= dist.thresh,]
        return(local.maxima)
    }
    return(NULL)
}

qtl.scan <- function( dir="pruned.permuted", chrs=c(1:19,"X"), low.thresh=4, dist.thresh=3.0e6, decreasing=FALSE, limit=NULL, pheno=NULL, qtl.file="qtls.RData", all.file="qtls.all.RData", fdr.thresh=0.1, mc.cores=30 ) {
                                        # identify QTLS from genome scan data and
                                        # compute false discovery rates
                                        # To be reported, QTLs must have a logP of at least low.thresh (this shoudl be set to a low value, just to eliminate obviously non-significant loci) and must be separated by at least dist.thresh from the nearest larger QTL

    qtls = NULL
    qtl.data = hash()
    
    if ( is.null(dir) & file.exists(all.file)) {
        load(all.file)
    }
    else {
        files = sort(list.files( path=dir, pattern=".RData", full=TRUE), decreasing=decreasing)
        if ( !is.null(pheno)) files = grep( pheno, files, value=TRUE, fixed=TRUE )
        if ( !is.null(limit)  ) files = files[limit[1]:limit[2]] # limit MUST be NULL or a two-number vector 
        
        
                                        #res = lapply( files , function( f,  low.thresh, dist.thresh ) {
        res = mclapply( files , function( f,  low.thresh, dist.thresh ) {
            load(f)
            cat(f,"\n")
            qtls = find.qtls.in.chr( data$df, data$phenotype, data$chr, low.thresh, dist.thresh )
            if ( !is.null(qtls)) {
                nq = nrow(qtls)
                qtls = data.frame( phenotype=rep(data$phenotype, nq), chr=rep(data$chr, nq), perm=rep(FALSE, nq), count=rep(1,nq), qtls )
            }
            permuted.qtls = find.permuted.qtls.in.chr( data, low.thresh, dist.thresh )
            if ( !is.null(permuted.qtls)) {
                nq = nrow(permuted.qtls)
                np = ncol(data$perm.logP)
                permuted.qtls = data.frame( phenotype=rep(data$phenotype, nq), chr=rep(data$chr, nq), perm=rep(TRUE, nq), count=rep(np, nq), permuted.qtls )
            }
            res = NULL
            res = rbind( qtls, permuted.qtls ) 
            return(res)
        },  low.thresh, dist.thresh, mc.cores=mc.cores )
        res = do.call( "rbind", res )
        
        if ( !is.null(all.file)) save( res, file=all.file )
    }

    phenotypes = unique( as.character(res$phenotype ))
    all.perm = res[res$perm==TRUE,]
    all.qtls = res[res$perm==FALSE,]    
    all.perm.qtls = NULL
    all.logP.perm = all.perm$logP
    
    nperm = all.perm$count[1]
    all.nperm = length(phenotypes)*nperm # total numbber of permutations

    qtls = mclapply(  phenotypes, function( p, all.perm, all.qtls, all.logP.perm, nperm, all.nperm ) {
        logP = all.qtls$logP[all.qtls$phenotype==p]
        len.logP = length(logP)
        cat( p, "\n")

        if ( len.logP>0 ) {
            logP.perm = all.perm$logP[all.perm$phenotype==p]
            r = order(logP, decreasing=TRUE)
            fd = rep(0, length(logP))
            fdr = rep(0, length(logP))
            fd.all = rep(0, length(logP))
            fdr.all = rep(0, length(logP))
            is.qtl = rep(FALSE, length(logP))
            is.qtl.all = rep(FALSE, length(logP))
            for( i in 1:length(logP)) {
                j = r[i]
                fd[j] = sum( logP.perm>=logP[j] )/nperm # expected number of false discoveries estimated from permutations
                fdr[j] = fd[j]/i # false discovery rate estimated as the number of false disoveries made at that the given threshold divided by the number actually observed in the real data
                if ( fdr[j] > 1 ) fdr[j] = 1
                fd.all[j] = sum( all.logP.perm>=logP[j] )/all.nperm # false discoveries estimated from all permutations assuming a common distribution across phenotypes
                fdr.all[j] = fd.all[j]/i
                if ( fdr.all[j] > 1 ) fdr.all[j] = 1
                cat( i, j, logP[j], fd[j], fd[j]/i, "\n")
            }
            
            fdr.adj = fdr
            fdr.all.adj = fdr.all
            r = order(logP, decreasing=FALSE)
            w=which( fdr.adj[r] <= fdr.thresh)
            if ( length(w) > 0 ) {
                lpr = logP[r]
                logP.thresh = lpr[w[1]]
                is.qtl = ifelse( logP >= logP.thresh, TRUE, FALSE )
            }
            w=which( fdr.all.adj[r] <= fdr.thresh)
            if ( length(w) > 0 ) {
                lpr = logP[r]
                logP.thresh = lpr[w[1]]
                is.qtl.all = ifelse( logP >= logP.thresh, TRUE, FALSE )
            }
            
            
            if ( len.logP > 1 ) {
                
                
                for( i in 2:len.logP) {
                    j = r[i]
                    j1 = r[i-1]
                    if ( fd[j] == fd[j1] ) fdr.adj[j] = fdr[j1]
                    if ( fd.all[j] == fd.all[j1] ) fdr.all.adj[j] = fdr.all[j1]
                }
            }
            
            
            qtls = all.qtls[all.qtls$phenotype==p,]
            qtls$fd = fd
            qtls$fdr = fdr
            qtls$fd.all = fd.all
            qtls$fdr.all=fdr.all
            qtls$fdr.adj = fdr.adj
            qtls$fdr.all.adj = fdr.all.adj
            qtls$is.qtl = is.qtl
            qtls$is.qtl.all = is.qtl.all
            
            r = order(logP, decreasing=TRUE)
            qtls = qtls[r,]
            
            return(qtls)
        }
    },  all.perm, all.qtls, all.logP.perm, nperm, all.nperm, mc.cores=mc.cores )
    qtls = do.call( "rbind", qtls )
    
    logP = all.qtls$logP
    logP.perm = all.perm$logP
    r = order(logP, decreasing=TRUE)
    qtls$fdr.global[r] = sapply( 1:length(logP), function(i, r, logP, logP.perm, nperm) {
        fd = sum( logP.perm>=logP[r[i]] )/nperm # expected number of false discoveries estimated from permutations
        return (fd/i) # false discovery rate estimated as the number of false disoveries made at that the given threshold divided by the number actually observed in the real data
    }, r, logP, logP.perm, nperm )
    
    qtls$from.bp = rep(NA, nrow(qtls))
    qtls$to.bp = rep(NA, nrow(qtls))
    
    save(qtls, file=qtl.file)
    return(qtls)
}

qtl.rescan <- function( qtl.file, mm.dir="/Net/sparse/data/scratch/MOTT/OUTBRED/kinship/", out.dir="./RESCAN/", fdr.thresh=0.1, window=1.5e6, per.chromosome=TRUE, mixed.model=TRUE, min.maf=0.001, filter=TRUE, mc.cores=40 ) {

    if ( ! file.exists( out.dir ) ) dir.create( out.dir)
    
    mm.files = list.files( path=mm.dir, pattern=".RData", full=TRUE)
    mm.phenotypes = sub("^.+emma.", "", mm.files )
    mm.phenotypes = sub(".RData$", "", mm.phenotypes )
    
    load(qtl.file)
    qtls = qtls[qtls$fdr.adj<fdr.thresh,]
    chrs = sort(unique( qtls$chr ))
    for( chr in chrs ) {
        ff.data = rwd.load.dosages.for.chr( chr, type="allele", min.maf=min.maf, filter=filter )
        which.qtls = which(qtls$chr == chr )
        if ( length(which.qtls) > 0 ) {
            qtls.chr = qtls[which.qtls,]
            phenotypes = sort(unique(qtls.chr$phenotype))
            mclapply( phenotypes, function( p, qtls.chr, mm.files, out.dir, window, ff.data, chr ) {
                which.p = which( qtls.chr$phenotype == p )
                load(mm.files[match(p,mm.phenotypes)])
                for( w in which.p ) {
                    qtl = qtls.chr[w,]
                    bp = qtl$bp
                    out.file = paste( out.dir, "/", p, ".", chr, ".", bp, ".RData", sep="")
                    cat(paste(qtl, sep="\t"),"\n")
                    if ( !file.exists(out.file)) {
                        if ( is.numeric( qtl$from.bp ) ) {from.bp = qtl$from.bp}
                        else {from.bp = bp-window}
                        if ( is.numeric( qtl$to.bp )) { to.bp = qtl$to.bp }
                        else {to.bp = bp+window}
                        dosages = load.complete.allele.dosages.for.locus( ff.data, chr, from.bp, to.bp )
                        did = dosages$ids
                        map = dosages$map
                        if ( per.chromosome ) { # use chromosome-specific adjustments
                            mm.use = mm$mm.chr[[chr]]
                        }
                        else {
                            mm.use = mm$mm.all # use a common adjustment
                        }
                        use = match( mm.use$ids, did, nomatch=0)
                        D = dosages$dosages
                        D = D[,use]
                        n = nrow(D)
                        if ( mixed.model == TRUE ) {
                            D = t(mm.use$multiplier) %*% t(D) # dimension nsub (rows) * nsnps (cols)
                            r = cor( D,  mm.use$y )
                            nsub = length(mm.use$y)
                        }
                        else {
                            D = t(D)
                            r = cor( D,  mm.use$y.original )
                            nsub = length(mm.use$y.original)
                        }
                        
                        r2 = r*r
                        r2 = ifelse(is.na(r2), 0, r2)
                        F = (r2*(nsub-2))/(1-r2+1.0e-20)
                        logP = -pf(F, 1, nsub-2, lower.tail = FALSE, log.p = TRUE)/log(10)
                        df = data.frame( bp=map$bp, logP=logP )
                        rescan.data = list( phenotype=mm$p, chr=chr, bp=bp, from=from.bp, to=to.bp, df=df, type="allele")
                        save(rescan.data,file=out.file)
                        cat("wrote ", out.file, "\n")
                    }
                }
            }, qtls.chr, mm.files, out.dir, window, ff.data, chr, mc.cores=mc.cores )
        }
    }
}

permutation.thresholds <- function ( in.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans" ) {
    files = sort(list.files( path=in.dir, pattern=".permute.RData", full=TRUE))
                                        #    files = files[1:100]
    res = mclapply( files, function( f ) {
        load(f)
        cat(f,"\n")
        mx=apply( data$perm.logP, 2, max )
        res = list( phenotype=data$phenotype, chr=data$chr, mx=mx )
        return(res) }, mc.cores=20 )
    df = data.frame( phenotype=rep(NA,length(res)), chromosome=rep(NA, length(res)) )
    mx = matrix( 0, nrow=length(res), ncol=length(res[[1]]$mx))
    i = 1
    for( x in res ) {
        df[i,] = c(res[[i]]$phenotype, res[[i]]$chr )
        mx[i,] = res[[i]]$mx
        i = i+1
    }
  
    phen = unique( df$phenotype )
    max.logP = aggregate( as.data.frame(mx), by=list(phenotype=df$phenotype), FUN=max )
    z = max.logP[,2:ncol(max.logP)]
    browser()
        max.logP[,2:ncol(max.logP)] = apply( z, 1, sort )
        
    return(max.logP)
}

genomewide.threshold.plots <- function( max.logP, pdf.file="permutation.ecdf.pdf" ) {
    pdf(pdf.file)
    n = 0
    z= as.matrix(max.logP[,2:ncol(max.logP)])
     e = ecdf( as.numeric(z ) )
    plot(e, xlab="logP", ylab="Pr(max>x)", col="red")
    apply( z, 1, function( x ) {
        e = ecdf( x)
        plot(e, col="black", add=TRUE, do.points=FALSE)
        n = n+1
    })
    plot(e, col="red", add=TRUE)
    dev.off()
    quantiles=c( 0.5, 0.9, 0.95, 0.99)
    for ( q in quantiles ) {
        y = uniroot( function( x, q, e) { e(x)-q} , q=q, e=e, lower=0.00001, upper=10 )
        cat( q, y$root, "\n")
    }
}
    
            
        
find.permuted.qtls.in.chr <- function( data, low.thresh, dist.thresh ) {
    df = data$df
    mx = apply( data$perm.logP, 2, find.qtls, low.thresh )

    if ( !is.null(mx)) {
        mx = do.call( "rbind", mx )
        dfw = df[mx[,1],]
        nw = nrow(dfw)
        dist = apply( data.frame(dfw$bp-df$bp[mx[,2]], df$bp[mx[,3]]-dfw$bp ), 1, min )
        local.maxima = data.frame(phenotype=rep(data$phenotype, nw), chr=rep(data$chr, nw), dfw, left.bp=df$bp[mx[,2]], right.bp=df$bp[mx[,3]], dist=dist )
        local.maxima$logP = mx[,4]
        local.maxima = local.maxima[local.maxima$dist >= dist.thresh,]
        return(local.maxima)
    }
    return(NULL)
}

# GSCANDB FUNCTIONS

chromosome.lengths <- function( chrs=c(1:19,"X") , genome.build=41, chicago=FALSE) {
    mx = c()
    for( c in chrs ) {
        data=load.pruned.allele.dosages.for.chr(c, chicago=chicago);
        mx = c( mx, max(data$map$bp))
    }
    paste( "insert into chromosome(name, genome_build_id, length)  values( \"", c, "\",", genome.build, "," , mx , ");", collapse = "\n", sep="")
}

make.gscandb.files <- function( in.dir, out.dir, logP.threshold=2, population="CFW", build="38", label="allele.logP", chrs=c(1:19,"X"), dosage.dir="/Net/dense/data/scratch/nicod/Robbie_run5.5_puned_maxgen/", suffix=".prunedgen.final.maf001.0.98.RData", chicago=FALSE,  mc.cores=20 ) {


    
    if ( ! file.exists(out.dir)) dir.create( out.dir )
    script.file=paste(out.dir, "/gscandb.upload.sh", sep="")
    
    if ( !is.null(chrs)) {
        marker.file = paste(out.dir, "/upload.markers.sh",sep="");
        
        fh = file( marker.file, open="w")
        cat( "#!/bin/sh\ncd /data/www/gscandb/upload\n", file=fh)
        txt = unlist( mclapply( chrs, function( chr, dosage.dir, suffix, out.dir, chicago ) {
            pruned.dosages = load.pruned.allele.dosages.for.chr( chr=chr, dir=dosage.dir, suffix=suffix, chicago=chicago )
            marker.file = paste(out.dir, "/", "marker.", chr, ".csv", sep="")
            nm = nrow(pruned.dosages$map)
            nulls = rep("",nm)
            m.name = paste(rep(chr, nm), ":", pruned.dosages$map$bp, sep="")
            df = data.frame( name=m.name,marker_type=rep("SNP", nm),leftseq=nulls,rightseq=nulls,alias=nulls)
            write.csv(df,marker.file, quote=F, row=F)
            cat(marker.file,"\n")
            
            mapping.file = paste(out.dir, "/", "marker_mapping.", chr, ".csv", sep="")
            mapping=data.frame(marker=m.name,genome_build=rep(build,nm),chromosome=rep(chr,nm),bp_position=pruned.dosages$map$bp,strand=rep(1,nm),cm=nulls)
            write.csv(mapping,mapping.file, quote=F, row=F)
            cat(mapping.file,"\n")
            return( paste( "sh gscanLoad.sh marker=", marker.file, " marker_mapping=", mapping.file, "\n", sep="") )
        },  dosage.dir, suffix, out.dir, chicago, mc.cores=mc.cores ))
        cat( paste( txt, sep=""), file=fh )
        close( fh )
    }
    
      
   
        
    markers = NULL
    marker.mapping = NULL
    files = list.files( path=in.dir, pattern=".X.permute.RData", full=FALSE)
    phenotypes = unique(sub( ".X.permute.RData", "", files ))
    files = list.files( path=in.dir, pattern=".RData", full=TRUE)
    out.files = unlist(mclapply( phenotypes , function( phenotype, files, out.dir, logP.threshold, population, build ) {
        
        phen.files = grep( phenotype, files, value=TRUE, fixed=TRUE )
        out.file = paste( out.dir, "/", phenotype, ".gscan", sep="")
        o.file = paste( phenotype, ".gscan", sep="")
        
        df = NULL
        for( file in phen.files ) {
            load(file)
#            cat(file, "\n")
            if ( data$phenotype != phenotype ) {
                warning( "incorrect phenotype ", data$phenotype, phenotype, "\n")
            }
            else {
                m.name = paste(data$chr, ":", data$df$bp, sep="")
                df = rbind ( df, data.frame( locus=m.name, logP=data$df$logP ))
            }
        }
        fh = file( out.file, open="w")
        cat( "SCAN_RESULTS\nPHENOTYPE ", phenotype, "\nPOPULATION ", population, "\nBUILD ", build, "\nUNIT logP\nSCAN_TYPE point\nBEGIN_SCAN_DATA\nlocus ", label, "\n", sep="", file=fh)
        write.table( df, col=FALSE, file=fh, quote=FALSE, row=FALSE )
        cat("END_SCAN_DATA\n", file=fh)
        close(fh)
        cat( out.file, o.file, "\n")
        return( o.file )
    }, files, out.dir, logP.threshold, population, build, mc.cores=mc.cores ))

    script = file( script.file, open="w" )
    cat( "#!/bin/sh\ncd /data/www/gscandb/upload\n", file=script)
    nf = length(out.files)
    cat( paste( "sh gscanDelete.sh pop=", population, " pheno=", phenotypes, " label=", label, collapse="\n", sep=""), "\n", file=script)
    cat( paste( "sh gscanLoad.sh dir=", out.dir, " gscan=", out.files, " labels=", label , "; rm ", out.dir, "/", out.files, "\n", collapse="\n", sep=""), "\n", file=script)
    close(script)
  
}

simulate.confidence.intervals <- function( qtl.file="/Net/dense/data/nicod/QTLs_areas_1Dec14.txt", n.sim=10000, window=5.0e6, out.dir="./", per.chromosome=TRUE, mixed.model=TRUE, mm.dir="/Net/dense/data/nicod/Ultimate_Analysis_R5.5/norm_scaled/", dir="/Net/dense/data/outbredmice/popgen/run5.5/prune_missing/pruned", suffix=".prunedgen.final.maf001.0.98.RData" , mc.cores1=10, mc.cores2=10, reverse=FALSE, chicago=chicago) {
  #JZ: replace w/ read table so that input is read in correctly and not shifted over by 1 col  
  #qtls = read.delim(qtl.file, stringsAsFactors=FALSE)
    qtls = read.table(qtl.file, header=T)
    print(qtls)
    
    if ( ! file.exists( out.dir ) ) dir.create( out.dir)
    
    chrs = sort(unique(qtls$chr) )
    if ( reverse ) chrs = rev(chrs)
    for( chr in chrs ) {
        
        qtls.chr = qtls[qtls$chr==chr,]
        files=paste(out.dir, "/", qtls.chr$phenotype, ".", chr, ".", qtls.chr$bp, ".CI.RData", sep="")
        file.existance = sapply( files, file.exists )
        qtls.chr = qtls.chr[!file.existance,]

        if ( nrow(qtls.chr) > 0 ) {
             dosages = load.pruned.allele.dosages.for.chr2( chr=chr, dir=dir, suffix=suffix, chicago=chicago  )
            mclapply( 1:nrow(qtls.chr), function( i, qtls.chr, out.dir, mm.dir, chr, dosages, n.sim, window, mc.cores2=mc.cores2) {
                p = qtls.chr$phenotype[i]
    
                file=paste(out.dir, "/", p, ".", chr, ".", qtls.chr$bp[i], ".CI.RData", sep="")
                if ( ! file.exists(file)) {
                    mm.file = paste(mm.dir, "/emma.", p, ".RData", sep="")
                    cat(mm.file, "\n")
                    load(mm.file)
                    mm$per.chromosome = per.chromosome
                    mm$mixed.model = mixed.model
                    
                    cat("attempting ", file, "\n")
                    sim = simulate.qtl.confidence.interval( chr, qtls.chr$bp[i], mm, dosages, n.sim=n.sim, window=window, out.dir=out.dir, mc.cores=mc.cores2 )
                    sim.data = list( qtl=qtls.chr[i,], n.sim=n.sim, window=window, sim=sim )
                    save(sim.data,file=file)
                }
                else {
                    cat( file, " exists\n")
                }
            },  qtls.chr, out.dir, mm.dir, chr, dosages, n.sim, window, mc.cores2=mc.cores2, mc.cores=mc.cores1 )
        }
    }
}
        
simulate.qtl.confidence.interval <- function( chr, bp, mm, dosages, n.sim=1000, window=5.0e6, out.dir="./",  mc.cores=48 ) {
    print("simulate.qtl.confidence.interval")
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
    use = match( mm.use$ids, did, nomatch=0)
    D = dosages$dosages
    D = D[,use]
    
    #JZ: find SNPs within the window (default = 5MB) and subset the dosages and map files
    snps.subset2 = which( abs(map$bp -bp) < window ) # outer subset defining the search window
    D.qtl = D[snps.subset2,]
    map.qtl = map[snps.subset2,]
    bp.idx = which( map.qtl$bp == bp )
    if ( length(bp.idx) == 0 ) {
        warning(paste("ERROR  coordinate" ,toString(chr), toString(bp) , "not found\n", sep=" "))
        return(NULL);
    }
    
    #JZ: n is number of SNPs within the window
    n = nrow(D.qtl)
    
    #JZ: what is the point of this?  compute the residuals?
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
    r2 = ifelse(is.na(r2), 0, r2)
    
    phen.r2 = r2[bp.idx]
    phen = lm.snp$y
    phen.resid = resid(lm.snp)
    
    #JZ: find all snps that are within window/2 and correlated with phenotype.  
    snps.subset1 = which( abs(map.qtl$bp -bp) < window/2 & r2>0 )
    print(dim(snps.subset1))
    
    #JZ: Select 1000 random snps to be causal and perform simulations with. 
    print("Selecting 1000 snps...")
    snps.sample = sort(sample(snps.subset1, n.sim, replace=TRUE)) # the snp sites to simulate from
    
    #JZ: compute summary statistics for 1000 snps
    var.snp = apply( D.qtl, 2, var, na.rm=TRUE ) # variability of the dosages at each snp
    snps.beta = r[bp.idx]*sqrt(var(phen))/sqrt(var.snp)
    snps.sample = c( bp.idx, snps.sample )
    
    print("Apply simulation to 1000 snps...")
    #JZ: apply simulation framework to 1000 snps.  Assume that snp is causal; 
    sims = mclapply(snps.sample, function( snp, snps.beta, var.snp, phen.resid, D.qtl, map.qtl, target.r2 ) {
        #            y = sample(phen.resid) + snps.beta[snp]*D.qtl[snp,]
        #JZ: simulate phenotype assuming that SNP is causal.  What is phen.resid?  why is it sampled but D.qtl is not? 
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
    
    print("concat results...")
    sims = do.call( "rbind", sims )
    colnames(sims) = c("snp", "bp", "snp.logP", "max.bp", "max.logP", "bp.delta", "snp.beta", "snp.var")
    return(t(sims))
    
}

evaluate.confidence.intervals <- function( dir = "./", pdf.file="CI.pdf", snp.window=1000) {
    files = sort(list.files( path=dir, pattern=".RData", full=TRUE))
#    files= files[1:100]
    results = mclapply( files, function( file, snp.window ) {
        load(file) # loads sim.data
        sim = t(sim.data$sim)
        qtl = sim.data$qtl
        snp.logP = mean( sim[,"snp.logP"])
        max.logP = mean( sim[,"max.logP"])
        delta = sort( sim[,"bp.delta"] )
        n = length(delta)
        region.CI.25 = delta[n*0.25]
        region.CI.75 = delta[n*0.75]
        region.CI.05 = delta[n*0.05]
        region.CI.95 = delta[n*0.95]
        region.CI.50 = delta[n*0.5]
        region.mean.bp = mean(abs(delta))

        r = order( abs( qtl$bp- sim[,"max.bp"] ) )
        w = r[1:snp.window]
        snp.n = length(w)
        delta.snp =  delta = sort( sim[w,"bp.delta"] )
        snp.CI.25 = delta.snp[snp.n*0.25]
        snp.CI.75 = delta.snp[snp.n*0.75]
        snp.CI.05 = delta.snp[snp.n*0.05]
        snp.CI.95 = delta.snp[snp.n*0.95]
        snp.CI.50 = delta.snp[snp.n*0.5]
        snp.mean.bp = mean(abs(delta.snp))

        res = list( qtl=qtl, CI = c( snp.CI.05, snp.CI.25, snp.CI.50, snp.CI.75, snp.CI.95, snp.mean.bp, snp.n, region.CI.05, region.CI.25, region.CI.50, region.CI.75, region.CI.95, region.mean.bp),  sim=sim )
        return(res)
    }, snp.window, mc.cores=10 )
    pdf( pdf.file )
    par(mfrow=c(2,2))

    df = NULL
    
    idx=1
    for( r in results ) {
        idx = idx+1
        if ( idx == -1 ) break 
        qtl = r$qtl
        sim = r$sim
        CI = r$CI
        w90=CI[5]-CI[1]
        w50=CI[4]-CI[2]
        rw = order( abs( qtl$bp- sim[,"max.bp"] ) )
        w = rw[1:snp.window]

        from.bp = min(sim[w,"max.bp"],na.rm=TRUE)
        to.bp = max(sim[w,"max.bp"],na.rm=TRUE)
        cat( length(w), unlist(qtl), "CI", CI, "from", from.bp, "to", to.bp,  "\n")

        sim.w = sim[w,]
        order.w = order(sim.w[,"bp"]- sim.w[,"max.bp"])
        sim.w = sim.w[order.w,]
        coords = sim.w[,"bp"]

        do.it = FALSE
        
        if ( do.it ) {
            delta.50 = floor(length(w)*0.5)
            n.50 = delta.50 -1
            m.50 = 1:n.50
            len.50 = coords[m.50+delta.50]-coords[m.50]
            w.50 = which.min( len.50 )
            start.50 = coords[w.50]
            end.50 = coords[w.50+delta.50]
            start.50.bp = sim.w[w.50,"max.bp"]
            end.50.bp = sim.w[w.50+delta.50,"max.bp"]
            
            delta.10 = floor(length(w)*0.1)
            delta.90 = floor(length(w)*0.90)
            n.10 = delta.10 -1
            m.10 = 1:n.10
            len.90 = coords[m.10+delta.90]-coords[m.10]
            w.90 = which.min( len.90 )
            start.90 = coords[w.90]
            end.90 = coords[w.90+delta.90]
            start.90.bp = sim.w[w.90,"max.bp"]
            end.90.bp = sim.w[w.90+delta.90,"max.bp"]
        }
        delta = sim[w,"bp.delta"]
        plot( density(delta/1.0e6),main = paste(qtl$phenotype, qtl$chr, qtl$bp, sprintf( "%.2f", qtl$logP)), sub=paste("50%:", w50/1.0e6, "90%:", w90/1.0e6), xlab="Error, Mb", ylab="frequency",cex=0.5, cex.main=0.8)
        
        
        abline(v=c(CI[1],CI[5])/1.0e6, col="red")
        abline(v=c(CI[2],CI[4])/1.0e6, col="green")
#        abline(v=c(start.50-start.50.bp,end.50-end.50.bp)/1.0e6, col="pink")
#        abline(v=c(start.90-start.90.bp,end.90-end.90.bp)/1.0e6, col="orange")
        sim = r$sim
        range.bp = range( sim[,"bp"], sim[,"max.bp"])/1.0e6
#        plot( sim[,"bp"]/1.0e6, sim[,"max.bp"]/1.0e6, pch=20, xlim=range.bp, ylim=range.bp, main=paste(qtl$phenotype, qtl$chr, qtl$bp, sprintf( "%.2f", qtl$logP)), cex=0.5,cex.main=0.8, xlab="True location Mb", ylab="predicted Mb" )
        high.var = sim[,"snp.var"] > 0.05
                                        #        points( sim[high.var,"bp"]/1.0e6, sim[high.var,"max.bp"]/1.0e6, pch=20,col="red",cex=0.4)
        plot( c(1),c(1), xlim=range.bp, ylim=range.bp, main=paste(qtl$phenotype, qtl$chr, qtl$bp, sprintf( "%.2f", qtl$logP)), cex=0.5,cex.main=0.8, xlab="True location Mb", ylab="predicted Mb", t="n")
        rect( xleft=(qtl$bp+CI[1])/1.0e6,ybottom=range.bp[1],xright=(qtl$bp+CI[5])/1.0e6,ytop=range.bp[2],col="red",border=NA)
        rect( xleft=(qtl$bp+CI[2])/1.0e6,ybottom=range.bp[1],xright=(qtl$bp+CI[4])/1.0e6,ytop=range.bp[2],col="green",border=NA)
#        rect( xleft=start.90/1.0e6,ybottom=range.bp[1],xright=end.90/1.0e6,ytop=range.bp[2],col="orange",border=NA)
                                        #        rect( xleft=start.50/1.0e6,ybottom=range.bp[1],xright=end.50/1.0e6,ytop=range.bp[2],col="pink",border=NA)
        abline(h=qtl$bp/1.0e6, col="black")
        smoothScatter( sim[,"bp"]/1.0e6, sim[,"max.bp"]/1.0e6,colramp = alpharamp("white", blues9),add=TRUE )
       
#        abline(h=from.bp/1.0e6, col="blue")
#        abline(h=to.bp/1.0e6, col="blue")
#        abline(v=c(qtl$bp+CI[1],qtl$bp+CI[5])/1.0e6, col="red")
#        abline(v=c(qtl$bp+CI[2],qtl$bp+CI[4])/1.0e6, col="green")
#        abline(v=c(start.50,end.50)/1.0e6, col="pink")
#        abline(v=c(start.90,end.90)/1.0e6, col="orange")
        df = rbind( df, unlist(c(qtl[1,1:5], CI[1:5]+qtl$bp, CI[6] )))
    }
    dev.off()
    df = data.frame(df)
    names(df) = c("phenotype", "chr" ,      "bp"    ,    "logP"  ,    "fdr.adj" ,     "CI.05"   ,       "CI.25"  ,        "CI.50",
 "CI.75"     ,     "CI.95"       ,   "mean.error")
    write.table(df, "CI.txt", quote=F, row=F, sep="\t")
    
} 
        

evaluate.confidence.intervals.logP <- function( qtl.dir = "/Net/sparse/data/scratch/MOTT/OUTBRED/CI.03032015/", scan.dir ="/Net/dense/data/nicod/Maxgen_Analysis_R5.5/norm_scaled249_scans/", pdf.file="CI.pdf", quantiles=c(0.25,0.5,0.75,0.9), dist.thresh=3.0e6, CI.file = "CI.txt", fine.dir= "/Net/dense/data/scratch/nicod/Robbie_run5.5_qtl.rescan_norm_scaled249_peaks/", fine.suffix="_Robbie_Run5.5_Maxgen_13Feb15_peak_annot.txt" , mc.cores=10) {
    files = sort(list.files( path=qtl.dir, pattern=".RData", full=TRUE))
#                                        files= files[1:100]
   
#     files = files[grep("Haem.RBC",files)]
    
    results = mclapply( files, function( file, dist.thresh, scan.dir ) {
        load(file) # loads sim.data
        cat("loading ", file, "\n")
        sim = t(sim.data$sim)
        qtl = sim.data$qtl
        delta.logP = abs(sim[,"max.logP"])-abs(sim[,"snp.logP"])
        delta.logP = ifelse( delta.logP < 0 , 0, delta.logP )
        delta.bp = abs(qtl$bp- sim[,"max.bp"])
        delta.bp = ifelse ( delta.logP == 0 , 0, delta.bp)
        r = order(delta.logP)
        delta.logP.r = delta.logP[r]
        q = round(length(r)*quantiles)
        CI.logP = delta.logP.r[q]

        
        scan.file = paste(scan.dir, "/", qtl$phenotype, ".", qtl$chr, ".permute.RData", sep="")
        if ( file.exists(scan.file)) {
            cat("loading ", scan.file, "\n")
            load(scan.file)
#            if ( qtl$phenotype == "Haem.RBC" & qtl$bp == 13940848 ) browser()
            df = data$df
            bp = qtl$bp
            lp = df$logP
            CI = matrix(NA, nrow=length(CI.logP), ncol=17)
            colnames(CI) = c("bp", "quantile", "delta.logP", "from.bp", "to.bp", "width.bp", "n.snps", "fine.bp", "fine.logP", "fine.from.bp", "fine.to.bp", "fine.width.bp", "fine.n.snps", "fine.CI.logP", "fine.CI.bp", "fine.CI.from.bp", "fine.CI.to.bp")
            CI[,1] = rep(bp,nrow(CI))
            CI[,2] = quantiles
            CI[,3] = CI.logP
            for( i in 1:length(CI.logP)) {
                w = which( (qtl$logP - df$logP <= CI.logP[i] ) & ( abs(df$bp-bp ) < dist.thresh ))
                if ( length(w) > 0 ) {
                    CI[i,4] = min(df$bp[w])
                    CI[i,5] = max(df$bp[w])
                    CI[i,6] = CI[i,5]-CI[i,4]+1
                    CI[i,7] = length(w)
                }
            }
            if ( !is.null(fine.dir)) { # Refine CIs using fine-mapping data with more SNPs
                fine.file = paste( fine.dir, "/", qtl$phenotype,"_chr", qtl$chr, ".", qtl$bp, fine.suffix, sep="")
                cat(fine.file,"\n")
               
                
                if ( file.exists(fine.file)) {
                    cat("loaded", fine.file,"\n")
                    y = read.delim(fine.file)
                    df = y[!is.na(y$logP),]
                    if ( qtl$chr == "X" ) {
                        subw = which(  df$score>0.4)
                    } else {
                        subw = which( df$score>0.4 & df$hwe > 1.0e-6) # quality filter
                    }
                    df = df[subw,]
                    w.max.logP = which.max( df$logP)
                    bp = df$bp[w.max.logP]
                    max.logP = df$logP[w.max.logP]
                    for( i in 1:length(CI.logP)) {
                        w = which( (max.logP - df$logP <= CI.logP[i] ) & ( abs(df$bp-bp ) < dist.thresh )) # max logp within dist.thresh
                        if ( length(w) > 0 ) {
                            CI[i,8] = bp
                            CI[i,9] = max.logP
                            CI[i,10] = min(df$bp[w])
                            CI[i,11] = max(df$bp[w])
                            CI[i,12] = CI[i,5]-CI[i,4]+1
                            CI[i,13] = length(w)
                        }

                        if ( qtl$chr == "X" ) {
                            ww = which(  df$bp >= CI[i,4] & df$bp <= CI[i,5] ) # max logp wihin the pre-existing CI
                        } else {
                            ww = which(  df$bp >= CI[i,4] & df$bp <= CI[i,5] ) # max logp wihin the pre-existing CI
                        }
                        if ( length(ww) > 0 ) {
                            dfww = df[ww,]
                            dfw = df[w,]
                            ww.max = which.max(dfww$logP)
                            max.bp = dfww$bp[ww.max]
                            CI[i,14] = dfww$logP[ww.max]
                            CI[i,15] = max.bp
                            w.max = which( dfw$bp == max.bp ) 
                            u = which( dfw$logP[w.max] - dfw$logP < CI.logP[i] ) # define CI allowing use of fine SNPs outside of pruned CI
                            u = which( CI[i,14] - dfw$logP < CI.logP[i] ) # define CI allowing use of fine SNPs outside of pruned CI
                            if( length(u) > 0 ) {
                                CI[i,16] = min( dfw$bp[u] )
                                CI[i,17] = max( dfw$bp[u] )
                            }
                        }
                    }
                } else {
                    cat ("ERROR", fine.file, "not found\n")
                }
            }
            
            CI.df = data.frame( phenotype = rep(qtl$phenotype,nrow(CI)),  logP=rep(qtl$logP,nrow(CI)), chr=rep(qtl$chr,nrow(CI)), CI)
            names(CI.df) = c( "phenotype",  "logP", "chr",colnames(CI))
        res = list( qtl=qtl, CI = CI.df, sim=sim)
            return(res)
        }
        else { return(NULL) }
    }, dist.thresh, scan.dir, mc.cores=mc.cores)

    CI = NULL

    for( r in results ) {
       if ( !is.null(r) )  CI = rbind( CI, r$CI )
    }
    write.table(CI, file=CI.file, quote=F, row=F, sep="\t")
} 
        
    
alpharamp<-function(c1,c2, alpha=128) {
    stopifnot(alpha>=0 & alpha<=256);
    function(n) paste(colorRampPalette(c(c1,c2))(n), format(as.hexmode(alpha), upper.case=T), sep="")}
