library(ff)
#library(emma)
library(parallel)
library(hash)
source("scripts/emma.R")
#JZ: can't download this library??
#library("RevoUtilsMath")
#setMKLthreads(1)

# read phenotypes from file


#  phenotype.file="/Net/dense/data/nicod/Harwell_all_measures_resids_no_duplicates_21Nov13.txt

read.phenotypes <- function( phenotype.file="/Net/dense/data/nicod/Harwell_all_measures_resids_Current_All_10Apr14.txt",normalise=FALSE) {
  phen=read.delim(phenotype.file);
  pnam = names(phen)
  resid = grep(".resid", pnam, fixed=TRUE, value = FALSE)
  ppnam = sub(".resid", "", pnam, fixed=TRUE )
  nam  = phen$Animal_ID
  pheno = hash()
  for( i in resid ) {
    y = phen[,i]
    names(y) = nam
    y = y[!is.na(y)]
    if ( normalise ) { #quantile normalise
      r = rank(y)
      y = qnorm(r/((length(y)+1)))
    }
    pheno[[ppnam[i]]] = y
  }
  return(pheno)
}

rwd.load.dosages <- function( type="allele", chrs=c(1:19)) {

    data = mclapply( chrs, rwd.load.dosages.for.chr, type, mc.cores=10 )
    names(data) = chrs
    
    map = NULL
    dosages = list()
    for ( chr in chrs ) {
        map = rbind(map,data.frame( chr=data[[chr]]$posSubset[,1], bp=as.numeric(data[[chr]]$posSubset[,2])))
    }
    ids = data[[1]]$ids

    return( list( chrs=chrs, ids=ids, map=map, dosages=data, type=type))
}


    
rwd.load.dosages.for.chr <- function( chr, type="allele", filter=TRUE, min.maf=0.05, min.infoCorrect=0.4, version="5.6.OF", hap.subset=0.05, range=NULL ) { # version="5.5"


    d =loadAncestralHaplotypeData(runVersion=version,chr)
    ids = sub("_recal.bam$", "", d$nameList, perl=TRUE)
    ids = sub("___", "/", ids, perl=TRUE)
    ids = sub("_recal.reheadered.bam$", "", ids, perl=TRUE )
    
    if ( type == "haplotype" ) {
        N=d$N # number of samples                                                                                                             
        T=d$T # number of SNPs on chromosome                                                                                                  
        K=4 # number of ancestral haplotypes
        map = data.frame( chr=d$pos[,1], bp=as.numeric(d$pos[,2]))
        use=seq(1,T,by=round(1/hap.subset)) # determine which to use                
        map = map[use,]
        dd=array(0,c(length(use),N,K)) # pre-declare                                                                                          
        off=0 # offset variable                                                                                                               
        for(i in 1:3) {
                                        # set dosages file                                                                                                                  
            if(i==1) {
                dosages=d$hapProbsA.ff
                if(is.null(dosages))
                    dosages=d$hapProbs.ff
            }
            if(i==2) dosages=d$hapProbsB.ff
            if(i==3) dosages=d$hapProbsC.ff
            if(!is.null(dosages)) {
                d.dim = dim(dosages)
                                        # determine where these SNPs are in chromosome                                                                                    
                local=seq(1,d.dim[1],1) + off
                                        # choose the right ones                                                                                                           
                tmp.use=is.na(match(local,use))==FALSE
                                        # now stick these SNPs in the pre-allocated array                                                                                 
                dd[is.na(match(use,local))==FALSE,,]= dosages[tmp.use,,]
                off = off + d.dim[1]
            }
        }
        
        tot = length(use)
        intervals = nrow(dd)
        cat("chr ", chr, "intervals:", intervals,  "total:", tot, tot/(N+1.0e-10), "\n")
        return( list( chr=chr, ids=ids, intervals=intervals, map=map, dosages=dd, type=type, hap.subset=hap.subset))
        
    } else if ( type == "allele") {
        dosages = d$dosages.ff
        
        infoCorrect = d$infoCorrect
        hwe = d$hwe
        af = d$af
        maf = ifelse( af > 0.5, 1-af, af )
        tot = nrow(dosages)
        map = data.frame( chr=d$pos[,1], bp=as.numeric(d$pos[,2]))
        
        if ( filter==TRUE ) {
            if ( chr != "X" ) {
                hwe.use = hwe > 1.0e-6
            } else {
                hwe.use = rep( TRUE, length(hwe) )
            }
#             hwe.use = hwe < 0.999
            info.use = infoCorrect>min.infoCorrect
            info.use = ifelse(is.na(info.use), FALSE, info.use)
            maf.use = maf > min.maf # maf.use > 0

            if ( !is.null(range)) {
                from.bp = range[1]
                to.bp = range[2]
                range.use = map$bp >= from.bp & map.bp <= to.bp
                use = which( hwe.use & info.use & maf.use & range.use )
            } else {
                use = which( hwe.use & info.use & maf.use)
            }
            
            tot = length(use)
            N = nrow(dosages)
            
            dosages = dosages[use,]
            map = map[use,]
            infoCorrect = infoCorrect[use]
            hwe = hwe[use]
            maf = maf[use]
            cat("chr ", chr, "version:", version, "snps:", N,  "hwe:", sum(hwe.use), "info:", sum(info.use, na.rm=TRUE), "maf:", sum(maf.use), "total:", tot, tot/(N+1.0e-10), "\n")
            return( list( chr=chr, ids=ids, Nsnps=tot, map=map, dosages=dosages, type=type, infoCorrect=infoCorrect,hwe=hwe, maf=maf, range=range, filter=filter))
            
        } else {
            tot = nrow(dosages)
            return( list( chr=chr, ids=ids, Nsnps=tot, map=map, dosages=dosages, type=type, infoCorrect=NULL,hwe=NULL, maf=NULL, range=range, filter=filter))
        }
    } else if ( type == "genotype") {
        dosages = d$genProbs.ff
    }
    else {
        stop( "unknown type ", type, "\n")
    }
    
    return(NULL)
}

# #JZ: Not used in either script?
# ### the following script loads the necessary data for a chromosome given run version, chr
# loadAncestralHaplotypeData=function(runVersion,ch)
# {
#   ### set output directory
#   outputdir=paste("/Net/dense/data/outbredmice/imputation/ancestral_haps/run",runVersion,"/",sep="")
#   ### load basic things from RData files
#   load(file=paste(outputdir,"RData/EMparameters.chr",ch,".RData",sep="")) # things like N, T
#   load(file=paste(outputdir,"RData/EMnecessary.chr",ch,".RData",sep="")) # components of the EM algorithm
#   load(file=paste(outputdir,"RData/info.chr",ch,".RData",sep="")) # info
#   load(file=paste(outputdir,"RData/pos.chr",ch,".RData",sep="")) # chrom and bp positions
#   load(file=paste(outputdir,"RData/hwe.chr",ch,".RData",sep="")) # hwe
#   load(file=paste(outputdir,"RData/nameList.chr",ch,".RData",sep=""))
#   ### get flat file references - dosages here
#   dosages.ff=ff(filename=paste(outputdir,"RData/dosages.ff.chr",ch,sep=""),vmode="double",dim=c(T,N),readonly=TRUE)
#   ### determine if to split based on whether hapProbs was too big to save as one file
#   splitNum=ceiling(as.numeric(T)*as.numeric(N*K)/2147483647)
#   #as on dense .Machine$integer.max = 2147483647
#   split=splitNum>1 # old variable
#   ### now get the hapProbs and genProbs files files
#   split=is.na(file.info(filename=paste(outputdir,"RData/hapProbs.ff.chr",ch,sep=""))$size)
#   if(splitNum==1) # small chromosome
#   {
#     hapProbs.ff=ff(filename=paste(outputdir,"RData/hapProbs.ff.chr",ch,sep=""),vmode="double",dim=c(T,N,K),readonly=TRUE)
#     genProbs.ff=ff(filename=paste(outputdir,"RData/genProbs.ff.chr",ch,sep=""),vmode="double",dim=c(T,N,3),readonly=TRUE)
#   }
#   if(splitNum==2) # medium chromosome
#   {
#     Ta=as.integer(ceiling(T/2))
#     Tb=as.integer(T-ceiling(T/2))
#     hapProbsA.ff=ff(filename=paste(outputdir,"RData/hapProbsA.ff.chr",ch,sep=""),vmode="double",dim=c(Ta,N,K),readonly=TRUE)
#     hapProbsB.ff=ff(filename=paste(outputdir,"RData/hapProbsB.ff.chr",ch,sep=""),vmode="double",dim=c(Tb,N,K),readonly=TRUE)
#     genProbsA.ff=ff(filename=paste(outputdir,"RData/genProbsA.ff.chr",ch,sep=""),vmode="double",dim=c(Ta,N,3),readonly=TRUE)
#     genProbsB.ff=ff(filename=paste(outputdir,"RData/genProbsB.ff.chr",ch,sep=""),vmode="double",dim=c(Tb,N,3),readonly=TRUE)
#   }
#   if(splitNum==3) # big chromosome
#   {
#     Ta=as.integer(ceiling(T/3))
#     Tb=as.integer(ceiling(T/3))  
#     Tc=as.integer(T-2*ceiling(T/3))
#     # hapProbs
#     hapProbsA.ff=ff(filename=paste(outputdir,"RData/hapProbsA.ff.chr",ch,sep=""),vmode="double",dim=c(Ta,N,K),readonly=TRUE)
#     hapProbsB.ff=ff(filename=paste(outputdir,"RData/hapProbsB.ff.chr",ch,sep=""),vmode="double",dim=c(Tb,N,K),readonly=TRUE)
#     hapProbsC.ff=ff(filename=paste(outputdir,"RData/hapProbsC.ff.chr",ch,sep=""),vmode="double",dim=c(Tc,N,K),readonly=TRUE)
#     # genProbs
#     genProbsA.ff=ff(filename=paste(outputdir,"RData/genProbsA.ff.chr",ch,sep=""),vmode="double",dim=c(Ta,N,3),readonly=TRUE)
#     genProbsB.ff=ff(filename=paste(outputdir,"RData/genProbsB.ff.chr",ch,sep=""),vmode="double",dim=c(Tb,N,3),readonly=TRUE)
#     genProbsC.ff=ff(filename=paste(outputdir,"RData/genProbsC.ff.chr",ch,sep=""),vmode="double",dim=c(Tc,N,3),readonly=TRUE)
#   }
#     ### now return variables
#   if(splitNum==1)
#     return(list(dosages.ff=dosages.ff,hapProbs.ff=hapProbs.ff,genProbs.ff=genProbs.ff,pos=pos,L=L,gen=gen,nameList=nameList,T=T,N=N,K=K, infoCorrect=infoCorrect,hwe=hwe, af=estimatedAlleleFrequency, splitNum=1))
#   if(splitNum==2)
#     return(list(dosages.ff=dosages.ff,hapProbsA.ff=hapProbsA.ff,hapProbsB.ff=hapProbsB.ff,genProbsA.ff=genProbsA.ff,genProbsB.ff=genProbsB.ff,pos=pos,L=L,gen=gen,nameList=nameList,T=T,N=N,K=K,Ta=Ta,Tb=Tb, infoCorrect=infoCorrect,hwe=hwe, af=estimatedAlleleFrequency, splitNum=2 ))
#   if(splitNum==3)
#     return(list(dosages.ff=dosages.ff,hapProbsA.ff=hapProbsA.ff,hapProbsB.ff=hapProbsB.ff,hapProbsC.ff=hapProbsC.ff,genProbsA.ff=genProbsA.ff,genProbsB.ff=genProbsB.ff,genProbsC.ff=genProbsC.ff,pos=pos,L=L,gen=gen,nameList=nameList,T=T,N=N,K=K,Ta=Ta,Tb=Tb,Tc=Tc, infoCorrect=infoCorrect,hwe=hwe, af=estimatedAlleleFrequency, splitNum=3))
# }
# 

make.rwd.kinship.matrix <- function (chrs=1:19, type="allele", downsample=0.05, filter=TRUE, file.prefix=paste("rwd.kinship", downsample, type, filter, sep="."), scale.dosages=TRUE) {

    K.chr = list()
    K.chr = mclapply( chrs, function( chr, type, downsample, filter ) {
                                        #    K = lapply( chrs, function( chr, type, downsample, filter ) {
        rwd.dosages = rwd.load.dosages.for.chr( type=type, chr=chr, filter=filter )
        D = rwd.dosages$dosages
        n = nrow(D)
        
        if ( downsample < 1 & downsample > 0 ) {
            size = round(downsample*n)
            use = sample(n, size=size, replace=FALSE )
        }
        D = D[use,]
        if ( scale.dosages==TRUE ) {
            tD = t(D)
            tD = scale(tD)
            D = t(tD)
            cat("rescaled dosages\n")
        }
        cat("chr ", chr, n, "downsampled:", size, size/(n+1.0e-10), "\n")
        
        k = cor(D)
        rownames(k) = as.character(rwd.dosages$ids)
        colnames(k) = as.character(rwd.dosages$ids)
        return( list( K=k, n=size))
    }, type, downsample, filter, mc.cores=20)
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
    kinship.matrices = list( K.all = K.all, K.n=K.n, K.chr=K.chr, type=type, filter=filter ) # Note output has been changed to give per chromosome files as well as total kinship
    save( kinship.matrices, file=file)
    cat("wrote kinship matrix to", file, "\n")
}

    

# Create GCTA input files from the kinship matrix

kinship2gcta <- function( K, file.prefix , subjects=NULL) {

  if ( is.null(subjects))
    subjects = colnames(K)
  if ( is.null(subjects))
    subjects = 1:nrow(K)

  w =  subjects != "Unknown"
  K = K[w,w]
  subjects = subjects[w]
  
  nr = nrow(K)
  nc = ncol(K)
  if ( nr == nc ) {
    U = upper.tri(K,diag=TRUE)
    w = which(U,arr.ind=TRUE)
    M = cbind(w[,2],w[,1],nc,K[U])
    con = gzfile(paste(file.prefix, ".grm.gz", sep=""),"w")
    write.table(M,file=con, quote=FALSE, row=FALSE, col=FALSE, sep="\t")
    close(con)
    write.table(cbind(rep("HS",length(subjects)),subjects),file=paste(file.prefix, ".grm.id",sep=""), quote=FALSE, row=FALSE, col=FALSE, sep="\t")
  }
}

 
#Mixed Models

mixed.model.heritabilities <- function( kinship.matrices, phenotypes=NULL, phenotype.file="/Net/dense/data/nicod/Harwell_all_measures_resids_Current_All_10Apr14.txt", out.dir=NULL, normalise=FALSE, out.file="heritabilities.txt", mc.cores=40 ) { # , mc.cores2=2 ) { # can't have nested mcapply() statements!


    K.all = NULL
    K.n = 0
    K.chr = kinship.matrices$K.chr
    for( chr in names(K.chr)) {
        if ( is.null(K.all)) {
            K.all = K.chr[[chr]]$K * K.chr[[chr]]$n
            K.n = K.chr[[chr]]$n
            #JZ add colnames/rownames so that intersection of samples in mixed.model.variances doesn't break
            colnames(K.all) <- colnames(K.chr[[chr]]$K)
            rownames(K.all) <- rownames(K.chr[[chr]]$K)
        }
        else {
			print(K.chr[[chr]]$n)
            K.all = K.all + K.chr[[chr]]$K * K.chr[[chr]]$n
            K.n = K.n + K.chr[[chr]]$n
        }
    }
    
    for( chr in names(K.chr)) {
        K = K.chr[[chr]]$K *  K.chr[[chr]]$n
        n = K.chr[[chr]]$n
        KK = (K.all-K)/(K.n - n)
        K.chr[[chr]]$K = (K.all - K)/(K.n-n)
    }
    kinship.matrices$K.chr = K.chr

	#JZ: added this line b/c I think it should be here
	K.all = K.all/K.n
#   	print(dim(K.all))
 
    phen = read.phenotypes(phenotype.file=phenotype.file, normalise=normalise)
   
    if ( is.null(phenotypes) )
        phenotypes = names(phen)

    if ( !is.null ( out.dir ) ) {
        if ( ! file.exists(out.dir)) 
            dir.create( out.dir )
    }
    else
        out.dir="./"

    df=mclapply( phenotypes , function( p, phen, kinship.matrices, out.dir ) {
        if ( !is.null( phen[[p]] )) {
            out.file=paste(out.dir, "/emma.", p, ".RData", sep="")
            cat(out.file, "\n")
            if ( !file.exists(out.file)) {
				print("starting mixed.model.variance()")
                mm = list()
                #JZ replace kinship.matrices$K.all with K.all
                #mm$mm.all = mixed.model.variance( phen[[p]], kinship.matrices$K.all )
                mm$mm.all = mixed.model.variance( phen[[p]], K.all )
                mm$p = p
                cat( p, mm$vg, mm$ve, mm$vg/(mm$vg+ mm$ve+1.0e-10),  "\n")
                mm.chr = list()
                pp = phen[[p]]
                print("here")
          
                mm.chr = list()
                for( chr in  names(kinship.matrices$K.chr) ) {
                    cat("chr", chr, p, "\n")
                    mm.chr[[chr]] = mixed.model.variance( phen[[p]], kinship.matrices$K.chr[[chr]]$K )
                }
               
                names(mm.chr) = names(kinship.matrices$K.chr)
                v.chr = NULL
                for( chr in names(kinship.matrices$K.chr) ) {
                    v.chr = rbind( v.chr, c( mm.chr[[chr]] $vg, mm.chr[[chr]]$ve, mm.chr[[chr]]$vg/(mm.chr[[chr]]$vg+ mm.chr[[chr]]$ve+1.0e-10)))
                }
                v.chr = data.frame(v.chr)
                names(v.chr) = c( "vg", "ve", "H2")
                mm$v.chr = data.frame( phen=rep( p, nrow(v.chr)), chr=names(kinship.matrices$K.chr), v.chr )
                print(mm$v.chr)
                mm$mm.chr = mm.chr
                save(mm, file=out.file)
                return( c(p, mm$mm.all$vg, mm$mm.all$ve, mm$mm.all$vg/(mm$mm.all$vg+ mm$mm.all$ve+1.0e-10)))
            }
        }
        return( NULL )
    #}, phen, kinship.matrices, out.dir, mc.cores=mc.cores)
    #JZ:remove mc.cores arg
    }, phen, kinship.matrices, out.dir)
    dff=NULL
    for( r in df )
        if ( !is.null(r)) dff = rbind( dff, r)

#    print(dff)
    if ( ! is.null(dff)) {
        dff = data.frame(dff)
        print(dim(dff)) #JZ
        names(dff) = c( "phenotype", "ve", "vg", "h" )
        write.table(dff, file=out.file, quote=F, row=F, sep="\t")
    }
}


mixed.model.variance<- function( y, K ) {
 # print("starting mixed.model.variance")
   # estimates the variance components for a mixed model applied to the phenotype y and kinship matrix K.
   # Takes care of missing data in phenotype and ensures kinship matrix is subsetted appropriately
   # Returns a list which includes the square root of the inverse variance matrix, for subsequent genome scanning.

    use =!is.na(y)
    yy = y[use]
    nam = names(yy)
    int = intersect(nam, colnames(K))
    yy = y[match(int, names(yy), nomatch=0)]
    n = length(yy)
    idx = match(int, colnames(K), nomatch=0)
    K = K[idx,idx]

    #JZ: change df to matrix so emma.REMLE will work
    yy = as.matrix(yy)
    K = as.matrix(K)
#	print(summary(yy))
#	print(summary(K))
    mixed.model.data=emma.REMLE(y=yy, X=matrix(1,ncol=1,nrow=n), K=K, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.R = NULL)
    print(mixed.model.data)
    mixed.model.data$v = mixed.model.data$vg + mixed.model.data$ve
    V = mixed.model.data$vg*K+mixed.model.data$ve*diag(n)
    #JZ
#    print("mixed.model.variance: dim V...")
#    print(dim(V))
#    print("mixed.model.variance: dim K...")
#    print(dim(K))
    
    mixed.model.data$V=V
    ##creates multiplier matrix to apply to the mixed models to transform them into linear models
    svd=svd(V)
    
    eigen_basis=svd$u
    eigen_values=svd$d
    teigen_basis=t(eigen_basis)

    sev = 1/sqrt(eigen_values)
    sqrt_lambda_matrix=diag(sev)
    e = t(apply( eigen_basis, 1, function(x,sev) { x*sev }, sev ))
    mixed.model.data$multiplier = e %*% teigen_basis
    mixed.model.data$y.original = yy
    mixed.model.data$y = mixed.model.data$multiplier %*% yy
    names(mixed.model.data$y) = int
    mixed.model.data$ids = int

return(mixed.model.data)

}

mixed.model.allele.dosage.scans <- function( mixed.model.dir, out.dir, dosages=NULL ) {

    if ( is.null(dosages)) dosages = rwd.load.dosages()
    files = list.files( path=mixed.model.dir , pattern=".RData", full.name=TRUE )
    lapply( files, mixed.model.allele.dosage.scan, dosages, out.dir )
}
                             
mixed.model.allele.dosage.scan <- function( mixed.model.data.file, allele.dosages, out.dir=NULL, mc.cores=20 ) {

    if ( !is.null ( out.dir ) & ! file.exists(out.dir))
        dir.create( out.dir )
    else
        out.dir="./"

    load( mixed.model.data.file) # loads mixed.model.data into object called mm. This is produced by a previsous call of mixed.model.heritabilities()
    mm.all = mm$mm.all
    did = dosages$ids
    use = match( mm.all$ids, did, nomatch=0)
#    logP = mclapply( names(allele.dosages), function( chr, mm, allele.dosages, use ) {
    dfs = lapply( allele.dosages$dosages, function( dose, mm.all,  use ) {

        
        D = dose$dosages
        D = D[,use]
        D = D %*% mm.all$multiplier
        r = cor(  mm.all$y, D )
        r2 = r*r
        r2 = ifelse(is.na(r2), 0, r2)
        F = (r2*(n-2))/(1-r2+1.0e-20)
        logP = -pf(F, 1, n-2, lower.tail = FALSE, log.p = TRUE)/log(10)
        df = data.frame( allele.dosages$map[snps], logP=logP )
        return(df)
    }, mm.all, use ) #, mc.cores=20)
    df = do.call( "rbind", dfs )
    save(df, file=paste( out.dir, "/", mm$p, ".scan.RData")) 
}

mixed.model.allele.dosage.scan.genome <- function( mixed.model.data.file, chrs=1:19, out.dir=NULL, filter=TRUE, downsample=1, permute=FALSE, mc.cores=20 ) {

    if ( !is.null ( out.dir ) ) {
        if (! file.exists(out.dir))
            dir.create( out.dir )
    }
    else
        out.dir="./"

    # note call to lapply - parallism is done at lower lvel to aoing IO bottleneck
    
    lapply( chrs, mixed.model.allele.dosage.scan.per.chr, mixed.model.data.file, out.dir=out.dir, filter=filter, downsample=downsample, permute=permute, mc.cores=mc.cores )
}

mixed.model.allele.dosage.scan.per.chr <- function( chr, mixed.model.data.file, out.dir=NULL, filter=TRUE, downsample=1, permute=FALSE, mc.cores=20 ) {

    
    load( mixed.model.data.file) # loads mixed.model.data into object called mm. This is produced by a previous call to mixed.model.heritabilities()
    file=paste(out.dir, "/", mm$p, ".", chr, ".RData", sep="")
    if ( permute ) {
        mm$y = sample( mm$y, replace=FALSE )
        file=paste(out.dir, "/", mm$p, ".", chr, "permute.RData", sep="")
    }
    
    if ( ! file.exists(file)) {
        dosages = rwd.load.dosages.for.chr( type="allele", chr=chr, filter=filter )
        did = dosages$ids
        use = match( mm$ids, did, nomatch=0)
        D = dosages$dosages
        D = D[,use]
        n = nrow(D)
        map = dosages$map
        if ( downsample < 1 & downsample > 0 ) {
            size = round(downsample*n)
            snps = sort(sample(n, size=size, replace=FALSE ))
            D = D[snps,]
            map = map[snps,]
            cat("chr ", chr, n, "downsampled:", size, size/(n+1.0e-10), "\n")
        }
       
        step = floor(nrow(D)/mc.cores)+1
        segs = seq( 1, nrow(D), step)
        D = t(D)
        
        results = mclapply( segs, function( seg1, step, D, mm ) {
            seg2 = seg1+step
            if ( seg2 > ncol(D) ) seg2 = ncol(D)
            D = t(mm$multiplier) %*% D[,seg1:seg2]
            r = cor( D,  mm$y )
            nsub = length(mm$y)
            r2 = r*r
            r2 = ifelse(is.na(r2), 0, r2)
            F = (r2*(nsub-2))/(1-r2+1.0e-20)
            logP = -pf(F, 1, nsub-2, lower.tail = FALSE, log.p = TRUE)/log(10)
            return(logP)
        },  step-1, D, mm,  mc.cores=mc.cores )

        logP = do.call( "c", results )
        df = data.frame( bp=map$bp, logP=logP )
        data = list( phenotype=mm$p, chr=chr, df=df)
        
        cat("writing", file, "\n")
        save(data, file=file)
    }
}

    
rwd.scan.phenotypes <- function( outdir="/Net/mus/data/www/OUTBRED/RWD.SCANS.dosagesSubset/",
                                phenotype.file="/Net/sparse/data/scratch/nicod/Harwell_all_measures_resids_11May13.txt",
                                phenotype.dir = "/Net/sparse/data/scratch/nicod/new.residuals_11May13", 
                                normalise=TRUE,covariate.file=NULL, type="allele", permute=FALSE, mc.cores=20, dosages=NULL, mixed.model=FALSE) {

  if ( !is.null(phenotype.dir))
   p = read.phenotype.dir(phenotype.dir=phenotype.dir, normalise=normalise)
  else
    p = read.phenotypes(phenotype.file=phenotype.file, normalise=normalise)

    if ( is.null(dosages)) {
    d = rwd.load.dosages(type=type)
  }
  else {
    d = dosages
  }
  fast = ifelse ( d$type == "allele", TRUE, FALSE)

  if ( !is.null(covariate.file)) {
    covariates = read.table(covariate.file)
  }
  else
    covariates=NULL
  if ( ! file.exists(outdir)) dir.create(outdir) 
  for ( pk in keys(p) ) {
    cat(pk,"\n")
    f = paste( outdir, "/", pk, ".scan.RData", sep="")
    if ( ! file.exists(f)) {
      df = try( rwd.scan.dosages( d, p, pk, covariates=covariates, permute=permute, fast=fast, mc.cores=mc.cores ) )
      if ( !inherits( df, "try-error" ) ) {
        if ( !is.null(df)) save(df,file=f)
        cat("wrote ", f, "\n")
      }
      else {
        cat( "FAILED ", f, "\n")
      }
    }
  }
}

# Per chromosome functions for Robbies October data


