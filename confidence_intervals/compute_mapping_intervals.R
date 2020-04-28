#### Mapping ####

source("/Net/mus/data/www/OUTBRED/outbred.pruned.R")      
setwd("/Net/dense/data/nicod/NoRel_Analysis_R5.6")

# 0. Generate kinship matrix:
# Has already been done for PCAs so I just copy it here:
# system("cp /Net/dense/data/nicod/New_Analysis_R5.6/pruned.scaled.kinship.allele.no_rel.RData /Net/dense/data/nicod/NoRel_Analysis_R5.6/pruned.scaled.kinship.allele.no_rel.RData")

# 1. Load the kinship matrix:
load("/Net/dense/data/nicod/NoRel_Analysis_R5.6/pruned.scaled.kinship.allele.no_rel.RData")     
# 2. compute the mixed model object for each phenotypes:
# system("mkdir /Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/")
mixed.model.heritabilities(kinship.matrices, phenotypes=NULL, out.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/", normalise=TRUE, out.file="heritabilities_alleles.txt",phenotype.file="/home/nicod/final_12Mar13/Harwell_all_measures_resids_For_Article_Sub_15Oct15.txt",mc.cores = 20)

# To check that I have mixed models for all 200 phenotypes:
length(list.files("/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/"))
# [1] 200     # correct!

# 3. Genome scans
# system("mkdir /Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans")     # only once
mixed.model.pruned.dosage.scan.wrapper( mm.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/", out.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans", mc.cores=1, allele.dir="/Net/dense/data/scratch/nicod/Robbie_run5.6_puned_no_rel", permute=100, decreasing=FALSE, per.chromosome=TRUE, type="allele" )
# To speed things up I run it from sparse, in reverse order
mixed.model.pruned.dosage.scan.wrapper( mm.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/", out.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans", mc.cores=1, allele.dir="/Net/dense/data/scratch/nicod/Robbie_run5.6_puned_no_rel", permute=100, decreasing=TRUE, per.chromosome=TRUE, type="allele" )

# To check that I have results for all 20 chromosomes for all 200 phenotypes:
length(list.files("/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans"))
# [1] 4000   
length(list.files("/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans"))/200
# [1] 20          # Good! there are 20 files for each phenotype (20 chromosomes)

# 4. Identify QTLs     (3Mb distance threshold!)
qtl.scan( dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans/", chrs=c(1:19,"X"), low.thresh=3, dist.thresh=3.0e6, decreasing=FALSE, limit=NULL, pheno=NULL, qtl.file="qtls_norm_scaled_oct.RData", all.file="qtls_norm_scaled_oct.all.RData", mc.cores=40, fdr.thresh = 0.1)
# fdr.threshold is used to categorise QTLs above the lowest QTL with an FDR above that threshold

load("/Net/dense/data/nicod/NoRel_Analysis_R5.6/qtls_norm_scaled_oct.RData")    
dim(qtls[qtls$is.qtl==TRUE,c(1,2,7,8,16,18)])
# [1] 422  6   # Was 486 before...
 
# A lot(!) on ChrX:
table(qtls$chr[qtls$is.qtl==TRUE])
10 11 12 13 14 15 16 17 18 19  1  2  3  4  5  6  7  8  9  X
14 73 16 21  7 13  1 55  8  8 29 34  6 14 28 34  6 10 20 25
 
dim(qtls[qtls$fdr.adj<0.1,c(1,2,7,8,16)])
# [1] 407   5
dim(qtls[qtls$fdr.adj<0.05,c(1,2,7,8,16)])
# [1] 267   5     # FDR 5%
dim(qtls[qtls$fdr.adj<0.01,c(1,2,7,8,16)])
# [1] 139   5

#### Estimate CI with pruned set on QTLs with FDR<5% (n=267) #####
# 0. Generate a table in the same format as /Net/dense/data/nicod/QTLs_areas_1Dec14.txt
load("/Net/dense/data/nicod/NoRel_Analysis_R5.6/qtls_norm_scaled_oct.RData")    
table=qtls[qtls$fdr.adj<0.05,c(1,2,7,8,16)]
dim(table)
# [1] 267 5     # correct
table$threshold=table$logP-1
table$start=table$bp-1e6
table$end=table$bp+1e6
table$dist=table$end-table$start
write.table(table,"/Net/dense/data/nicod/QTLs_areas_5.6_no_rel_oct_3Mb_mapping_logP_drop_23Oct15.txt",sep="\t",quote=F,row.names=F)

# 1. Simulations:
source("/Net/mus/data/www/OUTBRED/outbred.pruned.R")      
setwd("/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015")

simulate.confidence.intervals( qtl.file="/Net/dense/data/nicod/QTLs_areas_5.6_no_rel_oct_3Mb_mapping_logP_drop_23Oct15.txt", n.sim=10000, window=5.0e6, out.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015", per.chromosome=TRUE, mixed.model=TRUE, mm.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/", dir="/Net/dense/data/scratch/nicod/Robbie_run5.6_puned_no_rel/", suffix=".prunedgen.final.maf001.0.98.RData" , mc.cores1=10, mc.cores2=10, reverse=FALSE)
# run also from sparse to speed-up:
simulate.confidence.intervals( qtl.file="/Net/dense/data/nicod/QTLs_areas_5.6_no_rel_oct_3Mb_mapping_logP_drop_23Oct15.txt", n.sim=10000, window=5.0e6, out.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015", per.chromosome=TRUE, mixed.model=TRUE, mm.dir="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct/", dir="/Net/dense/data/scratch/nicod/Robbie_run5.6_puned_no_rel/", suffix=".prunedgen.final.maf001.0.98.RData" , mc.cores1=10, mc.cores2=10, reverse=TRUE)

length(list.files("/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015", pattern="RData"))
# [1] 267

# 2. Calculate pruned CIs only by using an empty fine.dir:
# (The returned .txt file has all fine columns as NA)
evaluate.confidence.intervals.logP( qtl.dir = "/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015", scan.dir ="/Net/dense/data/nicod/NoRel_Analysis_R5.6/norm_scaled_oct_scans/", pdf.file="CI.pdf", quantiles=c(0.25,0.5,0.75,0.9,0.95), dist.thresh=5.0e6, CI.file ="/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015/CI_pruned_23Oct15.txt", fine.dir= "/Net/dense/data/scratch/nicod/Robbie_run5.5_qtl.rescan_norm_scaled249_peaks_pdf/", fine.suffix="_Robbie_Run5.5_Maxgen_13Feb15_peak_annot.txt" , mc.cores=10)
ci=read.delim("/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015/CI_pruned_23Oct15.txt")
dim(ci)/5
# [1] 267 4     

#### Identify individual QTLs ####

# Identify overlapping QTLs (i.e. if a QTLs lies in the 95% CI of a neighbour):
ci=read.delim("/Net/dense/data/nicod/NoRel_Analysis_R5.6/CI.23102015/CI_pruned_23Oct15.txt")
ci=ci[order(-ci$logP),]
ci=ci[ci$quantile==0.95,]
ci$in_CI_reverse=NA

# n.snps=NA bug has entirely disappeared:
dim(ci[is.na(ci$n.snps),])
# [1] 0 21

# I test all QTLs starting from the strongest (highest logP): if there is another QTL in the 95% CI of the target QTL with a higher logP I remove (or flag) the target QTL
for (i in 1:length(ci$logP)){
     if (dim(ci[ci$phenotype==ci$phenotype[i] & ci$chr==ci$chr[i] & ci$bp>=ci$from.bp[i] & ci$bp<=ci$to.bp[i] & ci$logP>ci$logP[i],])[1]>0) ci$in_CI_reverse[i]="in_CI95_reverse"
     }
write.table(ci,"/Net/dense/data/nicod/QTLs_with_ci_1Jul15.txt",sep="\t",quote=F,row.names=F)
length(ci$logP[!is.na(ci$in_CI_reverse)])
# [1] 6     # It means there are still 6 QTLs I can remove because they are in the 95% CI of higher neighbouring QTL (there were 216 in the previous version of qtl.rescan())

# They are:
ci[!is.na(ci$in_CI_reverse),c(1:4,7:8)]
                    phenotype     logP chr        bp  from.bp     to.bp
875       Micronucleus.Mn.NCE 6.349747  17  39648533 34673930  40852666
930            Muscles.Gast.g 5.415120  11  90609160 89619475  95607942
375        Cardio.ECG.PR_main 5.363608   2  91856591 86858619  91986511
430        EPM.Total.Distance 5.226970  14  85648701 80674861  86007942
610  FACS.CD3posCD8posCD44pos 5.157208  11 101330963 96337925 101865204
1225               Sleep.Ampl 5.064155  11  37717794 37041286  41772403
# Add FDR to the ci table (from qtls)
load("/Net/dense/data/nicod/NoRel_Analysis_R5.6/qtls_norm_scaled_oct.RData")    
qtls=qtls[qtls$fdr.adj<0.05,c(1,2,7,8,16)]
dim(qtls)
# [1] 267  5     # correct
ci$fdr.adj=NA
for(i in 1:dim(ci)[1]){
    ci$fdr.adj[i]=qtls$fdr.adj[as.character(qtls$phenotype)==as.character(ci$phenotype[i]) & as.character(qtls$chr)==as.character(ci$chr[i]) & qtls$bp==ci$bp[i]]
}
# After removing the QTLs in 95%CI and keeping only FDR<5% I am left with 261 QTLs:
ci=ci[is.na(ci$in_CI_reverse) & ci$fdr.adj<0.05,]

# It turns out that the QTLs on chr2 for ECG.PR_main is large and detected as 2 peaks in this table. I remove it:
ci[ci$phenotype=="Cardio.ECG.PR_main",]
ci=ci[-228,]

# I prepare a table with these 260 QTLs:
write.table(ci,"/Net/dense/data/nicod/QTLs_list_260_26Oct15.txt",sep="\t",quote=F,row.names=F)

# Table with the 260 individual QTLs with FDR<5%:
table=read.delim("/Net/dense/data/nicod/QTLs_list_260_26Oct15.txt")
table$logP=round(table$logP,2)
table$fdr.adj=round(table$fdr.adj,2)
table=table[order(table$phenotype,table$chr,table$bp),]
