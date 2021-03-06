import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import gzip
import os
import sys

args = sys.argv
i = int(args[1]) - 1


ci = pd.read_table("/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/confidence_intervals/CI_pruned.txt")
ci = ci.dropna(subset=["from.bp", "to.bp"])
ci_95 = ci[ci["quantile"]==.95]
indices = ci_95.index
i = indices[i]

g_dir="/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90"
outdir="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/ld"

pheno = ci_95.loc[i, "phenotype"]
chrm = "chr"+str(ci_95.loc[i, "chr"])  
bp = str(ci_95.loc[i, "bp"])
qtl_id = "_".join([pheno, chrm, bp])
print(qtl_id)
if os.path.isfile(outdir+"/"+qtl_id+"_r2.txt"):
    sys.exit(1)
    
start_ci = ci_95.loc[i, "from.bp"]
end_ci = ci_95.loc[i, "to.bp"]

locus = []
dosages = []

with gzip.open(g_dir+"/combined."+str(chrm)+".oxinfo90.dosages.gz",'rt') as f:
    for line in f:

        tokens = line.split("\t")
        snp_tokens = tokens[0].split("_")
        pos = int(snp_tokens[1])
        if pos >= start_ci and pos <= end_ci:
            locus.append(tokens[0])
        
            dos = tokens[3:]
            dos = [float(d) for d in dos]
            dosages.append(dos)
        
    
dosages = pd.DataFrame(dosages)
dosages.index = locus

ref_snp = chrm+"_"+str(bp)

out = open(outdir+"/"+qtl_id+"_r2.txt", "w")
out.write("snp1\tsnp2\trsquare\n")
for snp in dosages.index:
    if snp != ref_snp:
        temp = dosages.loc[[ref_snp, snp], :]
        temp = temp.transpose()
        r = temp.corr()
        r = r.loc[ref_snp,snp]
        r2 = r*r
            
        out.write(snp+"\t"+ref_snp+"\t"+str(r2)+"\n")
out.close()
