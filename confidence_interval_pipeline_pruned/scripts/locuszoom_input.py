import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import gzip
import os
import sys

ciF = snakemake.input["ci"]
g_dir = snakemake.input["d_dir"]
qtl_id = snakemake.params["qtl"]
outdir = snakemake.params["outdir"]
quantile = snakemake.params["quantile"]


qtl_tokens = qtl_id.split(".")
pheno = ".".join(qtl_tokens[:-2])
chrm = int(qtl_tokens[-2])
bp = int(qtl_tokens[-1])

ci = pd.read_table(ciF)
ci = ci.dropna(subset=["from.bp", "to.bp"])
ci_95 = ci[ci["quantile"]==quantile]
i = None

for ind in ci_95.index:
	if ci_95.loc[ind, "phenotype"]==pheno and ci_95.loc[ind, "chr"] == chrm and ci_95.loc[ind, "bp"] == bp:
		i = ind
		break		
    
#start_ci = ci_95.loc[i, "from.bp"]
#end_ci = ci_95.loc[i, "to.bp"]

#use fixed window instead of just confidence interval window
window=6000000
start_ci = max(0, bp - window)
end_ci = bp + window

print(start_ci, end_ci)
locus = []
dosages = []

with gzip.open(os.path.join(g_dir,"combined.chr"+str(chrm)+".oxinfo90.dosages.gz"),'rt') as f:
    for line in f:

        tokens = line.split("\t")
        snp_tokens = tokens[0].split("_")
        pos = int(snp_tokens[1])
        if pos >= start_ci and pos <= end_ci:
            locus.append(tokens[0])
        
            dos = tokens[3:]
            dos = [float(d) for d in dos]
            dosages.append(dos)
        
print(len(locus))  
dosages = pd.DataFrame(dosages)
dosages.index = locus

ref_snp = "chr"+str(chrm)+"_"+str(bp)

out = open(os.path.join(outdir, qtl_id+"_r2.txt"), "w")
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
