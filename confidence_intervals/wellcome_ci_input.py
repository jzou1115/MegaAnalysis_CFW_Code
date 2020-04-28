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
ci_95 = ci[ci["quantile"]==.95]
indices = ci_95.index
i = indices[i]

g_dir="/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90"
outdir="/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/wellcome_trust"

pheno = ci_95.loc[i, "phenotype"]
chrm = "chr"+str(ci_95.loc[i, "chr"])  
bp = str(ci_95.loc[i, "bp"])
qtl_id = "_".join([pheno, chrm, bp])
print(qtl_id)
if os.path.isfile(outdir+"/"+qtl_id+"_r2.txt"):
    sys.exit(1)
    
dist_threshold=5000000
start_ci = max(0, ci_95.loc[i, "bp"] - dist_threshold)
end_ci = ci_95.loc[i, "bp"] + dist_threshold

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
locuszoom = pd.read_table("/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019_Summary_Statistics/COMBINED/"+pheno+".combined.assoc_locuszoom.txt.gz")
locuszoom.index=locuszoom["MarkerName"]

out = open(outdir+"/"+qtl_id+"_input.txt", "w")
out.write("SNP\tR2\tP\n")
missing = []
for snp in dosages.index:
	if snp != ref_snp:
		temp = dosages.loc[[ref_snp, snp], :]
		temp = temp.transpose()
		r = temp.corr()
		r = r.loc[ref_snp,snp]
		r2 = r*r
		try:
			p=locuszoom.loc[snp, "P-value"]
			out.write(snp+"\t"+str(r2)+"\t"+str(p)+"\n")     
		except:
			missing.append(snp)
#        out.write(snp+"\t"+ref_snp+"\t"+str(r2)+"\n")
print(len(missing))
out.close()
