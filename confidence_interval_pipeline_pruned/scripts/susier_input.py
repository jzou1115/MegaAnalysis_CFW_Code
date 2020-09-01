import sys
import pandas as pd
import numpy as np
import math
import os

#EDL.12.83885241
ci_f = snakemake.input[0]
qtl = snakemake.params[0]
outdir = snakemake.params[1]
quantile = snakemake.params[2]
out1 = snakemake.output[0]
out2 = snakemake.output[1]

ci= pd.read_table(ci_f)
ci = ci[ci["quantile"]==quantile]
#ci = pd.read_table("/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals/CI_decorrelatedThresholds_q95_genes.txt")
z_dir = "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019_Summary_Statistics/COMBINED"
g_dir = "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90"
#outdir = "/u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/input"

qtl_tokens = qtl.split(".")
phenotype = ".".join(qtl_tokens[0:len(qtl_tokens)-2])
chrm = int(qtl_tokens[-2])
bp = int(qtl_tokens[-1])

sub = ci[(ci["phenotype"]==phenotype) & (ci["chr"]==chrm) & (ci["bp"]==bp)]
i = sub.index[0]
start = ci.loc[i, "from.bp"]
end = ci.loc[i, "to.bp"]
print(start, end)

d_f = phenotype+".combined.assoc.txt.gz"
d = pd.read_table(os.path.join(z_dir, d_f))
d = d[d["chr"]==chrm]
print(d.shape)
#subset to snps in confidence interval region
d = d[(d["ps"]<= end) & (d["ps"]>= start)]
print(d.shape)
d.index = d["rs"]
d["z"] = d["beta"]/d["se"]
num_snps = d.shape[0]
d = d.dropna(subset=["z"])
print("Dropping %d snps with missing z" % (num_snps - d.shape[0]))
print(d.shape)

#prune snps
pruneF = "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals_pruned/snps/combined.chr"+str(chrm)+".prune.in" #includes pruned snps + qtls
prune = pd.read_table(pruneF, header=None)
snps = list(set(d["rs"]).intersection(set(prune[0])))

#subset to pruned snps
d = d.loc[snps, :]

#sort by position
d = d.sort_values(by="ps")
snps = d.index
z_out = d.loc[snps, ["rs", "z"]]
#z_out_f = "_".join([str(t) for t in [phenotype, str(chrm), str(bp)]])+"_pruned_sorted.z"
z_out.to_csv(os.path.join(outdir, out1), sep="\t", header=True, index=False)
print("output z")

snps = list(d.index)

#read genotypes and subset to region
g_f = "combined.chr"+str(chrm)+".oxinfo90.dosages.gz"
g = pd.read_table(os.path.join(g_dir, g_f), index_col = 0, header=None)
g = g.iloc[:, 2:] #skip allele columns
g = g.loc[snps,:]
g = g.transpose()
ld = g.corr(method="pearson")
#ld_out_f = "_".join([str(t) for t in [phenotype, str(chrm), str(bp)]])+"_pruned_sorted.ld"
ld.to_csv(os.path.join(outdir, out2), sep="\t", header=True, index=True, index_label=False)
    
