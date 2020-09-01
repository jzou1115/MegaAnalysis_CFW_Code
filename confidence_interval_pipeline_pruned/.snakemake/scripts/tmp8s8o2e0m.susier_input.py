
######## Snakemake header ########
import sys; sys.path.extend(['/u/home/j/jzou1115/project-ernst/software/anaconda3/lib/python3.7/site-packages', '/u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05Xi\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals/CI_pruned.txtq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0c(Xm\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/input/EDL.12.83885241_pruned_sorted.zq\rXn\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/input/EDL.12.83885241_pruned_sorted.ldq\x0ee}q\x0fh\x08}q\x10sbX\x06\x00\x00\x00paramsq\x11csnakemake.io\nParams\nq\x12)\x81q\x13(X\x0f\x00\x00\x00EDL.12.83885241q\x14XM\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/inputq\x15e}q\x16(h\x08}q\x17(X\x03\x00\x00\x00qtlq\x18K\x00N\x86q\x19X\x06\x00\x00\x00outdirq\x1aK\x01N\x86q\x1buh\x18h\x14h\x1ah\x15ubX\t\x00\x00\x00wildcardsq\x1ccsnakemake.io\nWildcards\nq\x1d)\x81q\x1eh\x14a}q\x1f(h\x08}q X\x01\x00\x00\x00qq!K\x00N\x86q"sh!h\x14ubX\x07\x00\x00\x00threadsq#K\x01X\t\x00\x00\x00resourcesq$csnakemake.io\nResources\nq%)\x81q&(K\x01K\x01e}q\'(h\x08}q((X\x06\x00\x00\x00_coresq)K\x00N\x86q*X\x06\x00\x00\x00_nodesq+K\x01N\x86q,uh)K\x01h+K\x01ubX\x03\x00\x00\x00logq-csnakemake.io\nLog\nq.)\x81q/}q0h\x08}q1sbX\x06\x00\x00\x00configq2}q3(X\x06\x00\x00\x00outdirq4XF\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_outq5X\x05\x00\x00\x00mmdirq6XB\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2q7X\x04\x00\x00\x00gdirq8XN\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90q9X\x04\x00\x00\x00kdirq:XK\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/grmsq;X\x05\x00\x00\x00cidirq<X[\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervalsq=X\x07\x00\x00\x00scandirq>XH\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/scansq?X\x08\x00\x00\x00fine_dirq@XG\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelatedqAX\x0b\x00\x00\x00chromosomesqBX\x14\x00\x00\x00data/chromosomes.txtqCX\n\x00\x00\x00phenotypesqDX\x13\x00\x00\x00data/phenotypes.txtqEX\x04\x00\x00\x00qtlsqFXe\x00\x00\x00/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txtqGX\x07\x00\x00\x00qtl_idsqHX!\x00\x00\x00data/qtls_decorrelated_info90.txtqIuX\x04\x00\x00\x00ruleqJX\x0c\x00\x00\x00susier_inputqKX\x0f\x00\x00\x00bench_iterationqLNX\t\x00\x00\x00scriptdirqMXY\x00\x00\x00/u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/scriptsqNub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/scripts/susier_input.py';
######## Original script #########
import sys
import pandas as pd
import numpy as np
import math
import os

#EDL.12.83885241
ci_f = snakemake.input[0]
qtl = snakemake.params[0]
outdir = snakemake.params[1]

ci= pd.read_table(ci_f)
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
d = d[(d["ps"]< end) & (d["ps"]> start)]
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
z_out_f = "_".join([str(t) for t in [phenotype, str(chrm), str(bp)]])+"_pruned_sorted.z"
z_out.to_csv(os.path.join(outdir, z_out_f), sep="\t", header=True, index=False)
print("output z")

snps = list(d.index)

#read genotypes and subset to region
g_f = "combined.chr"+str(chrm)+".oxinfo90.dosages.gz"
g = pd.read_table(os.path.join(g_dir, g_f), index_col = 0, header=None)
g = g.iloc[:, 2:] #skip allele columns
g = g.loc[snps,:]
g = g.transpose()
ld = g.corr(method="pearson")
ld_out_f = "_".join([str(t) for t in [phenotype, str(chrm), str(bp)]])+"_pruned_sorted.ld"
ld.to_csv(os.path.join(outdir, ld_out_f), sep="\t", header=True, index=True, index_label=False)
    
