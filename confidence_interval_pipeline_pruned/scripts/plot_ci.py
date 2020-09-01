import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob 
from os.path import basename
import math
import os


z_dir = snakemake.input["z_dir"]
f = snakemake.input["ld"]
ci_file = snakemake.input["ci"]
outFile = snakemake.output["outfile"]


ci = pd.read_table(ci_file)
ci = ci[ci["quantile"]==.95]

def plot_locus(ldFile, zFile, title, outFile):
    ld = pd.read_table(ldFile)
    z = pd.read_table(zFile, index_col=1)
    snps = set(ld["snp1"]).intersection(set(z.index))
    ld.index = ld["snp1"]
    ld = ld.loc[snps, :]
    z = z.loc[snps, :]
    z["logp"] = -np.log10(z["p_wald"])
    z["rsquare"] = ld["rsquare"]
    z["mb"] = z["ps"]/1000000
    
    plt.scatter("mb", "logp", c="rsquare", data=z)
    plt.colorbar()
    plt.xlabel("Position (MB)")
    plt.ylabel("-log10(p-value)")
    plt.title(pheno+"_"+chrm+"_"+str(bp))
    plt.savefig(outFile)
    plt.clf()
    
def plot_locus_ci(ldFile, zFile, title, outFile, start, end):
    ld = pd.read_table(ldFile)
    z = pd.read_table(zFile, index_col=1)
    snps = set(ld["snp1"]).intersection(set(z.index))
    ld.index = ld["snp1"]
    ld = ld.loc[snps, :]
    z = z.loc[snps, :]
    z["logp"] = -np.log10(z["p_wald"])
    z["rsquare"] = ld["rsquare"]
    z["mb"] = z["ps"]/1000000
    
    plt.scatter("mb", "logp", c="rsquare", data=z)
    plt.colorbar()
    plt.xlabel("Position (MB)")
    plt.ylabel("-log10(p-value)")
    plt.title(pheno+"_"+chrm+"_"+str(bp))
    plt.axvline(x=start/1000000)
    plt.axvline(x=end/1000000)
    plt.savefig(outFile)
    plt.clf()


baseF = f.split("/")[-1].split("_")[0]
tokens = baseF.split(".")
pheno = ".".join(tokens[0:len(tokens)-2])
chrm = tokens[-2]
bp = tokens[-1]
print(pheno, chrm, bp)

zFile = z_dir+"/"+pheno+".combined.assoc.txt.gz"
title = pheno+"_"+chrm+"_"+str(bp)


ci_qtl = ci[(ci["phenotype"]==pheno) & (ci["bp"]==int(bp)) ]

    
i = ci_qtl.index[0]
start_ci = ci.loc[i, "from.bp"]
end_ci = ci.loc[i, "to.bp"]
if math.isnan(ci.loc[i, "from.bp"]):
    plot_locus(f, zFile, title, outFile)
else:
    plot_locus_ci(f, zFile, title, outFile, start_ci, end_ci)
