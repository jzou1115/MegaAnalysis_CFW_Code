import pandas as pd 
import numpy as np
import glob
import sys

args = sys.argv

n_pheno=int(args[1])

num_phenotypes=37

phenotype_table = pd.read_table("/u/home/j/jzou1115/project-zarlab/CFW/phenotypes/pheno.names_shyam.txt", header=None)



def formatGemmaOutput(n):
    phenotype = phenotype_table.loc[n-1, 0]
    print(phenotype, n)
    gwas_dir="/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association/gwas_out_cov"
    ox1_gwas_data = []
    ox2_gwas_data = []
    for chrm in range(1,20):
        f = pd.read_table(gwas_dir+"/ox1_pheno_"+str(n)+"_chr"+str(chrm)+".assoc.txt")
        ox1_gwas_data.append(f)
    
        f2 = pd.read_table(gwas_dir+"/ox2_pheno_"+str(n)+"_chr"+str(chrm)+".assoc.txt")
        ox2_gwas_data.append(f2)

    ox1_gwas_data = pd.concat(ox1_gwas_data, ignore_index=True)
    ox2_gwas_data = pd.concat(ox2_gwas_data, ignore_index=True)
    
    ox1_gwas_data["z"] = ox1_gwas_data["beta"]/ ox1_gwas_data["se"]
    ox2_gwas_data["z"] = ox2_gwas_data["beta"]/ ox2_gwas_data["se"]
    
    n_ox1=1037
    n_ox2=1036
    
    snps = list(set(ox1_gwas_data["rs"]).intersection(set(ox2_gwas_data["rs"])))
    ox1_gwas_data.index = ox1_gwas_data["rs"]
    ox2_gwas_data.index = ox2_gwas_data["rs"]
    
    ox1_gwas_data = ox1_gwas_data.loc[snps,:]
    ox2_gwas_data = ox2_gwas_data.loc[snps,:]
    
    out = {}
    out["s1"] = ox1_gwas_data["z"]
    out["s2"] = ox2_gwas_data["z"]
    out["n1"] = n_ox1
    out["n2"] = n_ox2
    out = pd.DataFrame(out)
    out.index = snps
    
    outdir = "/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association/model_in"
    out.to_csv(outdir+"/ox1_ox2_"+phenotype+".txt", sep="\t", index_label=False)


def formatGemmaOutput90(n):
    phenotype = phenotype_table.loc[n-1, 0]
    print(phenotype, n)
    gwas_dir="/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association/gwas_out_cov"
    ox1_gwas_data = []
    ox2_gwas_data = []
    for chrm in range(1,20):
        f = pd.read_table(gwas_dir+"/ox1_pheno_"+str(n)+"_chr"+str(chrm)+".assoc.txt")
        ox1_gwas_data.append(f)
    
        f2 = pd.read_table(gwas_dir+"/ox2_pheno_"+str(n)+"_chr"+str(chrm)+".assoc.txt")
        ox2_gwas_data.append(f2)

    ox1_gwas_data = pd.concat(ox1_gwas_data, ignore_index=True)
    ox2_gwas_data = pd.concat(ox2_gwas_data, ignore_index=True)
    
    ox1_gwas_data["z"] = ox1_gwas_data["beta"]/ ox1_gwas_data["se"]
    ox2_gwas_data["z"] = ox2_gwas_data["beta"]/ ox2_gwas_data["se"]
    
    n_ox1=1037
    n_ox2=1036
    
    snps = set(ox1_gwas_data["rs"]).intersection(set(ox2_gwas_data["rs"]))
    snps = list(snps.intersection(snps90))
    
    ox1_gwas_data.index = ox1_gwas_data["rs"]
    ox2_gwas_data.index = ox2_gwas_data["rs"]
    
    ox1_gwas_data = ox1_gwas_data.loc[snps,:]
    ox2_gwas_data = ox2_gwas_data.loc[snps,:]
    
    out = {}
    out["s1"] = ox1_gwas_data["z"]
    out["s2"] = ox2_gwas_data["z"]
    out["n1"] = n_ox1
    out["n2"] = n_ox2
    out = pd.DataFrame(out)
    out.index = snps
    
    outdir = "/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association/model_in"
    out.to_csv(outdir+"/ox1_ox2_"+phenotype+"_90.txt", sep="\t", index_label=False)

#get set of snps meeting info90 threshold
z_dir = "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/oxinfo90/assoc"
snps90 = set([])
for n in range(4, 38):
    pheno=phenotype_table.loc[n-1, 0]
    uc = pd.read_table(z_dir+"/"+pheno+".uc.assoc.txt.gz")
    ox = pd.read_table(z_dir+"/"+pheno+".ox.assoc.txt.gz")

    snps = set(uc["rs"]).union(set(ox["rs"]))

    snps90.update(snps)
    
print(len(snps90))
print("snps90")
formatGemmaOutput90(n_pheno)
print("unfiltered")
formatGemmaOutput(n_pheno)
