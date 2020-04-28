import pandas as pd

geno_dir = "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50"

for chrm in range(1, 20):
    
    genotype_f = geno_dir+"/"+"ox1.chr"+str(chrm)+".oxinfo50.dosages.gz"
    f = pd.read_table(genotype_f, header=None)
    f.head()
    
    out = open(geno_dir+"/"+"ox1.chr"+str(chrm)+".oxinfo50.dosages.annotation.txt", "w")
    for i in f.index:
        snp = f.loc[i, 0]
        snp_tokens = snp.split("_")
        c = snp_tokens[0]
        pos = snp_tokens[1]
        
        out.write(snp+"\t"+pos+"\t"+c+"\n")
    out.close()

for chrm in range(1, 20):
    
    genotype_f = geno_dir+"/"+"ox2.chr"+str(chrm)+".oxinfo50.dosages.gz"
    f = pd.read_table(genotype_f, header=None)
    f.head()
    
    out = open(geno_dir+"/"+"ox2.chr"+str(chrm)+".oxinfo50.dosages.annotation.txt", "w")
    for i in f.index:
        snp = f.loc[i, 0]
        snp_tokens = snp.split("_")
        c = snp_tokens[0]
        pos = snp_tokens[1]
        
        out.write(snp+"\t"+pos+"\t"+c+"\n")
    out.close()

