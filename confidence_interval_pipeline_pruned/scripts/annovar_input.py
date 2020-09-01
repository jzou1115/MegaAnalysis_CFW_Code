import pandas as pd
import os
import gzip
import sys

#        snps=expand(os.path.join(fine_dir, "output", "{q}_snps.txt"), q=qtlids),
#        r2=expand(os.path.join(os.path.join(outdir, "locuszoom", "{q}_r2.txt"), q=qtlids)

snpF = snakemake.input["snps"]
r2F = snakemake.input["r2"]
output = snakemake.output["outfile"]

#contains snps identified in susier and those in LD
snps = []
if os.stat(snpF).st_size != 0:
	d = pd.read_table(snpF, header=None)
	snps.extend(list(d[0]))

	d2 = pd.read_table(r2F)
	d2 = d2[d2["rsquare"]>.99]
	snps.extend(list(d2["snp1"]))

else:
	out = open(output, "w")
	out.close()
	sys.exit(0)
	
#snps = []
#for f in snp_files:
#	if os.stat(f).st_size != 0:
#		d = pd.read_table(f, header=None)
#		snps.extend(d[0])

snps = set(snps)

print(len(snps))
data = {}
data["chr"] = [t.split("_")[0] for t in snps]
data["pos"] = [t.split("_")[1] for t in snps]
data = pd.DataFrame(data)
data.index = snps
data = data.sort_values(by="chr")

gdir = "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90"
out = []
for c in set(data["chr"]):
	sub = data[data["chr"]==c]
	snps_c = set(sub.index)
    
	a1 = []
	a2 = []
	with gzip.open(os.path.join(gdir, "combined."+str(c)+".oxinfo90.dosages.gz"),'rt') as f:
		for line in f:
			tokens = line.split()
			if tokens[0] in snps_c:
				a1.append(tokens[1])
				a2.append(tokens[2])
                
	sub["allele1"] = a1
	sub["allele2"] = a2
    
	out.append(sub)


out = pd.concat(out, ignore_index=True)

out["end"] = out["pos"]

out.to_csv(output, header=False, index=False, sep="\t", columns = ["chr", "pos", "end", "allele1", "allele2"])
