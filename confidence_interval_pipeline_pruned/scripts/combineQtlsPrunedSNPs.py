import pandas as pd

prunedF = snakemake.input[0]
qtlF = snakemake.input[1]
outF = snakemake.output[0]
chrm = snakemake.params[0]

qtls = pd.read_table(qtlF, delim_whitespace=True)
qtlids = []
for i in qtls.index:
	if str(qtls.loc[i,"chr"]) == chrm:
		qtlids.append("chr"+str(qtls.loc[i,"chr"])+"_"+str(qtls.loc[i,"bp"]))

pruned = pd.read_table(prunedF, header=None)
out = open(outF, "w")
snpids = set([])
for i in pruned[0]:
	snpids.add(i)
	out.write(i+"\n")

for i in qtlids:
	if i not in snpids:
		out.write(i+"\n")
out.close()

