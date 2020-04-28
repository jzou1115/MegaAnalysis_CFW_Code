import pandas as pd
import numpy as np
import sys

ind = int(sys.argv[1])
outdir = sys.argv[2]

mega = pd.read_table("/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/pval_corrected_perm/CombinedToAll.txt")
pval = mega.loc[ind, "combined.logp"]


nullF="/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019_Summary_Statistics/permutations/perm.p_wald.combined.txt"
null = pd.read_table(nullF, header=None)
null_dist = np.array(null[0])

def computePermPval(pval, null):
    pval = 10**(-1*pval)
    count = sum(pval>null)
    print(count*1.0/len(null))
    return count*1.0/len(null)

p_adj = computePermPval(pval, null_dist)


out = open(outdir+"/perm_p_"+str(ind)+".txt" , "w")
out.write(str(p_adj)+"\n")
out.close()

