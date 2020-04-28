#!/bin/sh
#$ -N gemma
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=20G,highp,h_rt=23:59:59
#$ -m a
#$ -M jzou1115
# -t 4-4:1
#$ -t 4-37:1

pheno=$SGE_TASK_ID

for chrm in {1..19}; do 
outdir=/u/home/j/jzou1115/project-zarlab/CFW/ox1_ox2_association

#OX1 analysis
genotypes=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/ox1.chr${chrm}.oxinfo50.dosages.gz
phenotypes=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/pheno.onlygenoed.ox1.oxinfo50.meta.txt
anno=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/ox1.chr${chrm}.oxinfo50.dosages.annotation.txt
out=ox1_pheno_${pheno}_chr${chrm}
cov=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/covars.ox1.oxinfo50.meta.txt
/u/home/j/jzou1115/project-ernst/software/gemma-0.98.1-linux-static -g $genotypes -p $phenotypes -a $anno -lm 4 -o $out -n $pheno -c $cov


#OX2 analysis
genotypes=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/ox2.chr${chrm}.oxinfo50.dosages.gz
phenotypes=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/pheno.onlygenoed.ox2.oxinfo50.meta.txt
anno=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/ox2.chr${chrm}.oxinfo50.dosages.annotation.txt
out=ox2_pheno_${pheno}_chr${chrm}
cov=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo50/covars.ox2.oxinfo50.meta.txt
/u/home/j/jzou1115/project-ernst/software/gemma-0.98.1-linux-static -g $genotypes -p $phenotypes -a $anno -lm 4 -o $out -n $pheno -c $cov

done
