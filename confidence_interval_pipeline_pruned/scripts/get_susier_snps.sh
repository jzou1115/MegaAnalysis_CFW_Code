qtlfile=data/qtls_decorrelated_info90.txt

while read q ; do 

Rscript scripts/get_susier_snps.R /u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/output/${q}.RData /u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/input/${q}_pruned_sorted.z /u/home/j/jzou1115/project-zarlab/CFW/Finemapping_Analysis_decorrelated/output/${q}

done < $qtlfile
