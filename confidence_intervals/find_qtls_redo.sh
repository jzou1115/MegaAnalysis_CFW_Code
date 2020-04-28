dir=/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/confidence_intervals
while read qtl ; do 
	phenotype=`echo $qtl | awk '{print $1}'`
	chr=`echo $qtl | awk '{print $2}'`
	bp=`echo $qtl | awk '{print $3}'`
	file=${phenotype}.${chr}.${bp}.CI.RData 

	if [ ! -f ${dir}/${file} ] ; then
		echo $qtl >> /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_perm_redo.txt
	fi
done < /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_perm.txt
