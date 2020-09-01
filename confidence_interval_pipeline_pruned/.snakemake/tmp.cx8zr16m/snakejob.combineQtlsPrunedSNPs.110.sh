#!/bin/sh
# properties = {"type": "single", "rule": "combineQtlsPrunedSNPs", "local": false, "input": ["/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/oxinfo90/LDprunedSNPs/combined.chr10.prune.in", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt"], "output": ["/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/combined.chr10.oxinfo90.txt"], "wildcards": {"chr": "10"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 110, "cluster": {"queue": "highp", "memory": "30G", "name": "combineQtlsPrunedSNPs", "output": "logs/combineQtlsPrunedSNPs.chr=10.out", "error": "logs/combineQtlsPrunedSNPs.chr=10.err", "mail": "abe", "time": "23:59:59", "cores": "1"}}
 cd /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned && \
/u/home/j/jzou1115/project-ernst/software/anaconda3/bin/python3.7 \
-m snakemake /u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/genotypes/combined.chr10.oxinfo90.txt --snakefile /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/.snakemake/tmp.cx8zr16m /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/oxinfo90/LDprunedSNPs/combined.chr10.prune.in /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt --latency-wait 600 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules combineQtlsPrunedSNPs --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/.snakemake/tmp.cx8zr16m/110.jobfinished || (touch /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/.snakemake/tmp.cx8zr16m/110.jobfailed; exit 1)

