#!/bin/sh
# properties = {"type": "single", "rule": "annovar_input", "local": false, "input": ["/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/finemapping/output/soleus.13.9242435_snps.txt", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/locuszoom/soleus.13.9242435_r2.txt"], "output": ["/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/finemapping/output/soleus.13.9242435_annovar_input.txt"], "wildcards": {"q": "soleus.13.9242435"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 231, "cluster": {"queue": "highp", "memory": "5G", "name": "annovar_input", "output": "logs/annovar_input.soleus.13.9242435.out", "error": "logs/annovar_input.soleus.13.9242435.err", "mail": "a", "time": "00:59:59", "cores": "1"}}
 cd /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned && \
/u/home/j/jzou1115/project-ernst/software/anaconda3/bin/python3.7 \
-m snakemake /u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/finemapping/output/soleus.13.9242435_annovar_input.txt --snakefile /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/.snakemake/tmp.q_fxjqco /u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/finemapping/output/soleus.13.9242435_snps.txt /u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out_pruned/locuszoom/soleus.13.9242435_r2.txt --latency-wait 600 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules annovar_input --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/.snakemake/tmp.q_fxjqco/231.jobfinished || (touch /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline_pruned/.snakemake/tmp.q_fxjqco/231.jobfailed; exit 1)

