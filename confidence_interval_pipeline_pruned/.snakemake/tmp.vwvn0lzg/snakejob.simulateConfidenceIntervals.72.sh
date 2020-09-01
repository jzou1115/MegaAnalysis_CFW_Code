#!/bin/sh
# properties = {"type": "single", "rule": "simulateConfidenceIntervals", "local": false, "input": ["/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr1.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr2.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr3.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr4.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr5.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr6.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr7.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr8.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr9.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr10.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr11.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr12.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr13.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr14.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr15.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr16.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr17.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr18.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr19.oxinfo90.Rdata", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity30.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity5.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.begin.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.end.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.middle.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.total.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.basal.activity.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.a.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.n.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.decay.activity.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.EDL.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.baseline.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.corr.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.baseline.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.corr.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.corr.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.gastroc.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.glucose.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.diff.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.ratio.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.plantaris.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp12.ppi.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp6.ppi.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.soleus.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.startle.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tail.length.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.TA.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tibia.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.bmi.tibia.RData", "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.RData"], "output": ["/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals/bmd.a.11.95761488.CI.RData"], "wildcards": {"q": "bmd.a.11.95761488"}, "params": {"qtl": "bmd.a.11.95761488", "g_dir": "/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90", "mm_dir": "/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2", "ci_out": "/u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals"}, "log": [], "threads": 1, "resources": {}, "jobid": 72, "cluster": {"queue": "highp", "memory": "25G", "name": "simulateConfidenceIntervals", "output": "logs/simulateConfidenceIntervals.bmd.a.11.95761488.out", "error": "logs/simulateConfidenceIntervals.bmd.a.11.95761488.err", "mail": "ae", "time": "23:59:59", "cores": "1"}}
 cd /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline && \
/u/home/j/jzou1115/project-ernst/software/anaconda3/bin/python3.7 \
-m snakemake /u/home/j/jzou1115/project-zarlab/CFW/confidence_interval_pipeline_out/confidence_intervals/bmd.a.11.95761488.CI.RData --snakefile /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/.snakemake/tmp.vwvn0lzg /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/data/qtls/qtls_decorrelated_threshold2.txt /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr1.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr2.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr3.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr4.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr5.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr6.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr7.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr8.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr9.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr10.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr11.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr12.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr13.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr14.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr15.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr16.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr17.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr18.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/chr19.oxinfo90.Rdata /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity30.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity5.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.begin.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.end.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.middle.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.activity.total.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.basal.activity.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.a.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.n.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.bmd.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.decay.activity.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.EDL.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.baseline.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.corr.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.context.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.baseline.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.corr.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.cue.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.corr.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.fc.uncond.freeze.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.gastroc.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.glucose.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.diff.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.habit.ratio.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.plantaris.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp12.ppi.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.pp6.ppi.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.soleus.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.startle.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tail.length.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.TA.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.tibia.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.bmi.tibia.RData /u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.weight.RData --latency-wait 600 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules simulateConfidenceIntervals --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/.snakemake/tmp.vwvn0lzg/72.jobfinished || (touch /u/project/zarlab/jzou1115/CFW/MegaAnalysis_CFW_Code/confidence_interval_pipeline/.snakemake/tmp.vwvn0lzg/72.jobfailed; exit 1)

