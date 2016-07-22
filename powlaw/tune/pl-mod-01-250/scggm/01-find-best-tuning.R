#Find the best tuning selection from step 01-cv-xy-yy-scggm
#the base directory
dir_prev <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/scggm/01-cv-xy-yy-scggm/"
ddirs <- list.dirs(path = dir_prev, full.names = TRUE, recursive = FALSE)
#omit the log directory
ddirs <- ddirs[!grepl(pattern = "logs", x = ddirs)]
ddirs <- ddirs[c(1:24,26:29, 31:100)]
#summarise best tunings
source("/home/cconley/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/scggm/report-top-scggm-tunings.R")
#datasets 25,30 each  had one tune id/ one fold combination(1/360) that failed. This is okay to proceed. 
#we do not need to retune this
res_prev <- report_cv_score(ddirs = ddirs, accept_partial = FALSE)
saveRDS(res_prev, file  = file.path(dir_prev, "cv_score_scggm_summary.rds"))
