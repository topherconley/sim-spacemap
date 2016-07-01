#Find the best tuning selection from step 01-cv-xy-yy-scggm
#the base directory
dir01 <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/01-cv-xy-yy-scggm"
ddirs <- list.dirs(path = dir01, full.names = TRUE, recursive = FALSE)
#omit the log directory
ddirs <- ddirs[!grepl(pattern = "logs", x = ddirs)]
#summarise best tunings
source("~/repos/sim-spacemap/hub11-20/tune/hub11-20-250/scggm/report-top-scggm-tunings.R")
res01 <- report_cv_score(ddirs)
saveRDS(res01, file  = file.path(dir01, "cv_score_scggm_summary.rds"))
