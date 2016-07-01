#Find the best tuning selection from step 03-cv-xy-yy-scggm
#the base directory
dir03 <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/03-cv-xy-yy-scggm"
ddirs <- list.dirs(path = dir03, full.names = TRUE, recursive = FALSE)
#omit the log directory
ddirs <- ddirs[!grepl(pattern = "logs", x = ddirs)]
#summarise best tunings
source("~/repos/sim-spacemap/hub11-20/tune/hub11-20-250/scggm/report-top-scggm-tunings.R")
res03 <- report_cv_score(ddirs)
saveRDS(res03, file  = file.path(dir03, "cv_score_scggm_summary.rds"))
