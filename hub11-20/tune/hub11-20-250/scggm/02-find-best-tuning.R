#Find the best tuning selection from step 02-cv-xy-yy-scggm
#the tuning grid used
tsc02 <- expand.grid(lam1 = exp(seq(log(0.06), log(0.33), length  = 10)),
                     lam2 = exp(seq(log(0.04), log(0.33), length = 10)))
#the base directory
dir02 <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/02-cv-xy-yy-scggm"
ddirs <- list.dirs(path = dir02, full.names = TRUE, recursive = FALSE)
#omit the log directory
ddirs <- ddirs[!grepl(pattern = "logs", x = ddirs)]
#summarise best tunings
source("~/repos/sim-spacemap/hub11-20/tune/hub11-20-250/scggm/report-top-scggm-tunings.R")
res02 <- report_cv_score(ddirs)
saveRDS(res02, file  = file.path(dir02, "cv_score_scggm_summary.rds"))
