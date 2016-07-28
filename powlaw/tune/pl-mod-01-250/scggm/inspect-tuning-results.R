
basedir <- "/home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/scggm/"
i01 <- readRDS(file.path(basedir, "01-cv-xy-yy-scggm", "cv_score_scggm_summary.rds"))
i02 <- readRDS(file.path(basedir, "02-cv-xy-yy-scggm", "cv_score_scggm_summary.rds"))

default_grid <- c(0.32, 0.16, 0.08, 0.04, 0.02, 0.01)
tsc01 <- expand.grid(lam1 = default_grid,
                     lam2 = default_grid)

#find the top tuning parameter for this specific dataset. 
did <- c(1:24, 26:29, 31:100)
library(foreach)
tsc02 <- foreach(dataid = did) %do% { 
  top_tsc01_id <- i01$top$best_id[i01$top$dataid == dataid]
  top_lam1 <- tsc01[top_tsc01_id,"lam1"]
  top_lam2 <- tsc01[top_tsc01_id,"lam2"]
  #create new grid
  tsc <- expand.grid(lam1 = seq(top_lam1*.60, top_lam1*1.90, length  = 6),
                     lam2 = seq(top_lam2*.30, top_lam2*1.90, length  = 6))
}
names(tsc02) <- did