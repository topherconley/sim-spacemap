#Objective to assure that the cvVote procedure is not biased by a tolerance-ignorant cvVote procedure.
basedir <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/03-cv-xy-yy-spacemap/"
#foreach dataset
library(spacemap)
foreach(dataid = 1:100) %do% { 
  ddir <- file.path(basedir, paste0("d", sprintf("%03d", dataid)))
  #tfiles <- list.files(path = ddir, pattern = "tuneid_")
  #tid <- unique(as.integer(regmatches(tfiles, regexpr(pattern = "(?<=tuneid_)(.*?)(?=_fold)", 
  #                                          text = tfiles, perl = TRUE))))
  cvint <- readRDS(file = paste0(ddir, "/cv_object_dataid_", dataid, ".rds"))
  tid <- cvint$minIndices["rss"]                                 
  tf <- list.files(path = ddir, pattern = paste0("tuneid_", sprintf("%03d", tid)))
  cvrds <- cvVoteRDS(mod_lambda = tf, tol = 1e-6, major_thresh= 0.5, method = "spacemap")
  c(identical(cvrds$xy, cvint$cvFits$rss$vote$Gamma + 0),
    identical(cvrds$yy, cvint$cvFits$rss$vote$ParCor + 0))
}


#The prior cvVote procedure is biased. 
