#find best models : outputs the <cv_selected> object seen below. 
source("id-best-iteration.R")
best_holdout_files <- lapply(cv_selected$dataid, get_best_holdout_files, cv_selected = cv_selected)

##################################
#obtain true model
library(spacemap)
library(igraph)
mod <- readRDS("~/scratch-data/sim-spacemap/powlaw/pl-mod-01/pl-mod-01-model.rds")
truth <- list(xy = mod$model_sig_xy, yy = mod$model_sig_yy)

###################################
#identify X variable types
cgpag <- mod$graph
#identify X confounders
xhubtoall <- induced_subgraph(graph = cgpag, vids = unlist(ego(graph = cgpag, order = 5, nodes = V(cgpag)[type %in% c("xhubs")])))
xhubtoxbg <- induced_subgraph(graph = xhubtoall, vids = V(xhubtoall)[type %in% c("xhubs", "xbg")])
xhub_to_xbg_edges <- E(xhubtoxbg)[ V(xhubtoxbg)[type %in% c("xhubs")] %--% V(xhubtoxbg)[type %in% c("xbg")]]
xbg_to_xbg_edges <- E(xhubtoxbg)[ V(xhubtoxbg)[type %in% c("xbg")] %--% V(xhubtoxbg)[type %in% c("xbg")]]
#confounders (background X's connected to X hubs)
cnfdrs <- V(xhubtoxbg)$name[as.integer(V(xhubtoxbg)[type %in% "xbg"])]
# extra background X's not connected to X hubs
extrax <- setdiff(V(cgpag)$name[as.integer(V(cgpag)[type %in% "xbg"])], cnfdrs)
# X hubs
xhubs <- V(cgpag)$name[as.integer(V(cgpag)[type %in% "xhubs"])]
#direct edge confounders
direct_cnfdrs <- unlist(lapply(
  strsplit(x = attr(xhub_to_xbg_edges, which = "vnames"), 
           split = "|", fixed = TRUE), function(x) x[2]))
#assure these are direct confounders
stopifnot(length(intersect(direct_cnfdrs, V(xhubtoxbg)[type %in% c("xhubs")])) == 0)
#indirect confounders
indirect_cnfdrs <- setdiff(cnfdrs, direct_cnfdrs)

xtypes <- list(xhubs = xhubs, cnfdrs = cnfdrs, extrax = extrax)

xtype_format_fdr_perf <- function(fit, truth, xtype, tol) { 
  xtypeTruth <- truth$xy[xtype,]
  xtypeFit <- fit$xy[xtype,]
  reportPerf(fitParCor = xtypeFit, trueParCor =  xtypeTruth, YY = F, tol = 0)["fp"] 
}

library(foreach)
tol <- 0.0
voted_fdr_perf <- foreach(bhf = best_holdout_files, .packages = "spacemap") %do% {
  fit <- cvVoteRDS(mod_lambda = bhf, tol = tol, major_thresh = 0.5, method = "space")
  rownames(fit$xy) <- rownames(truth$xy)
  as.data.frame(t(sapply(xtypes, function(xtype) 
    xtype_format_fdr_perf(fit = fit, truth = truth, 
                          xtype = xtype, tol = tol))))
}

library(data.table)
voted <- rbindlist(voted_fdr_perf)
voted[,list(mean_xhubs = mean(xhubs.fp), mean_cnfdrs = mean(cnfdrs.fp), mean_extrax = mean(extrax.fp))]
voted[,list(sd_xhubs = sd(xhubs.fp), sd_cnfdrs = sd(cnfdrs.fp), sd_extrax = sd(extrax.fp))]
saveRDS(voted, file = file.path(basedir, "space_top_cv_voting_fdr_by_xtype.rds"))
save.image(file = file.path(basedir, "space_top_cv_voting_fdr_by_xtype.rda"))

