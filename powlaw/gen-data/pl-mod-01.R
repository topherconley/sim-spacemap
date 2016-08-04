

#DESCRIPTION:
#Create the network: pl-mod-01

#IMPORTANT!!
#Use the laplacian or SPACE recipe for the covariance matrix creation. 
gen_under_laplacian <- TRUE
#reproducibility
set.seed(26225)

library(igraph)
source("~/repos/sim-spacemap/powlaw/gen-data/gen-power-law-module.R")

#generate 5 modules and stitch together node and attribute information
plgmods <- lapply(1:5, gen_powlaw_module, R = 100, S = 40, ntophubs = 4, prop_xhubs = 0.5, 
                  xbg_clusters = 5, prob_xbg_edge = 0.3)
gmods <- lapply(plgmods, function(x) x$mod)
cgpag <- do.call(what = union, args = gmods)

#label vertex types
vertex_groups <- c("xhubs", "yhubs", "xbg", "ybg") #, "xsubhubs", "xsubbg")
names(vertex_groups) <- vertex_groups
id_groups <- lapply(vertex_groups, function(x) Reduce('c',  sapply(plgmods, function(y) y[[x]])))
for (i in seq_along(id_groups)) { 
  cgpag <- set_vertex_attr(cgpag, name = "type", index = V(cgpag)[name %in% id_groups[[i]]], value =  names(id_groups)[i])
}
#vertex_attr_names(cgpag)

#generate weights among X--Y 
unif_x2y = c(0.5, 1)
xy_edges <- E(cgpag)[ V(cgpag)[type %in% "xhubs"] %--% V(cgpag)[type %in% c("yhubs", "ybg")]]
wgts <- runif(n =  length(xy_edges), min = unif_x2y[1], max =  unif_x2y[2])
#laplacian generation does not use negative weights currently
if (!gen_under_laplacian)  { 
  wgts <- ifelse(test = rbinom(n = length(xy_edges), size = 1, prob = 0.5), wgts, -1*wgts)
}
cgpag <- set_edge_attr(graph = cgpag, name = "weight", index = xy_edges, value = wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "type", index = xy_edges, value = "xy")

#generate weights among Y--Y edges 
unif_y2y = c(0.5, 1)
yy_edges <- E(cgpag)[ V(cgpag)[type %in% c("yhubs", "ybg")] %--% V(cgpag)[type %in% c("yhubs", "ybg")]]
wgts <- runif(n =  length(yy_edges), min = unif_y2y[1], max =  unif_y2y[2])
if (!gen_under_laplacian) { 
  wgts <- ifelse(test = rbinom(n = length(yy_edges), size = 1, prob = 0.5), wgts, -1*wgts)
}
cgpag <- set_edge_attr(graph = cgpag, name = "weight", index = yy_edges, value = wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "type", index = yy_edges, value = "yy")

#generate weights among X--X 
unif_x2x = c(0.2, 0.5)
xx_edges <- E(cgpag)[ V(cgpag)[type %in% c("xhubs", "xbg")] %--% V(cgpag)[type %in% c("xhubs", "xbg")]]
wgts <- runif(n =  length(xx_edges), min = unif_x2x[1], max =  unif_x2x[2])
cgpag <- set_edge_attr(graph = cgpag, name = "weight", index = xx_edges, value = wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "type", index = xx_edges, value = "xx")

#assure that all edges have been accounted for with a weight
stopifnot(ecount(cgpag) == sum(sapply(c(xx_edges, xy_edges, yy_edges), length)))
edge_attr_names(cgpag)

fit_power_law(x = degree(cgpag))


#set_edge_attr(graph = cgpag, name = "")
vertex_attr_names(cgpag)

#High modularity score
ceb <- cluster_edge_betweenness(graph = cgpag, directed = FALSE)
modularity(x = cgpag, membership = ceb$membership)

#cytoscape visualization
write.graph(graph = cgpag, file = "~/tmp/cg-module-power-law.gml", format = "gml")

hist(edge_attr(graph = cgpag, name = "weight"))

#Define the X and Y index sets
iy <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("yhubs", "ybg")])
ix <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("xhubs", "xbg")])
ixhubs <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("xhubs")])
ixbg <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("xbg")])
#assure non-overlapping
stopifnot(intersect(iy, ix) == integer(0))

#identify X confounders
xhubtoall <- induced_subgraph(graph = cgpag, vids = unlist(ego(graph = cgpag, order = 5, nodes = V(cgpag)[type %in% c("xhubs")])))
xhubtoxbg <- induced_subgraph(graph = xhubtoall, vids = V(xhubtoall)[type %in% c("xhubs", "xbg")])
xhub_to_xbg_edges <- E(xhubtoxbg)[ V(xhubtoxbg)[type %in% c("xhubs")] %--% V(xhubtoxbg)[type %in% c("xbg")]]
# xhub_to_xbg_edges <- E(cgpag)[ V(cgpag)[type %in% c("xhubs")] %--% V(cgpag)[type %in% c("xbg")]]
# xhub_to_xbg_edges
# ecount(xhubtoxbg)
# xhub_to_xbg_edges
# V(xhubtoxbg)$type
# plot(xhubtoxbg)
###############################################################
#Create the positive defininte Covariance/Precision matrices 
#and the Partial Correlation Matrix
tol <- 1e-6
if (gen_under_laplacian) { 
  
  #Define the laplacian := Precision Matrix
  # It is positive semi-definite and can be converted to 
  # positive definite (See sCGGM paper simulation).
  precision_xxyy <- laplacian_matrix(graph = cgpag, normalized = F, sparse = FALSE)
  #verify properties of laplacian
  library(matrixcalc)
  is.positive.definite(precision_xxyy)
  is.positive.semi.definite(precision_xxyy)
  is.symmetric.matrix(precision_xxyy)
  #Add small positive values to diagonal
  diag(precision_xxyy) <- diag(precision_xxyy) + 0.2
  stopifnot(is.positive.definite(precision_xxyy))
  
  #inspect precision diagonal and make sure that it is not all 1. 
  #if use normalized Laplacian then the diagonal is all 1's. 
  #We don't use the normalized Laplacian to make the 
  #estimation more challenging for Spacemap. 
  summary(diag(precision_xxyy))
  
  #Broken down by x and y blocks
  precision_xy <- precision_xxyy[ix,iy]
  precision_yy <- precision_xxyy[iy,iy]
  is.positive.definite(precision_yy)
  
  model_sig_xy <- (abs(precision_xy) > tol) + 0
  model_sig_yy <- (abs(precision_yy) > tol) + 0
  
  ###########################################
  #Validate the assumptions of the network topology with precision matrix
  Adj <- as_adj(graph = cgpag, type = "both", sparse = FALSE)
  diag(Adj) <- 1
  identical(abs(Adj) > tol, abs(precision_xxyy) > tol)
  
  library(spacemap)
  #Y to Y
  nonZeroUpper(as.matrix(precision_xxyy)[iy,iy], tol)
  #All X  to Y
  nonZeroWhole(as.matrix(precision_xxyy)[ix,iy], tol)
  #X hubs to Y
  nonZeroWhole(as.matrix(precision_xxyy)[ixhubs,iy], tol)
  #X->Y edges with multiple hits
  sum(degree(cgpag)[ixhubs])
  #no edges from X background to Y
  nonZeroWhole(as.matrix(precision_xxyy)[ixbg,iy], tol)
  
  #########################################
  #Convert Precision Matrix to Sigma Matrix
  sigma_xxyy <- chol2inv(chol(precision_xxyy))  
  colnames(sigma_xxyy) <- colnames(precision_xxyy)
  rownames(sigma_xxyy) <- rownames(precision_xxyy)
  stopifnot(is.positive.definite(sigma_xxyy))
  #check numerical stability: error is below machine tolerance
  kappa(sigma_xxyy)
  kappa(precision_xxyy)
  eye <- sigma_xxyy %*%precision_xxyy
  diag(eye) <- 0
  summary(abs(eye[upper.tri(eye)]))  
  #further evidence that solve function is numerically stable: reversible inversion
  #concentration matrix  < === inv(Sigma) < === inv(Precision)
  #Concentration matrix ought to equal precision
  conc_xxyy <- solve(sigma_xxyy)
  colnames(conc_xxyy) <- colnames(precision_xxyy)
  rownames(conc_xxyy) <- rownames(precision_xxyy)
  stopifnot(all.equal(conc_xxyy, precision_xxyy))

  library("corpcor")
  parcor_xxyy <- cor2pcor(m = sigma_xxyy)
  colnames(parcor_xxyy) <- colnames(precision_xxyy)
  rownames(parcor_xxyy) <- rownames(precision_xxyy)
  
  #partial correlation not numerically symmetric
  #but differences are well below tolerance considered in evaluation 
  #stopifnot(is.symmetric.matrix(parcor_xxyy))
  tmp <- parcor_xxyy - t(parcor_xxyy)
  summary(abs(tmp[upper.tri(tmp)]))
  #hist(log10(abs(tmp[upper.tri(tmp)])))

  #1 to 1 precision_xxyy and Partial Correlation
  stopifnot(identical(abs(parcor_xxyy) > tol, abs(conc_xxyy) > tol))
  stopifnot(identical(abs(parcor_xxyy) > tol, abs(precision_xxyy) > tol))
  
  
# Spacemap and corpccor package agree
#   parcor_xxyy <- conc2parcor(conc_xxyy)
#   colnames(parcor_xxyy) <- colnames(precision_xxyy)
#   rownames(parcor_xxyy) <- rownames(precision_xxyy)
#   diag(parcor_xxyy) <- 1
#   stopifnot(all.equal(parcor_xxyy, parcor_xxyy2))
  
  #Distribution on Partial Correlations
  pc <- parcor_xxyy
  diag(pc) <-  0
  nzpc <- (abs(pc) > tol)
  ut <- upper.tri(pc)
  pc[!nzpc] <- 0.0 
  summary(pc[nzpc & ut])
  hist(pc[nzpc & ut], 20)
  library(ggplot2)
  library(scales)
  qplot(x = pc[nzpc & ut]) + geom_histogram() + theme_bw() + 
    scale_x_continuous(breaks = pretty_breaks(15)) + 
    xlab("Partial Correlations")
  ggsave("~/Dropbox/Chris_Conley/ms-spacemap/figures/pl-mod-01-partial-correlations-hist.png")
  
  #degree distribution 
  xd <- degree(graph = cgpag, v= V(cgpag)[type %in% c("xhubs", "xbg")])
  yd <- degree(graph = cgpag, v= V(cgpag)[type %in% c("yhubs", "ybg")])
  dlab <- c(rep("X", length(xd)), rep("Y", length(yd)))
  ggplot(data= data.frame(degree = c(xd,yd), type = dlab), aes(x = degree)) + 
    facet_grid(type ~ . ) + geom_histogram() + theme_bw()
  ggsave("~/Dropbox/Chris_Conley/ms-spacemap/figures/pl-mod-01-degree-hist.png")
  
  # Alternate Method for generating multivariate normal: 
  #(Precision Parameterization -> Cholesky Decomposition -> Backsolve  )
  #
  # R <- chol(precision_xxyy)
  # Lqr <- qr(x = precision_xxyy)
  # 
  # rmvnorm_precision <- function(R) {
  #   z <- rnorm(n = nrow(R), mean = 0, sd = 1)
  #   v <- backsolve(r = R, x = z)
  # }
  # dat_precision <- t(replicate(n, rmvnorm_precision(R)))
  # dim(dat_precision)
  # 
  # summary(apply(dat_precision, 2, sd))
  # 
  # 
}



###############################################################

#generate simulation under covariance matrix
datadir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/datasets/n250/"
simdir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/"
system(paste("mkdir -p ", datadir))

#set up the partition
n <- 250
library(caret)
testSetIds <- createFolds(y = 1:n, k  = 10, returnTrain = FALSE)
trainSetIds <- lapply(testSetIds, function(x) setdiff(1:n, x))
saveRDS(list(testSetIds = testSetIds, trainSetIds = trainSetIds), file = file.path(datadir, "pl-mod-01-n250-train-test-partition.rds"))

#no overlap between train and test sets
stopifnot(any(sapply(1:10, function(i) length(intersect(trainSetIds[[i]], testSetIds[[i]])) == 0)))
                
library(mvtnorm)
for(i in 1:100) { 
  XXYY <- rmvnorm(n = n, mean  = rep(0.0, ncol(sigma_xxyy)), sigma = sigma_xxyy)
  XXYY <- scale(XXYY)
  colnames(XXYY) <- colnames(sigma_xxyy)
  out <- list(XY = XXYY, Y = XXYY[,iy], X = XXYY[,ix])
  saveRDS(object = out, file = file.path(datadir, paste0("pl-mod-01-dataid_", sprintf("%03d", i), ".rds")))
}

#Save the model input
cgmodel <- list(graph = cgpag, sigma_xxyy = sigma_xxyy, 
                precision_xxyy = precision_xxyy, parcor_xxyy = parcor_xxyy,
                model_sig_xy = model_sig_xy, model_sig_yy = model_sig_yy, tol = tol)
saveRDS(object = cgmodel, file  = file.path(simdir, "pl-mod-01-model.rds"))


###############################################################
#Generate covariance matrix under SPACE recipe

# } else { 
#   
#   #######Diagonal dominance(SPACE, Jasa, 2009, recipe)#########
#   
#   # tmp <- as.matrix(get.adjacency(graph = cgpag, type = "both", attr = "weight"))
#   # is.symmetric.matrix(tmp)
#   # is.positive.definite(tmp)
#   
#   dd <- 1
#   tmp <- as.matrix(get.adjacency(graph = cgpag, type = "upper", attr = "weight"))
#   rs <- rowSums(abs(tmp))
#   rescale <- ifelse(rs == 0, 1, 6.4*rs)
#   #tmp <- t(t(tmp)/ (2*rescale))
#   tmp <- t(sapply(1:nrow(tmp), function(i) tmp[i,] / rescale[i]))
#   tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
#   A <- tmp
#   A <- ( tmp + t(tmp) )  / 2
#   diag(A) <- dd
#   is.symmetric.matrix(A)
#   tmp <- A
#   tmp[lower.tri(tmp, diag = T)] <- 0.0
#   #diagonally dominant
#   offdiag <- rowSums(abs(tmp)) 
#   stopifnot(max(offdiag) < dd)
#   stopifnot(is.positive.definite(A))
#   #Ainv <- solve(A) 
#   #is.symmetric.matrix(Ainv) 
#   # FALSE
#   Ainv <- chol2inv(chol(A))
#   stopifnot(is.symmetric.matrix(Ainv))
#   kappa(A)
#   kappa(Ainv)
#   #numerical issue here
#   #NegSqrtDiagAinv <- diag(x = sqrt(1/diag(Ainv)))
#   #vsigma_xxyy <- NegSqrtDiagAinv %*% Ainv %*% NegSqrtDiagAinv
#   #is.positive.definite(vsigma_xxyy)
#   #FALSE 
#   
#   #diagonal is all ones
#   U <- 5*(R + S)
#   sigma_xxyy <- matrix(NA, U, U)
#   for(i in 1:U) for(j in 1:U) sigma_xxyy[i,j] <- (Ainv[i,j] / sqrt(Ainv[i,i]*Ainv[j,j]))
#   #all.equal(sigma_xxyy, vsigma_xxyy)
#   # [1] TRUE
#   stopifnot(all.equal(rep(1.0, 5*(R +S)), diag(sigma_xxyy)))
#   stopifnot(is.positive.definite(sigma_xxyy))
#   #very stable condition number
#   kappa(sigma_xxyy)
#   
#   
#   # cor_xxyy <- cov2cor(sigma_xxyy)
#   # kappa(cor_xxyy)
#   # is.positive.definite(cor_xxyy)
#   
#   #calculate partial correlations
#   library("corpcor")
#   parcor_xxyy <- cor2pcor(m = sigma_xxyy)
#   conc_xxyy <- solve(sigma_xxyy)
#   summary(diag(conc_xxyy))
# #   #Buggy Spacemap way
# #   parcor_xxyy <- conc2parcor(conc_xxyy)
# #   diag(parcor_xxyy) <- 1
# #   stopifnot(all.equal(parcor_xxyy, parcor_xxyy2))
#   pc <- parcor_xxyy
#   
#   stopifnot(identical(abs(parcor_xxyy) > tol, abs(conc_xxyy) > tol))
#   #stopifnot(identical(abs(parcor_xxyy) > tol, abs(conc_xxyy) > tol))
#   
# 
#   diag(pc) <-  0
#   nzpc <- abs(pc) > tol
#   summary(pc[nzpc])
#   ut <- upper.tri(pc)
#   hist(pc[nzpc & ut], 20)
#   summary(abs(pc[nzpc]))
#   
#   #Validate the assumptions of the network topology
# 
#   library(spacemap)
#   #Y to Y
#   nonZeroUpper(as.matrix(pc)[iy,iy], tol)
#   #All X  to Y
#   nonZeroWhole(as.matrix(pc)[ix,iy], tol)
#   #X hubs to Y
#   nonZeroWhole(as.matrix(pc)[ixhubs,iy], tol)
#   #27 X->Y edges with multiple hits (330 - 313)
#   sum(degree(cgpag)[ixhubs])
#   #no edges from X background to Y
#   nonZeroWhole(as.matrix(pc)[ixbg,iy], tol)
#   #trueParCor <- list(xy = precision_xy, yy = precision_yy)
# }


#Other graphs to consider

# #sample hierarchical SBM
# set.seed(15)
# #number of nodes
# n <- 1000
# #number of vertices per block
# m <- sample(x = 75:200, size = 7)
# rho <- runif(n = 7, min = 0.1, max = 0.2)


# ?sample_hierarchical_sbm
# ?sample_sbm
# ?sample_motifs
# ?sample_smallworld
# ?sample_pa
