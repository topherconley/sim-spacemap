
library(igraph)

gen_powlaw_module <- function(id, R, S, ntophubs, prop_xhubs, 
                              xbg_clusters = 5, prob_xbg_edge = 0.3) { 
  library(igraph)
  #out-degree of new vertices being attached
  out.seq <- pmax(1,rpois(n = R, lambda = 2.2)) #2.2
  #simulate a Barbasi-Albert scale-free preferential attachment network with R nodes
  pag <- sample_pa(n = R, power = 1.3, out.seq =  out.seq, out.pref = TRUE, directed = FALSE)
  #the number of edges
  #ecount(pag)
  #the number of vertices
  #vcount(pag)
  deg_pag <- degree(pag)
  #degree distribution summary
  #summary(deg_pag)
  #boxplot(degree(pag))
  #Inspect the power law fit
  plf <- power.law.fit(x = deg_pag)
  
  #Identify hubs for X and Y
  hub_pag <- hub.score(pag, scale = TRUE, weights=NULL, options = igraph.arpack.default)
  #plot(hub_pag$vector, deg_pag, pch = 19)
  #top hubs split in half among X and Y
  tophubs_index <- sort(hub_pag$vector, index.return = TRUE, decreasing = TRUE)$ix[1:ntophubs]
  xhubs_index <- sample(tophubs_index, size = ceiling(length(tophubs_index)*prop_xhubs))
  yhubs_index <- setdiff(tophubs_index, xhubs_index)
  #hist(degree(pag)[xhubs_index])
  #hist(degree(pag)[yhubs_index])
  
  #Identify the non-hub Y nodes
  ybg_index <- setdiff(1:R, c(xhubs_index, yhubs_index))
  
  #How many x-x hub edges exist? 
  #very few edges  ut of the whole graph  and the total X hub edges, okay to efffectively delete X--X hub edges
  xx_edges <- E(pag)[ xhubs_index %--% xhubs_index]
  #delete these edges from the graph as we are trying to create no geographic (genomic) correlation between true hubs
  cgpag <- delete_edges(pag, xx_edges)
  #ecount(cgpag)
  #new degrees for X hubs
  #degree(cgpag, v = xhubs_index)

  #Still remains a power-law with 2 < alpha < 3 in the range expected for biological networks
  noxx_plf <- power.law.fit(x = degree(cgpag))
  #noxx_plf$alpha
  
  #Assert the module has no disconnected components among the original power-law distribution
  stopifnot(count_components(cgpag) == 1)
  
  #add additional X vertices of size @param<S> which will have some edges between them. 
  cgpag <- add_vertices(graph = cgpag, nv = S)
  xbg_index <- (R + 1):(R + S)
  
  #Create X--X subgraphs of size @param<xbg_clusters> from an Erdos-Reny model, which has the potential to include an Xhub. 
  start <- R + 1
  ostart <- start
  for (xhub in xhubs_index) { 
    xid <- c(xhub, start:(start + xbg_clusters - 1))
    tmpg <- sample_gnp(n = length(xid), p  = prob_xbg_edge)
    tmpg <- set_vertex_attr(graph = tmpg, name = "name", value = xid)
    cgpag <- add_edges(cgpag, edges = as.vector(t(get.edgelist(tmpg))))
    start <- start + 5
  }
  #set of X's in X clusters with potential connection to xhubs
  xset1 <- ostart:(start - 1)
  
  #Create X--X subgraphs of size @param<xbg_clusters> from an Erdos-Reny model but do not include an X hub. 
  for(xbs in seq(from = start, to = (R + S), by = xbg_clusters)) {
    xid <- xbs:(xbs + xbg_clusters - 1)
    tmpg <- sample_gnp(n = length(xid), p  = prob_xbg_edge)
    tmpg <- set_vertex_attr(graph = tmpg, name = "name", value = xid)
    cgpag <- add_edges(cgpag, edges = as.vector(t(get.edgelist(tmpg))))
  }
  #set of X's in X clusters with no connection to xhubs
  xset2 <- start:(R + S)
  
  #unique vertex id
  #set_vertex_attr(cgpag, name = "name", index = V(cgpag), value = paste(id, V(cgpag), sep = ":"))
  V(cgpag)$name <- paste(id, V(cgpag), sep = ":")
  #cgpag <- set_vertex_attr(cgpag, name = "module_id", index = V(cgpag), value = id)
  list(mod = cgpag, 
       xhubs = paste(id, xhubs_index, sep = ":"), 
       yhubs = paste(id, yhubs_index, sep = ":"),
       xbg  = paste(id, xbg_index, sep = ":"),
       ybg = paste(id, ybg_index, sep = ":"), 
       xsubhubs = paste(id, xset1, sep = ":"), 
       xsubbg = paste(id, xset2, sep = ":"))
}

#number of nodes with the potential for an edge
R <- 100
#number of X nodes with no edge
S <- 40
#reproducibility
set.seed(26225)

#generate 5 modules and stitch together information
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

V(cgpag)[type %in% "xsubbg"]
V(cgpag)[type %in% vertex_groups[1]]
V(cgpag)[name %in% id_groups[["xhubs"]]]
#

#generate weights among X--Y 
unif_x2y = c(0.5, 0.8)
xy_edges <- E(cgpag)[ V(cgpag)[type %in% "xhubs"] %--% V(cgpag)[type %in% c("yhubs", "ybg")]]
wgts <- runif(n =  length(xy_edges), min = unif_x2y[1], max =  unif_x2y[2])
#wgts <- ifelse(test = rbinom(n = ecount(cgpag), size = 1, prob = 0.5), wgts, -1*wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "weight", index = xy_edges, value = wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "type", index = xy_edges, value = "xy")

#generate weights among Y--Y edges 
unif_y2y = c(0.5, 0.8)
yy_edges <- E(cgpag)[ V(cgpag)[type %in% c("yhubs", "ybg")] %--% V(cgpag)[type %in% c("yhubs", "ybg")]]
wgts <- runif(n =  length(yy_edges), min = unif_y2y[1], max =  unif_y2y[2])
#wgts <- ifelse(test = rbinom(n = ecount(cgpag), size = 1, prob = 0.5), wgts, -1*wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "weight", index = yy_edges, value = wgts)
cgpag <- set_edge_attr(graph = cgpag, name = "type", index = yy_edges, value = "yy")

#generate weights among X--X 
unif_x2x = c(0.5, 0.8)
xx_edges <- E(cgpag)[ V(cgpag)[type %in% c("xhubs", "xbg")] %--% V(cgpag)[type %in% c("xhubs", "xbg")]]
wgts <- runif(n =  length(xx_edges), min = unif_x2x[1], max =  unif_x2x[2])
#wgts <- ifelse(test = rbinom(n = ecount(cgpag), size = 1, prob = 0.5), wgts, -1*wgts)
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


write.graph(graph = cgpag, file = "~/tmp/cg-module-power-law.gml", format = "gml")
write.graph(graph = induced.subgraph(cgpag, vids = 1:R), file = "~/tmp/no-xx-cg-power-law.gml", format = "gml")

hist(edge_attr(graph = cgpag, name = "weight"))



################################################################################



#Define the laplacian := Precision Matrix
# It is positive semi-definite and can be converted to 
# positive definite (See sCGGM paper simulation).
#D <- diag(x = degree(cgpag))
#precision_xxyy <- D - A
precision_xxyy <- laplacian_matrix(graph = cgpag, normalized = TRUE, sparse = FALSE)
#verify properties of laplacian
library(matrixcalc)
is.positive.definite(precision_xxyy)
is.positive.semi.definite(precision_xxyy)
is.symmetric.matrix(precision_xxyy)
#Add small positive values to diagonal
diag(precision_xxyy) <- diag(precision_xxyy) + 0.05
stopifnot(is.positive.definite(precision_xxyy))
#Define the X and Y index sets


iy <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("yhubs", "ybg")])
ix <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("xhubs", "xbg")])
ixhubs <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("xhubs")])
ixbg <- vertex_attr(cgpag, name = "name", index = V(cgpag)[type %in% c("xbg")])
#assure non-overlapping
stopifnot(intersect(iy, ix) == integer(0))
#Broken down by x and y blocks
L_xy <- precision_xxyy[ix,iy]
L_yy <- precision_xxyy[iy,iy]
is.positive.definite(L_yy)

#Validate the assumptions of the network topology
tol <- 1e-6
library(spacemap)
#Y to Y
nonZeroUpper(as.matrix(precision_xxyy)[iy,iy], tol)
#All X  to Y
nonZeroWhole(as.matrix(precision_xxyy)[ix,iy], tol)
#X hubs to Y
nonZeroWhole(as.matrix(precision_xxyy)[ixhubs,iy], tol)
#X->Y edges with multiple hits
sum(degree(cgpag)[xhubs_index])
#no edges from X background to Y
nonZeroWhole(as.matrix(precision_xxyy)[ixbg,iy], tol)


#Convert precision_xxyy Matrix to Sigma Matrix
sigma_xxyy <- chol2inv(chol(precision_xxyy))  
stopifnot(is.positive.definite(sigma_xxyy))
summary(diag(sigma_xxyy))
#generate simulation
library(mvtnorm)
n <- 250
dat <- rmvnorm(n = n, mean  = rep(0.0, ncol(sigma_xxyy)), sigma = sigma_xxyy)

summary(apply(dat, 2, sd))

kappa(sigma_xxyy)
kappa(precision_xxyy)


library("corpcor")
parcor_xxyy2 <- cor2pcor(m = sigma_xxyy)
#Spacemap way
conc_xxyy <- solve(sigma_xxyy)
colnames(parcor_xxyy) <- colnames(precision_xxyy)
rownames(parcor_xxyy) <- rownames(precision_xxyy)
kappa(conc_xxyy)
summary(diag(conc_xxyy))
parcor_xxyy <- conc2parcor(conc_xxyy)
colnames(parcor_xxyy) <- colnames(precision_xxyy)
rownames(parcor_xxyy) <- rownames(precision_xxyy)
diag(parcor_xxyy) <- 1
stopifnot(all.equal(parcor_xxyy, parcor_xxyy2))
pc <- parcor_xxyy

diag(pc) <-  0
tol <- 1e-6
nzpc <- abs(pc) > tol
summary(pc[nzpc])
hist(pc[nzpc], 20)
summary(abs(pc[nzpc]))

#1 to 1 precision_xxyy and Partial Correlation
stopifnot(identical(abs(parcor_xxyy) > tol, abs(conc_xxyy) > tol))
#1 to 1 precision_xxyy and Partial Correlation
stopifnot(identical(abs(parcor_xxyy) > tol, abs(precision_xxyy) > tol))

# #Method 1 (Precision Parameterization -> Cholesky Decomposition -> Backsolve  )
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



#######Diagonal dominance(SPACE, Jasa, 2009, recipe)#########

# tmp <- as.matrix(get.adjacency(graph = cgpag, type = "both", attr = "weight"))
# is.symmetric.matrix(tmp)
# is.positive.definite(tmp)

# dd <- 1
# tmp <- as.matrix(get.adjacency(graph = cgpag, type = "upper", attr = "weight"))
# rs <- rowSums(abs(tmp))
# rescale <- ifelse(rs == 0, 1, 1.4*rs)
# #tmp <- t(t(tmp)/ (2*rescale))
# tmp <- t(sapply(1:nrow(tmp), function(i) tmp[i,] / rescale[i]))
# tmp[lower.tri(tmp)] <- tmp[upper.tri(tmp)]
# A <- tmp
# A <- ( tmp + t(tmp) )  / 2
# diag(A) <- dd
# is.symmetric.matrix(A)
# tmp <- A
# tmp[lower.tri(tmp, diag = T)] <- 0.0
# #diagonally dominant
# offdiag <- rowSums(abs(tmp)) 
# #stopifnot(max(offdiag) < dd)
# stopifnot(is.positive.definite(A))
# #Ainv <- solve(A) 
# #is.symmetric.matrix(Ainv) 
# # FALSE
# Ainv <- chol2inv(chol(A))
# stopifnot(is.symmetric.matrix(Ainv))
# kappa(A)
# kappa(Ainv)
# #numerical issue here
# #NegSqrtDiagAinv <- diag(x = sqrt(1/diag(Ainv)))
# #vsigma_xxyy <- NegSqrtDiagAinv %*% Ainv %*% NegSqrtDiagAinv
# #is.positive.definite(vsigma_xxyy)
# #FALSE 
# 
# #diagonal is all ones
# U <- (R + S)
# sigma_xxyy <- matrix(NA, U, U)
# for(i in 1:U) for(j in 1:U) sigma_xxyy[i,j] <- (Ainv[i,j] / sqrt(Ainv[i,i]*Ainv[j,j]))
# #all.equal(sigma_xxyy, vsigma_xxyy)
# # [1] TRUE
# stopifnot(all.equal(rep(1.0, R +S), diag(sigma_xxyy)))
# stopifnot(is.positive.definite(sigma_xxyy))
# #very stable condition number
# kappa(sigma_xxyy)
# 
# 
# # cor_xxyy <- cov2cor(sigma_xxyy)
# # kappa(cor_xxyy)
# # is.positive.definite(cor_xxyy)
# 
# #calculate partial correlations
# library("corpcor")
# parcor_xxyy2 <- cor2pcor(m = sigma_xxyy)
# #Spacemap way
# conc_xxyy <- solve(sigma_xxyy)
# summary(diag(conc_xxyy))
# parcor_xxyy <- conc2parcor(conc_xxyy)
# diag(parcor_xxyy) <- 1
# stopifnot(all.equal(parcor_xxyy, parcor_xxyy2))
# pc <- parcor_xxyy
# 
# 
# 
# 
# 
# diag(pc) <-  0
# nzpc <- abs(pc) > tol
# summary(pc[nzpc])
# hist(pc[nzpc], 20)
# summary(abs(pc[nzpc]))
# 
# 
# # #Define the X and Y index sets
# iy <- setdiff(1:R, xhubs_index)
# ix <- xset
# 
# #Validate the assumptions of the network topology
# tol <- 1e-6
# library(spacemap)
# #Y to Y
# nonZeroUpper(as.matrix(pc)[iy,iy], tol)
# #All X  to Y
# nonZeroWhole(as.matrix(pc)[ix,iy], tol)
# #X hubs to Y
# nonZeroWhole(as.matrix(pc)[xhubs_index,iy], tol)
# #27 X->Y edges with multiple hits (330 - 313)
# sum(degree(cgpag)[xhubs_index])
# #no edges from X background to Y
# nonZeroWhole(as.matrix(pc)[xbg_index,iy], tol)
# trueParCor <- list(xy = L_xy, yy = L_yy)
# 
# 
# 
# library(mvtnorm)
# dat <- rmvnorm(n = n, mean  = rep(0.0, ncol(sigma_xxyy)), sigma = sigma_xxyy)
# 
# 


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
