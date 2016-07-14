
library(igraph)

gen_powlaw_module <- function(R, S, id) { 
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
  plf <- power.law.fit(x = degree(pag))
  
  #Identify hubs for X and Y
  hub_pag <- hub.score(pag, scale = TRUE, weights=NULL, options = igraph.arpack.default)
  #plot(hub_pag$vector, deg_pag, pch = 19)
  #top hubs split in half among X and Y
  tophubs_index <- sort(hub_pag$vector, index.return = TRUE, decreasing = TRUE)$ix[1:4]
  xhubs_index <- sample(tophubs_index, size = ceiling(length(tophubs_index)/2))
  yhubs_index <- setdiff(tophubs_index, xhubs_index)
  #hist(degree(pag)[xhubs_index])
  #hist(degree(pag)[yhubs_index])
  
  #How many x-x edges exist? 
  #very few edges (14), out of the whole graph (1185) and the total X hub edges (341), okay to efffectively delete
  xx_edges <- E(pag)[ xhubs_index %--% xhubs_index]
  #delete these edges from the graph as we are trying to create no geographic (genomic) correlation between true hubs
  cgpag <- delete_edges(pag, xx_edges)
  #14 have been deleted
  #ecount(cgpag)
  #new degrees for X hubs
  #degree(cgpag, v = xhubs_index)

  #Still remains a power-law with alpha  = 2.59 in the range expected for biological networks
  noxx_plf <- power.law.fit(x = degree(cgpag))
  #noxx_plf$alpha
  
  #the graph  has no disconnected components among the original power-law distribution
  stopifnot(count_components(cgpag) == 1)
  
  #generate weights among X--Y and Y--Y edges 
  wgts <- runif(n =  ecount(cgpag), min = 0.2, max =  1)
  #wgts <- ifelse(test = rbinom(n = ecount(cgpag), size = 1, prob = 0.5), wgts, -1*wgts)
  cgpag <- set_edge_attr(graph = cgpag, name = "weight", value = wgts)
  
  #ceb <- cluster_edge_betweenness(graph = cgpag, directed = FALSE)
  #modularity(x = cgpag, membership = ceb$membership)
  
  #add additional 200 vertices that have some correlation structure wth each hub
  cgpag <- add_vertices(graph = cgpag, nv = S)
  xbg_index <- (R + 1):(R + S)
  
  
  #For each X hub create a Erdos-Renyi random graph with 5 other x vertices that have no edges to y vertices.
  start <- R + 1
  ostart <- start
  for (xhub in xhubs_index) { 
    xid <- c(xhub, start:(start + 4))
    tmpg <- sample_gnp(n = 6, p  = 0.3)
    tmpg <- set_vertex_attr(graph = tmpg, name = "name", value = xid)
    cgpag <- add_edges(cgpag, edges = as.vector(t(get.edgelist(tmpg))))
    start <- start + 5
  }
  
  xset1 <- c(xhubs_index, ostart:(start - 1))
  xset1g <- subgraph.edges(cgpag, E(cgpag)[xset1 %--% xset1])
  
  #For every 5 other background X's, create an Erdos-Reny random graph in the same fashion, but do not include an X hub. 
  for(xbs in seq(from = start, to = (R + S), by = 5)) {
    xid <- xbs:(xbs + 4)
    tmpg <- sample_gnp(n = 5, p  = 0.3)
    tmpg <- set_vertex_attr(graph = tmpg, name = "name", value = xid)
    cgpag <- add_edges(cgpag, edges = as.vector(t(get.edgelist(tmpg))))
  }
  
  xset2 <- start:(R + S)
  xset2g <- subgraph.edges(cgpag, E(cgpag)[xset2 %--% xset2])
  #no overlap 
  stopifnot(length(intersect(xset1, xset2)) == 0)
  
  #generate weights among X--X edges (all positive)
  xset <- unique(c(xset1, xset2))
  wgts <- runif(n =  ecount(xset1g) + ecount(xset2g), min = 0.2, max = 1)
  xset <- unique(c(xset1, xset2))
  cgpag <- set_edge_attr(graph = cgpag, name = "weight", index = E(cgpag)[xset %--% xset], value = wgts)
  
  #Label types for visualization
  cgpag <- set_vertex_attr(cgpag, name = "type", index = xhubs_index, value = "xhub")
  cgpag <- set_vertex_attr(cgpag, name = "type", index = yhubs_index, value = "yhub")
  cgpag <- set_vertex_attr(cgpag, name = "type", index = xbg_index, value = "x")
  ybg_index <- setdiff(1:R, c(xhubs_index, yhubs_index))
  cgpag <- set_vertex_attr(cgpag, name = "type", index = ybg_index, value = "y")
  
  #unique vertex id
  #set_vertex_attr(cgpag, name = "name", index = V(cgpag), value = paste(id, V(cgpag), sep = ":"))
  V(cgpag)$name <- paste(id, V(cgpag), sep = ":")
  cgpag <- set_vertex_attr(cgpag, name = "module_id", index = V(cgpag), value = id)
  list(mod = cgpag, 
       xhubs_index = paste(id, xhubs_index, sep = ":"), 
       yhubs_index = paste(id, yhubs_index, sep = ":"),
       xbg_index  = paste(id, xbg_index, sep = ":"),
       ybg_index = paste(id, ybg_index, sep = ":"))
}


#number of nodes with the potential for an edge
R <- 100
#number of X nodes with no edge
S <- 40

#reproducibility
set.seed(26225)

#generate 5 modules and stitch together information
plgmods <- lapply(1:5, function(id) gen_powlaw_module(R = R, S = S, id = id))
gmods <- lapply(plgmods, function(x) x$mod)
cgpag <- do.call(what = union, args = gmods)
xhubs_index <- Reduce('c', lapply(plgmods, function(x) x$xhubs_index))
xbg_index <- Reduce('c', lapply(plgmods, function(x) x$xbg_index))
yhubs_index <- Reduce('c', lapply(plgmods, function(x) x$yhubs_index))
ybg_index <- Reduce('c', lapply(plgmods, function(x) x$ybg_index))

fit_power_law(x = degree(cgpag))


#set_edge_attr(graph = cgpag, name = "")
vertex_attr_names(cgpag)
edge_attr_names(cgpag)

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
#Precision <- D - A
Precision <- laplacian_matrix(graph = cgpag, normalized = TRUE, sparse = FALSE)
#verify properties of laplacian
library(matrixcalc)
is.positive.definite(Precision)
is.positive.semi.definite(Precision)
is.symmetric.matrix(Precision)
#Add small positive values to diagonal
diag(Precision) <- diag(Precision) + 0.2
stopifnot(is.positive.definite(Precision))
#Define the X and Y index sets
iy <- setdiff(1:R, xhubs_index)
ix <- xset
#assure non-overlapping
stopifnot(intersect(iy, ix) == integer(0))
#Broken down by x and y blocks
L_xy <- Precision[ix,iy]
L_yy <- Precision[iy,iy]
is.positive.definite(L_yy)

#Validate the assumptions of the network topology
tol <- 1e-6
library(spacemap)
#Y to Y
nonZeroUpper(as.matrix(Precision)[iy,iy], tol)
#All X  to Y
nonZeroWhole(as.matrix(Precision)[ix,iy], tol)
#X hubs to Y
nonZeroWhole(as.matrix(Precision)[xhubs_index,iy], tol)
#X->Y edges with multiple hits
sum(degree(cgpag)[xhubs_index])
#no edges from X background to Y
nonZeroWhole(as.matrix(Precision)[xbg_index,iy], tol)


#Convert Precision Matrix to Sigma Matrix
sigma_xxyy <- chol2inv(chol(Precision))  
stopifnot(is.positive.definite(sigma_xxyy))
summary(diag(sigma_xxyy))
#generate simulation
library(mvtnorm)
n <- 250
dat <- rmvnorm(n = n, mean  = rep(0.0, ncol(sigma_xxyy)), sigma = sigma_xxyy)

summary(apply(dat, 2, sd))

kappa(sigma_xxyy)
kappa(Precision)


library("corpcor")
parcor_xxyy2 <- cor2pcor(m = sigma_xxyy)
#Spacemap way
conc_xxyy <- solve(sigma_xxyy)
summary(diag(conc_xxyy))
parcor_xxyy <- conc2parcor(conc_xxyy)
diag(parcor_xxyy) <- 1
stopifnot(all.equal(parcor_xxyy, parcor_xxyy2))
pc <- parcor_xxyy

diag(pc) <-  0
tol <- 1e-6
nzpc <- abs(pc) > tol
summary(pc[nzpc])
hist(pc[nzpc], 20)
summary(abs(pc[nzpc]))

#1 to 1 Precision and Partial Correlation
identical(abs(parcor_xxyy) > tol, abs(conc_xxyy) > tol)

# #Method 1 (Precision Parameterization -> Cholesky Decomposition -> Backsolve  )
# R <- chol(Precision)
# Lqr <- qr(x = Precision)
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
