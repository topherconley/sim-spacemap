

xx_banded_net <- function(rho, px) { 
  #banded correlation matrix
  library(Matrix)
  #create list of matrices
  block <- function(r) { 
    xx_sigma <- matrix(rho[r], px[r],px[r])
    for (i  in 1:px)
      for (j in 1:px) { 
        xx_sigma[i,j] <-  xx_sigma[i,j]^(abs(i - j))
      }
    xx_sigma
  }
  lmat <- lapply(seq_along(rho), block)
  xx_sigma <- as.matrix(bdiag(lmat))
  #get inverse
  xx_sigma_inv <- solve(xx_sigma)
  #make symmetric after small round-off errors in inverse
  xx_sigma_inv_sym <- (xx_sigma_inv + t(xx_sigma_inv))/2
  #drop the small error
  xx_sigma_inv_sym[abs(xx_sigma_inv_sym) < .Machine$double.eps] <- 0
  #assure that the differences are negligible after symmetry step
  stopifnot(all.equal(xx_sigma_inv_sym, xx_sigma_inv))
  #assure precision/sigma are p.d.
  library(matrixcalc)
  stopifnot(is.positive.definite(xx_sigma_inv_sym))
  stopifnot(is.positive.definite(xx_sigma))

  xx_adj <- (abs(xx_sigma_inv_sym) > 0) + 0
  diag(xx_adj) <- 0
  #Summary
  nxx <- sum(abs(xx_sigma_inv_sym[upper.tri(xx_sigma_inv_sym)]) > .Machine$double.eps)
  message("Number of conditionally dependent edges between X--X: ", nxx)
  message("Min. diagonal of X--X precision: ", min(diag(xx_sigma_inv_sym)))
  message("Max. diagonal of X--X precision: ", max(diag(xx_sigma_inv_sym)))
  
#   library(Matrix)
#   image(abs(band(xx_sigma_inv_sym, -100, 100)) > .Machine$double.eps )

  list(xx_precision = xx_sigma_inv_sym, xx_sigma = xx_sigma, 
       xx_adj = xx_adj)
} 

set.seed(94627)
rho <- sample(seq(0.4, 0.6,length = 10))
px <- rep(50, 10)
xx_net <- xx_banded_net(rho = rho, px = px)
library(igraph)

library(Matrix)
image(abs(band(xx_net$xx_adj, -100, 100)) > .Machine$double.eps )

ixx <- igraph::graph_from_adjacency_matrix(adjmatrix = xx_net$xx_adj, mode  = "undirected")
#igraph::count_components(ixx)

igraph::components(ixx)

hist(degree(ixx))
