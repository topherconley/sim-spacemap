para_linear <- function(Y, adj_matrix){
##Y: n by p matrix, dependent variable
##adj_matrix: p by matrix, adjacency matrix
##estimate mean, variance by linear regression based on adjacency matrix
##return mean vector, variance vector, and coefficient matrix
##written in 10-12-2012

  n <- nrow(Y)
  p <- ncol(Y)
  if(nrow(adj_matrix) != ncol(adj_matrix))  
    stop("adjacency matrix is not square matrix!")
  if(nrow(adj_matrix) != p)
    stop("adjacency matrix does not have same dimension as Y")

  mu.est <- numeric(p)
  sigma.est <- numeric(p)
  beta.est <- matrix(0, p, p)

  mean.Y <- colMeans(Y)

  for(i in 1:p){
    par_index <- which(adj_matrix[,i] == 1)
    if(length(par_index) > 0){
      X <- matrix(0,n,(length(par_index) + 1))
      X[,1] <-1
      X[, 2: ncol(X)] <- Y[, par_index]

      temp <- solve(t(X) %*% X, t(X) %*% Y[,i])
      mu.est[i] <- temp[1]
      beta.est[i, par_index] <- temp[2: ncol(X)]
      sigma.est[i] <- sum((Y[, i] - X %*% temp) ^ 2) / (n - 1)
    }
    else{
      mu.est[i] <- mean.Y[i]
      sigma.est[i] <- sum((Y[, i] - mean.Y[i])^2) / (n - 1)
    }
  }
  return(list(mu.tilta = mu.est, sigma.est = sigma.est, beta.est = beta.est))
}

cal_order<- function(adj_matrix){
##calculate topological order of nodes in directed acyclic graph
##written in 10-12-2012

  p <- nrow(adj_matrix)
  if(nrow(adj_matrix) != ncol(adj_matrix))  
    stop("adjacency matrix is not square matrix!")

  node_order <- NULL
  ind_p <- numeric(p)

  csums <- colSums(adj_matrix)
  npar_index <- which(csums == 0)
  hpar_index <- which(csums != 0)
  ind_p[npar_index] <- 1

  node_order <- c(node_order, npar_index)
  current <- hpar_index
  while(length(current) > 0){
    j <- current[1]
    par_index <- which(adj_matrix[,j] == 1)
    if(all(ind_p[par_index] == 1)){
      node_order <- c(node_order,j)
      current <- current[-1]
      ind_p[j] <- 1
    }
    else{
      current <- current[-1]
      current <- c(current,j)
    }
  }

  return(node_order)
}

para_est <- function(mu, sigma, beta, adj_matrix, node_order){
##mu and sigma: mean and variance, p-length vector
##beta: p by p matrix, coefficient matrix
##estimate mean, covariance matrix based on adjacency matrix, topological order and estimated mean, variance and coefficients.
##written in 10-12-2012

  p <- nrow(adj_matrix)
  if(nrow(adj_matrix) != ncol(adj_matrix))  
    stop("adjacency matrix is not square matrix!")

  mu.est <- numeric( p)
  sigma.est <- matrix(0,p,p)
  B <- matrix(0,p,p)

  for(i in node_order){
    par_index <- which(adj_matrix[, i] == 1)
    if(length(par_index) == 0){
      B[i, i] <- 1
    }
    else{
      temp <- matrix(beta[i, par_index], 1, length(par_index))
      B[i, ] <- temp %*% B[par_index, ]
      B[i, i] <- B[i, i] + 1
    }
  }
  mu.est <- B %*% matrix(mu, p, 1)
  sigma.est <- B %*% diag(sigma) %*% t(B)

  return(list(mu.est = mu.est,sigma.est = sigma.est, B = B))
}

para_est2 <- function(mu, sigma, beta, adj_matrix, node_order){
##another version of para_est function
##written in 10-12-2012

  p <- nrow(adj_matrix)
  if(nrow(adj_matrix) != ncol(adj_matrix))  
    stop("adjacency matrix is not square matrix!")

  mu.est <- numeric(p)
  sigma.est <- matrix(0, p, p)
  B <- matrix(0, p, p)

  for(i in node_order){
    par_index <- which(adj_matrix[, i] == 1)
    if(length(par_index) == 0){
      B[i, i] <- 1
      mu.est[i] <- mu[i]
    }
    else{
      temp <- matrix(beta[i, par_index], 1, length(par_index))
      B[i, ] <- temp %*% B[par_index, ]
      B[i, i] <- B[i, i] + 1
      mu.est[i] <- mu[i] + sum(beta[i, par_index] * mu.est[par_index])
    }
  }
  sigma.est <- B %*% diag(sigma) %*% t(B)

  return(list(mu.est = mu.est, sigma.est = sigma.est, B = B))
}






