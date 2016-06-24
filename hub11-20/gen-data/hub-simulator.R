#'@title Generate model 1 simulation for spacemap development. 
#'@param mdegree the lower bound for the degree of master predictors
#'@param pnone the number of X variables independent of Y.
#'@param tol the convergence tolerance
#'@param randxy Logical indicating whether non-zero degree (X,Y) are randomly assigned.
#'@return A list of model parameters including the true adjacency matrices for (X,Y).
simSpaceMapModel1 <- function(ParCor.true = NULL, SigmaS = NULL, randxy = FALSE,
                              mdegree = 9, pnone = 20L, tol = 1e-6) {  
  
  #Check parameter types
  stopifnot(is.integer(pnone), is.integer(mdegree))
  
  if (!is.null(ParCor.true) & is.null(SigmaS)) {
    SigmaSInv <- (-1)*ParCor.true
    #Put partial correlation and Sig.Inv on the same scale
    diag(SigmaSInv) <- 1
    SigmaS <- solve(SigmaSInv)
  } else if (is.null(ParCor.true) & !is.null(SigmaS)) {
    library(spacemap)
    ParCor.true <- cov2ParCor(SigmaS)
  } else {
    message("Must specify either SigmaS or ParCor, but not both.")
    return(NULL)
  }
  
  #degree distribution calculation
  #theory tells us the next two lines are equivalent. 
  #true.adj <- abs(SigmaSInv) > tol
  true.adj <- abs(ParCor.true) > tol
  
  #do not count self-adjacency 
  diag(true.adj) <- FALSE
  degree <- rowSums(true.adj)
  #Simulation may already have noise variables
  YPoorIndex <- which(degree == 0)
  
  if (pnone > 0L) {
    library(Matrix)
    SigmaG <- as.matrix(bdiag(SigmaS, diag(pnone)))    ### add "pnone" number of independent N(0,1) variables 
    #update degree distribution with poor predictors. 
    degree <- c(degree, vector(mode = "numeric", length = pnone))
  } else {
    SigmaG <- SigmaS
  }
  
  XPoorIndex <- setdiff(which(degree == 0), YPoorIndex)
  XMasterIndex <- which(degree >= mdegree)
  Xindex <- union(XMasterIndex, XPoorIndex)
  Yindex <- union(which(degree >= 1 & degree < mdegree), YPoorIndex)
  #Completely separate
  stopifnot(length(intersect(Xindex,Yindex)) == 0)
  
  #randomly assign X nodes and Y nodes that do not have degree zero.
  #in the case of the hub (p=500) simulation, this is all non-noise added
  #variables. 
  if(randxy) {
    remainingNodes <- seq_along(degree)
    remainingNodes <- setdiff(remainingNodes, XPoorIndex)
    XMasterIndex <- sample(remainingNodes, length(XMasterIndex))
    Xindex <- union(XMasterIndex, XPoorIndex)
    remainingNodes <- setdiff(remainingNodes, XMasterIndex)
    stopifnot(length(Yindex) == length(remainingNodes))
    Yindex <- remainingNodes
  }
  #True Parameters
  trueXY <- rbind(ParCor.true[XMasterIndex,Yindex],
                  matrix(0,nrow = pnone, ncol = length(Yindex)))
  trueYY <- ParCor.true[Yindex,Yindex]
  
  list(Sigma = SigmaG, trueYY = trueYY, trueXY = trueXY, 
       Xindex = Xindex, XPoorIndex = XPoorIndex, XMasterIndex = XMasterIndex, 
       Yindex = Yindex)
}

#'@title Generate random date according to Model 1
#'@note A companion function to simSpaceMapModel1(...)
#'@param n the sample size
#'@param model the model object returned from simSpaceMapModel1(...)
#'@return A list of data (X,Y) according to Model 1.
simSpaceMapData1 <- function(n = 250, model) {
  #generate simulation
  library(mvtnorm)
  dat <- rmvnorm(n = n, mean  = rep(0.0, ncol(model$Sigma)), sigma = model$Sigma)
  list(XY = dat, X = dat[,model$Xindex], Y = dat[,model$Yindex], 
       Xindex = model$Xindex, XPoorIndex = model$XPoorIndex, XMasterIndex = model$XMasterIndex, 
       Yindex = model$Yindex)
}
