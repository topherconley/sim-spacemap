#reproducibility of data
set.seed(5317)

#################################
# Distinguising Simulation Paramters of hub-11-20-250
masterDegree <- 11L ## minimum out-degree of master predictors 
indepX <- 20L ## number of independent X variables. 
n <- 250L ## sample size
randxy = FALSE
tol <- 1e-6 ##convergence tolerance

##get the template partial correlation matrix from space package example
library(space) 
data(spaceSimu)
ParCor.true=spaceSimu$ParCor.true

#degree distribution is the same with orders of magnitude precision
identical(abs(ParCor.true) > 1e-2, abs(ParCor.true) > 1e-6)

#function to generate model
source("~/repos/sim-spacemap/hub11-20/gen-data/hub-simulator.R")
## generate true covariance matrix, and the X, Y indices 
#evaluate model parameters
model.c <- simSpaceMapModel1(ParCor.true, mdegree = masterDegree, 
                             randxy = randxy,
                             pnone = indepX, tol = tol)
# Data generation with the correlation matrix instead of the covariance: 
# The purpose is to make the diagonal elements of the precision matrix to be different from 
# 1 so that the estimates, sig.fit, are not set to the truth by default. 
model.c$Sigma <- cov2cor(model.c$Sigma)
# diag(solve(model.c$Sigma))
# summary(diag(solve(cov2cor(model.c$Sigma))))
# identical(solve(model.c$Sigma) > tol, solve(cov2cor(model.c$Sigma)) > tol)

#useful for repeated evaluation of performance against truth
trueParCor <- list(xy = model.c$trueXY, yy = model.c$trueYY)

##Make 100 datasets
datasets_path <- "~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250"
ndatasets <- 100L
for(dataid in seq_len(ndatasets)) { 
  ## generate data
  data.c=simSpaceMapData1(n=n, model=model.c)
  
  #standardize the data
  data.c$XY <- scale(data.c$XY)
  data.c$X <- scale(data.c$X)
  data.c$Y <- scale(data.c$Y)
  
    
  #list of objects to save
  data_out <- list(id = dataid, 
                   XY = data.c$XY, 
                   X = data.c$X,
                   Y = data.c$Y,
                   Xindex = data.c$Xindex, 
                   Yindex = data.c$Yindex, 
                   trueParCor = trueParCor, 
                   tol = 1e-3)
  file_out <- file.path(datasets_path, paste0("hub11-20-250-dataid_", sprintf("%03d",dataid), ".rds"))
  saveRDS(data_out, file = file_out)
}

#The particular model 1 parameter settings generated a true graph comprised of: 
length(data.c$Yindex) ## number of Y nodes
length(data.c$Xindex) ## number of X nodes
length(data.c$Xindex) - indepX ## number of master predictors.


# library(spacemap)
# system.time({tmp <- space.joint(Y.m = data.c$XY, lam1 = 85, lam2 = 0, iter = 3, tol = 1e-3)})
# 
# system.time({tmp <- spacemap(Y.m = data.c$Y, X.m = data.c$X, slasso = 75, sridge = 0, rlasso = 30, rgroup = 20, tol = 1e-3)})

