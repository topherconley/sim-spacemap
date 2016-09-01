#not found
#source("~/repos/sim-spacemap/linmech/space_score_functions_10082011.txt")  ##space.shd, and space.shd.adj
load("~/tmp/dense_p500_adj_matrix.Rdata")
true.v=vstructures(true.dir)
############################# part II: generate Guassin linear mechanism based on pancr.adj 
p=nrow(true.dir)
n=100

nrep=10
source("~/repos/sim-spacemap/linmech/network_gene_function_10172012.R")
data_generation<-function(n,p,pancr.adj){

  Y.n=NULL
  panc.simu.Data=gene.pancr.label(n, pancr.adj,SN=runif(nrow(pancr.adj), min=0.5, max=1.5))

  Y=panc.simu.Data$data
  Y.gm=apply(Y, 2, mean)
  Y.gsd=apply(Y, 2, sd)
  Y.n=(Y-matrix(Y.gm, n, p, byrow=T))/matrix(Y.gsd, n, p, byrow=T)
  return(Y.n=Y.n)
}


##look at gene.pancr.label

set.seed(2000)
seeds=sample(2:1000,10)
for(rep in 1:nrep){
  print(rep)
  set.seed(seeds[rep])
  Y.n=data_generation(n,p,true.dir)
  file.name=paste("dense_p500_",n,"_rep",rep,".Rdata",sep='')
  save(Y.n,file=file.name)
}
