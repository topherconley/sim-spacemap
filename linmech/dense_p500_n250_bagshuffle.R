
#============================== Sample code to run on Gauss =============================#
source("score_bic_11152012_dposv.txt") ##score.bic
source("score_shd_09192012.txt")
source("summary_functions.txt")
source("score_halfshd_05082012.txt")
source("vstruct_score_09102012.txt")

#######
load("dense_p500_adj_matrix.Rdata")

n=250
p=nrow(true.ske)
step.max=1000
n.B=100              #number of bootstrap 
boot.method ="hc"    ##method to get bootstrapped networks
seed.u=100           ##seed for bootstrap resample

nrep=100


#========== Steup for running on Gauss ==========#


"ppaste" <- function(...){paste(...,sep="")}

args <- commandArgs(TRUE)

cat(ppaste("Command-line arguments:\n"))
print(args)

###################
sim_start <- 0
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  sinkit <- FALSE
} else {
  # SLURM can use either 0- or 1-indexing...
  sinkit <- TRUE
  sim_num <- sim_start + as.numeric(args[1])
}

i <- sim_num
sinkfile <- paste("Output",i,".txt",sep="")

cat(paste("\nAnalyzing dataset number ",i,"...\n\n",sep=""))

gauss <- TRUE

if (sinkit){
  cat(paste("Sinking output to: ",sinkfile,"\n",sep=""))
  sink(sinkfile)
}  



#========== Run the simulation study ==========#
rep=i
  file.name=paste("dense_p500_",n,"_rep",rep,"_indep.Rdata",sep='')
  load(file.name)

  adjscore.baggingShuff=array(0,c(p,p,n.B))

   Y=Y.n
   n=nrow(Y)
   p=ncol(Y)

   temp=NULL
  set.seed(4)
  for(b in 1:n.B){
    y.index=sample(1:n,n,replace=TRUE)
    Y.B=Y[y.index,]

    temp[[b]]=score_bic(Y.B, threshold=0, step=step.max,standardize=TRUE, shuffle=TRUE)
    adjscore.baggingShuff[,,b]=temp[[b]]$adj.matrix
   }
   
score.bagshf=score_bic(Y, threshold=0, step=step.max,standardize=TRUE, shuffle=TRUE)

    shd.bagshf=score_shd(score.type="SHD",adjscore.baggingShuff,step=step.max,threshold=0.5)
    adjshd.bagshf=score_shd(score.type="SHD.adj",adjscore.baggingShuff,step=step.max,threshold=0.5)
    source("avgCPSF_10102012_test.txt")
    avgCPSF.bagshf=pair_CPSF(boot.adj=adjscore.baggingShuff,threshold=0,selfreq.threshold=0,step=step.max)
    source("avgCPSF_largeNeigh_10192012.txt")
    LNCPSF.bagshf=pair_CPSF(boot.adj=adjscore.baggingShuff,threshold=0,selfreq.threshold=0,step=step.max)
    thalf.score.bagshf=score_halfshd(score.type="SHD_thalf",boot.adj=adjscore.baggingShuff,step=step.max,blacklist=NULL,whitelist=NULL,threshold=0.5,print=TRUE)
    half.score.bagshf=score_halfshd(score.type="SHD_half",boot.adj=adjscore.baggingShuff,step=step.max,blacklist=NULL,whitelist=NULL,threshold=0.5,print=TRUE)
    vscore.bagshf=score_vstruct(boot.adj=adjscore.baggingShuff, step = step.max, blacklist = NULL, whitelist = NULL, threshold = 0, skeleton_unit=1, vstruct_unit=0, v.factor=0, print = TRUE)

      file.name=paste("bagshf_dense_p500_n_",n,"_rep",rep,".Rdata",sep='')

    save(temp,adjscore.baggingShuff,score.bagshf,shd.bagshf,adjshd.bagshf,avgCPSF.bagshf,LNCPSF.bagshf,thalf.score.bagshf,half.score.bagshf,vscore.bagshf,file=file.name)

