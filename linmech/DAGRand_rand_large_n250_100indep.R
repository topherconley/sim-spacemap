######DAGrand_rand
source("score_shd_09192012.txt") 
source("summary_functions.txt")
source("simpleRF_seq_02072013.txt")  
#######
load("dense_p500_adj_matrix.Rdata")
true.v=vstructures(true.dir)

n=250
p=nrow(true.ske)
step.max=1000
nrep=100
adjscore.indep.rand=array(0,c(p,p,nrep))
score.indep.rand=NULL
for(rep in 1:nrep){
  file.name=paste("dense_p500_",n,"_rep",rep,"_indep.Rdata",sep='')
   load(file.name)
   Y=Y.n
   n=nrow(Y)
   p=ncol(Y)
   rf.cutoff=runif(1,0,1)
   score.indep.rand[[rep]]=score_bic(Y, true.ske=true.ske,threshold=0, step=step.max,random.forest=TRUE,rf.cutoff=rf.cutoff, shuffle=TRUE)
   adjscore.indep.rand[,,rep]=score.indep.rand[[rep]]$adj.matrix
}

shd.indep.rand=score_shd(score.type="SHD",adjscore.indep.rand,step=step.max,threshold=0.5)
    adjshd.indep.rand=score_shd(score.type="SHD.adj",adjscore.indep.rand,step=step.max,threshold=0.5)
    source("avgCPSF_10032012_test.txt")
    avgCPSF.indep.rand=pair_CPSF(boot.adj=adjscore.indep.rand,threshold=0,selfreq.threshold=0,step=step.max)
    source("avgCPSF_largeNeigh_10192012.txt")
    LNCPSF.indep.rand=pair_CPSF(boot.adj=adjscore.indep.rand,threshold=0,selfreq.threshold=0,step=step.max)
    file.name=paste("DAGrand_rand_large_",n,"_indep.Rdata",sep='')
    save(adjscore.indep.rand,score.indep.rand, shd.indep.rand, adjshd.indep.rand, avgCPSF.indep.rand, LNCPSF.indep.rand,file=file.name)


####
result_skeleton(shd.indep.rand$adj.matrix,true.ske)
temp=result_vstruct(shd.indep.rand$adj.matrix,true.v) 
c(temp$total.num,temp$corr.num)

result_skeleton(adjshd.indep.rand$adj.matrix,true.ske)
temp=result_vstruct(adjshd.indep.rand$adj.matrix,true.v) 
c(temp$total.num,temp$corr.num)

result_skeleton(avgCPSF.indep.rand$adj.matrix,true.ske)
temp=result_vstruct(avgCPSF.indep.rand$adj.matrix,true.v) 
c(temp$total.num,temp$corr.num)

result_skeleton(LNCPSF.indep.rand$adj.matrix,true.ske)
temp=result_vstruct(LNCPSF.indep.rand$adj.matrix,true.v) 
c(temp$total.num,temp$corr.num)

