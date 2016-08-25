source("score_GIC_04202013_dposv.txt") ##score.bic
source("score_shd_09192012.txt")
source("summary_functions.txt")

#######
load("dense_p500_adj_matrix.Rdata")
true.v=vstructures(true.dir)

n=100
p=nrow(true.ske)
step.max=1000
n.B=100
nrep=10

res.bagshf.gic=matrix(0,2,nrep)
res.bagshf.mgic=matrix(0,2,nrep)
vs.bagshf.gic=matrix(0,2,nrep)
vs.bagshf.mgic=matrix(0,2,nrep)

for(rep in 1:nrep){
   adjscore.gic=array(0,c(p,p,n.B))
   adjscore.mgic=array(0,c(p,p,n.B))

   file.name=paste("dense_p500_",n,"_rep",rep,".Rdata",sep='')
   load(file.name)
   Y=Y.n
   n=nrow(Y)
   p=ncol(Y)

   set.seed(4)
   for(b in 1:n.B){
     y.index=sample(1:n,n,replace=TRUE)
     Y.B=Y[y.index,]

     temp=score_gic(Y.B, threshold=0, step=step.max,score.type="GIC",standardize=TRUE, shuffle=TRUE)
     adjscore.gic[,,b]=temp$adj.matrix

     temp2=score_gic(Y.B, threshold=0, step=step.max,score.type="MGIC",standardize=TRUE, shuffle=TRUE)
     adjscore.mgic[,,b]=temp2$adj.matrix
   }

   score.gic=score_gic(Y, threshold=0, step=step.max,score.type="GIC", standardize=TRUE, shuffle=TRUE)
   score.mgic=score_gic(Y, threshold=0, step=step.max,score.type="MGIC", standardize=TRUE, shuffle=TRUE)


shd.gic=score_shd(score.type="SHD",adjscore.gic,step=step.max,threshold=0.5)
shd.mgic=score_shd(score.type="SHD",adjscore.mgic,step=step.max,threshold=0.5)  

  adjshd.gic=score_shd(score.type="SHD.adj",adjscore.gic,step=step.max,threshold=0.5)

adjshd.mgic=score_shd(score.type="SHD.adj",adjscore.mgic,step=step.max,threshold=0.5)

    source("avgCPSF_10032012_test.txt")
    avgCPSF.gic=pair_CPSF(boot.adj=adjscore.gic,threshold=0,selfreq.threshold=0,step=step.max)

avgCPSF.mgic=pair_CPSF(boot.adj=adjscore.mgic,threshold=0,selfreq.threshold=0,step=step.max)
    source("avgCPSF_largeNeigh_10192012.txt")
    LNCPSF.gic=pair_CPSF(boot.adj=adjscore.gic,threshold=0,selfreq.threshold=0,step=step.max)

LNCPSF.mgic=pair_CPSF(boot.adj=adjscore.mgic,threshold=0,selfreq.threshold=0,step=step.max)

    file.name=paste("gic_mgic_dense_n",n,"_rep",rep,".Rdata",sep='')
    save(adjscore.gic,adjscore.mgic,score.gic,score.mgic,shd.gic, shd.mgic, adjshd.gic, adjshd.mgic, avgCPSF.gic, avgCPSF.mgic, LNCPSF.gic, LNCPSF.mgic,file=file.name)
}

#########
res.bagshf.bic=matrix(0,2,nrep)
res.bagshf.shd=matrix(0,2,nrep)
res.bagshf.adjshd=matrix(0,2,nrep)
res.bagshf.avgCPSF=matrix(0,2,nrep)
res.bagshf.LNCPSF=matrix(0,2,nrep)

vs.bagshf.bic=matrix(0,2,nrep)
vs.bagshf.shd=matrix(0,2,nrep)
vs.bagshf.shdadj=matrix(0,2,nrep)
vs.bagshf.avgCPSF=matrix(0,2,nrep)
vs.bagshf.LNCPSF=matrix(0,2,nrep)

for(rep in 1:nrep){
    file.name=paste("gic_mgic_dense_n",n,"_rep",rep,".Rdata",sep='')
  load(file.name)
  res.bagshf.bic[,rep]=result_skeleton(score.gic$adj.matrix,true.ske)
  res.bagshf.shd[,rep]=result_skeleton(shd.gic$adj.matrix,true.ske)
  res.bagshf.adjshd[,rep]=result_skeleton(adjshd.gic$adj.matrix,true.ske)
  res.bagshf.avgCPSF[,rep]=result_skeleton(avgCPSF.gic$adj.matrix,true.ske)
  res.bagshf.LNCPSF[,rep]=result_skeleton(LNCPSF.gic$adj.matrix,true.ske)

   temp=result_vstruct(score.gic$adj.matrix,true.v) 
   vs.bagshf.bic[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(shd.gic$adj.matrix,true.v) 
   vs.bagshf.shd[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(adjshd.gic$adj.matrix,true.v) 
   vs.bagshf.shdadj[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(avgCPSF.gic$adj.matrix,true.v) 
   vs.bagshf.avgCPSF[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(LNCPSF.gic$adj.matrix,true.v) 
   vs.bagshf.LNCPSF[,rep]=c(temp$total.num,temp$corr.num)
}

avgske.bic=rowMeans(res.bagshf.bic)
avgvs.bic=rowMeans(vs.bagshf.bic)
avgske.shd=rowMeans(res.bagshf.shd)
avgske.shdadj=rowMeans(res.bagshf.adjshd)
avgske.avgCPSF=rowMeans(res.bagshf.avgCPSF)
avgske.LNCPSF=rowMeans(res.bagshf.LNCPSF)
avgvs.shd=rowMeans(vs.bagshf.shd)
avgvs.shdadj=rowMeans(vs.bagshf.shdadj)
avgvs.avgCPSF=rowMeans(vs.bagshf.avgCPSF)
avgvs.LNCPSF=rowMeans(vs.bagshf.LNCPSF)

avg_table=matrix(0,5,9)
avg_table[,1]=c(0,0,0,0,0)
avg_table[,2]=c(avgske.bic[2],avgske.shd[2],avgske.shdadj[2],avgske.avgCPSF[2],avgske.LNCPSF[2])

avg_table[,3]=c(apply(res.bagshf.bic,1,sd)[2],apply(res.bagshf.shd,1,sd)[2],apply(res.bagshf.adjshd,1,sd)[2],apply(res.bagshf.avgCPSF,1,sd)[2],apply(res.bagshf.LNCPSF,1,sd)[2])

avg_table[,4]=c(avgske.bic[1],avgske.shd[1],avgske.shdadj[1],avgske.avgCPSF[1],avgske.LNCPSF[1])

avg_table[,5]=c(apply(res.bagshf.bic,1,sd)[1],apply(res.bagshf.shd,1,sd)[1],apply(res.bagshf.adjshd,1,sd)[1],apply(res.bagshf.avgCPSF,1,sd)[1],apply(res.bagshf.LNCPSF,1,sd)[1])

avg_table[,6]=c(avgvs.bic[2],avgvs.shd[2],avgvs.shdadj[2],avgvs.avgCPSF[2],avgvs.LNCPSF[2])
avg_table[,7]=c(apply(vs.bagshf.bic,1,sd)[2],apply(vs.bagshf.shd,1,sd)[2],apply(vs.bagshf.shdadj,1,sd)[2],apply(vs.bagshf.avgCPSF,1,sd)[2],apply(vs.bagshf.LNCPSF,1,sd)[2])
avg_table[,8]=c(avgvs.bic[2],avgvs.shd[1],avgvs.shdadj[1],avgvs.avgCPSF[1],avgvs.LNCPSF[1])

avg_table[,9]=c(apply(vs.bagshf.bic,1,sd)[1],apply(vs.bagshf.shd,1,sd)[1],apply(vs.bagshf.shdadj,1,sd)[1],apply(vs.bagshf.avgCPSF,1,sd)[1],apply(vs.bagshf.LNCPSF,1,sd)[1])

ll=nrow(avg_table)
for(i in 1:ll){
  tmp=paste(i,"&",avg_table[i,1],"&",avg_table[i,2],"(",round(avg_table[i,3],2),")","&",avg_table[i,4],"(",round(avg_table[i,5],2),")","&",avg_table[i,6],"(",round(avg_table[i,7],2),")","&",avg_table[i,8],"(",round(avg_table[i,9],2),")","\\",sep='')
  print(tmp)
}


##########
res.bagshf.bic=matrix(0,2,nrep)
res.bagshf.shd=matrix(0,2,nrep)
res.bagshf.adjshd=matrix(0,2,nrep)
res.bagshf.avgCPSF=matrix(0,2,nrep)
res.bagshf.LNCPSF=matrix(0,2,nrep)

vs.bagshf.bic=matrix(0,2,nrep)
vs.bagshf.shd=matrix(0,2,nrep)
vs.bagshf.shdadj=matrix(0,2,nrep)
vs.bagshf.avgCPSF=matrix(0,2,nrep)
vs.bagshf.LNCPSF=matrix(0,2,nrep)

for(rep in 1:nrep){
    file.name=paste("gic_mgic_dense_n",n,"_rep",rep,".Rdata",sep='')
  load(file.name)
  res.bagshf.bic[,rep]=result_skeleton(score.mgic$adj.matrix,true.ske)
  res.bagshf.shd[,rep]=result_skeleton(shd.mgic$adj.matrix,true.ske)
  res.bagshf.adjshd[,rep]=result_skeleton(adjshd.mgic$adj.matrix,true.ske)
  res.bagshf.avgCPSF[,rep]=result_skeleton(avgCPSF.mgic$adj.matrix,true.ske)
  res.bagshf.LNCPSF[,rep]=result_skeleton(LNCPSF.mgic$adj.matrix,true.ske)

   temp=result_vstruct(score.mgic$adj.matrix,true.v) 
   vs.bagshf.bic[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(shd.mgic$adj.matrix,true.v) 
   vs.bagshf.shd[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(adjshd.mgic$adj.matrix,true.v) 
   vs.bagshf.shdadj[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(avgCPSF.mgic$adj.matrix,true.v) 
   vs.bagshf.avgCPSF[,rep]=c(temp$total.num,temp$corr.num)

   temp=result_vstruct(LNCPSF.mgic$adj.matrix,true.v) 
   vs.bagshf.LNCPSF[,rep]=c(temp$total.num,temp$corr.num)
}

avgske.bic=rowMeans(res.bagshf.bic)
avgvs.bic=rowMeans(vs.bagshf.bic)
avgske.shd=rowMeans(res.bagshf.shd)
avgske.shdadj=rowMeans(res.bagshf.adjshd)
avgske.avgCPSF=rowMeans(res.bagshf.avgCPSF)
avgske.LNCPSF=rowMeans(res.bagshf.LNCPSF)
avgvs.shd=rowMeans(vs.bagshf.shd)
avgvs.shdadj=rowMeans(vs.bagshf.shdadj)
avgvs.avgCPSF=rowMeans(vs.bagshf.avgCPSF)
avgvs.LNCPSF=rowMeans(vs.bagshf.LNCPSF)

avg_table=matrix(0,5,9)
avg_table[,1]=c(0,0,0,0,0)
avg_table[,2]=c(avgske.bic[2],avgske.shd[2],avgske.shdadj[2],avgske.avgCPSF[2],avgske.LNCPSF[2])

avg_table[,3]=c(apply(res.bagshf.bic,1,sd)[2],apply(res.bagshf.shd,1,sd)[2],apply(res.bagshf.adjshd,1,sd)[2],apply(res.bagshf.avgCPSF,1,sd)[2],apply(res.bagshf.LNCPSF,1,sd)[2])

avg_table[,4]=c(avgske.bic[1],avgske.shd[1],avgske.shdadj[1],avgske.avgCPSF[1],avgske.LNCPSF[1])

avg_table[,5]=c(apply(res.bagshf.bic,1,sd)[1],apply(res.bagshf.shd,1,sd)[1],apply(res.bagshf.adjshd,1,sd)[1],apply(res.bagshf.avgCPSF,1,sd)[1],apply(res.bagshf.LNCPSF,1,sd)[1])

avg_table[,6]=c(avgvs.bic[2],avgvs.shd[2],avgvs.shdadj[2],avgvs.avgCPSF[2],avgvs.LNCPSF[2])
avg_table[,7]=c(apply(vs.bagshf.bic,1,sd)[2],apply(vs.bagshf.shd,1,sd)[2],apply(vs.bagshf.shdadj,1,sd)[2],apply(vs.bagshf.avgCPSF,1,sd)[2],apply(vs.bagshf.LNCPSF,1,sd)[2])
avg_table[,8]=c(avgvs.bic[2],avgvs.shd[1],avgvs.shdadj[1],avgvs.avgCPSF[1],avgvs.LNCPSF[1])

avg_table[,9]=c(apply(vs.bagshf.bic,1,sd)[1],apply(vs.bagshf.shd,1,sd)[1],apply(vs.bagshf.shdadj,1,sd)[1],apply(vs.bagshf.avgCPSF,1,sd)[1],apply(vs.bagshf.LNCPSF,1,sd)[1])

ll=nrow(avg_table)
for(i in 1:ll){
  tmp=paste(i,"&",avg_table[i,1],"&",avg_table[i,2],"(",round(avg_table[i,3],2),")","&",avg_table[i,4],"(",round(avg_table[i,5],2),")","&",avg_table[i,6],"(",round(avg_table[i,7],2),")","&",avg_table[i,8],"(",round(avg_table[i,9],2),")","\\",sep='')
  print(tmp)
}

