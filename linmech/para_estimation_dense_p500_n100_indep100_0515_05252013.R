setwd("~/repos/spacemap/sim/utils/linmech/")
source("summary_functions.R")
load("dense_p500_adj_matrix.Rdata")  ## true adj 
true.v=vstructures(true.dir)
#####################
nrep=100
n=100
p=nrow(true.ske)


######## get true covariance : gene.pancr.label, gene.pancr.label.ber
source("network_gene_function_10172012.R")
data_generation<-function(n,p,pancr.adj){

  Y.n=NULL
  panc.simu.Data=gene.pancr.label(n, pancr.adj,SN=runif(nrow(pancr.adj), min=0.5, max=1.5))

  Y=panc.simu.Data$data
  Y.gm=apply(Y, 2, mean)
  Y.gsd=apply(Y, 2, sd)
  Y.n=(Y-matrix(Y.gm, n, p, byrow=T))/matrix(Y.gsd, n, p, byrow=T)
  return(list(Y.n=Y.n,beta.true=t(panc.simu.Data$beta),sd.true=panc.simu.Data$esd))
}

set.seed(2948570)
seeds=sample(2:10295800,nrep)
for(rep in 1:nrep){
  print(rep)
  set.seed(seeds[rep])
  para=data_generation(n,p,true.dir)
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  save(para,file=file.name)
}


############ parameter estimation 
source("linear_mechanism_mu_sig_beta_estimation.R")

nrep=100
n=100
p=nrow(true.dir)
mu.true=matrix(0,nrep,p)
sigma.true=array(0,c(p,p,nrep))
con.true=array(0,c(p,p,nrep))

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')  ## load data 
  load(file.name)

  mu.tilta=rep(0,p)
  beta.true=para$beta.true
  sigma.tilta=(para$sd.true)^2
  
  node.order.true=cal_order(true.dir)
  temp=para_est2(mu.tilta, sigma.tilta, beta.true, true.dir, node.order.true)
  mu.true[rep,]=temp$mu.est
  sigma.true[,,rep]=temp$sigma.est
  con.true[,,rep]=solve(sigma.true[,,rep])
}


############################
#####summary and loss 
KL.loss.sample=numeric(nrep)
entropy.loss.sample=numeric(nrep)
con.l2.sample=numeric(nrep)
sigma.l2.sample=numeric(nrep)
for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)
  sample.cov=cov(para$Y.n)

  sigma.l2.sample[rep]=sqrt(sum((sample.cov-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
  entropy.loss.sample[rep]=sum(diag(sample.cov%*%con.true[,,rep]))-log(det(sample.cov%*%con.true[,,rep]))-p
}
mean(sigma.l2.sample)
mean(entropy.loss.sample)

sd(sigma.l2.sample)
sd(entropy.loss.sample)


###true adj_matrix
mu.est.true=matrix(0,nrep,p)
sigma.est.true=array(0,c(p,p,nrep))
con.est.true=array(0,c(p,p,nrep))
sigma.l2.true=numeric(nrep)
con.l2.true=numeric(nrep)
KL.loss.true=numeric(nrep)
entropy.loss.true=numeric(nrep)

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)

  Y=para$Y.n
  node.order1=cal_order(true.dir)
  temp1=para_linear(Y,true.dir)
  res1=para_est(temp1$mu.tilta, temp1$sigma.est, temp1$beta.est, true.dir, node.order1)
  mu.est.true[rep,]=res1$mu.est
  sigma.est.true[,,rep]=res1$sigma.est
  con.est.true[,,rep]=solve(sigma.est.true[,,rep])
    sigma.l2.true[rep]=sqrt(sum((sigma.est.true[,,rep]-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
  con.l2.true[rep]=sqrt(sum((con.est.true[,,rep]-con.true[,,rep])^2))/sqrt(sum(con.true[,,rep]^2))
  KL.loss.true[rep]=sum(diag(con.est.true[,,rep]%*%sigma.true[,,rep]))-log(det(con.est.true[,,rep]%*%sigma.true[,,rep]))-p
  entropy.loss.true[rep]=sum(diag(sigma.est.true[,,rep]%*%con.true[,,rep]))-log(det(sigma.est.true[,,rep]%*%con.true[,,rep]))-p

}
mean(sigma.l2.true)
mean(con.l2.true)
mean(KL.loss.true)
mean(entropy.loss.true)

sd(sigma.l2.true)
sd(con.l2.true)
sd(KL.loss.true)
sd(entropy.loss.true)

###BIC
mu.est.bic=matrix(0,nrep,p)
sigma.est.bic=array(0,c(p,p,nrep))
con.est.bic=array(0,c(p,p,nrep))

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)

  file.name=paste("bagshf_dense_p500_n_",n,"_rep",rep,".Rdata",sep='')
  load(file.name)

  Y=para$Y.n
  adj.matrix1=score.bagshf$adj.matrix
  node.order1=cal_order(adj.matrix1)
  temp1=para_linear(Y,adj.matrix1)
  res1=para_est(temp1$mu.tilta, temp1$sigma.est, temp1$beta.est, adj.matrix1, node.order1)
  mu.est.bic[rep,]=res1$mu.est
  sigma.est.bic[,,rep]=res1$sigma.est
  con.est.bic[,,rep]=solve(sigma.est.bic[,,rep])
  
}

mu.est.shd=matrix(0,nrep,p)
sigma.est.shd=array(0,c(p,p,nrep))
con.est.shd=array(0,c(p,p,nrep))

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)

  file.name=paste("bagshf_dense_p500_n_",n,"_rep",rep,".Rdata",sep='')
  load(file.name)

  Y=para$Y.n
  adj.matrix2=shd.bagshf$adj.matrix
  node.order2=cal_order(adj.matrix2)
  temp1=para_linear(Y,adj.matrix2)
  res1=para_est(temp1$mu.tilta, temp1$sigma.est, temp1$beta.est, adj.matrix2, node.order2)
  mu.est.shd[rep,]=res1$mu.est
  sigma.est.shd[,,rep]=res1$sigma.est
  con.est.shd[,,rep]=solve(sigma.est.shd[,,rep])
  
}

mu.est.adjshd=matrix(0,nrep,p)
sigma.est.adjshd=array(0,c(p,p,nrep))
con.est.adjshd=array(0,c(p,p,nrep))

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)

  file.name=paste("bagshf_dense_p500_n_",n,"_rep",rep,".Rdata",sep='')     
  load(file.name)

  Y=para$Y.n
  adj.matrix3=adjshd.bagshf$adj.matrix
  node.order3=cal_order(adj.matrix3)
  temp1=para_linear(Y,adj.matrix3)
  res1=para_est(temp1$mu.tilta, temp1$sigma.est, temp1$beta.est, adj.matrix3, node.order3)
  mu.est.adjshd[rep,]=res1$mu.est
  sigma.est.adjshd[,,rep]=res1$sigma.est
  con.est.adjshd[,,rep]=solve(sigma.est.adjshd[,,rep])
  
}

mu.est.cpsf=matrix(0,nrep,p)
sigma.est.cpsf=array(0,c(p,p,nrep))
con.est.cpsf=array(0,c(p,p,nrep))

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)

  file.name=paste("bagshf_dense_p500_n_",n,"_rep",rep,".Rdata",sep='')    
  load(file.name)

  Y=para$Y.n
  adj.matrix4=avgCPSF.bagshf$adj.matrix
  node.order4=cal_order(adj.matrix4)
  temp1=para_linear(Y,adj.matrix4)
  res1=para_est(temp1$mu.tilta, temp1$sigma.est, temp1$beta.est, adj.matrix4, node.order4)
  mu.est.cpsf[rep,]=res1$mu.est
  sigma.est.cpsf[,,rep]=res1$sigma.est
  con.est.cpsf[,,rep]=solve(sigma.est.cpsf[,,rep])
  
}

mu.est.lncpsf=matrix(0,nrep,p)
sigma.est.lncpsf=array(0,c(p,p,nrep))
con.est.lncpsf=array(0,c(p,p,nrep))

for(rep in 1:nrep){
  file.name=paste("dense_p500_n",n,"_rep",rep,"_indep100_para.Rdata",sep='')
  load(file.name)

  file.name=paste("bagshf_dense_p500_n_",n,"_rep",rep,".Rdata",sep='')   
  load(file.name)

  Y=para$Y.n
  adj.matrix5=LNCPSF.bagshf$adj.matrix
  node.order5=cal_order(adj.matrix5)
  temp1=para_linear(Y,adj.matrix5)
  res1=para_est(temp1$mu.tilta, temp1$sigma.est, temp1$beta.est, adj.matrix5, node.order5)
  mu.est.lncpsf[rep,]=res1$mu.est
  sigma.est.lncpsf[,,rep]=res1$sigma.est
  con.est.lncpsf[,,rep]=solve(sigma.est.lncpsf[,,rep])
  
}

save.image("para_estimation_p500_n100_indep100_0515_aggregated.Rdata")


##mu part
mu.l2.bic=mean(sqrt(rowSums(mu.est.bic^2))/p)  ##2.418124e-18
mu.l2.shd=mean(sqrt(rowSums(mu.est.shd^2))/p)  ##2.21983e-18
mu.l2.adjshd=mean(sqrt(rowSums(mu.est.adjshd^2))/p)  ##2.303366e-18
mu.l2.cpsf=mean(sqrt(rowSums(mu.est.cpsf^2))/p)
mu.l2.lncpsf=mean(sqrt(rowSums(mu.est.lncpsf^2))/p)

##sigma part
sigma.l2.bic=numeric(nrep)
sigma.l2.shd=numeric(nrep)
sigma.l2.adjshd=numeric(nrep)
sigma.l2.cpsf=numeric(nrep)
sigma.l2.lncpsf=numeric(nrep)

for(rep in 1:nrep){
  sigma.l2.bic[rep]=sqrt(sum((sigma.est.bic[,,rep]-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
  sigma.l2.shd[rep]=sqrt(sum((sigma.est.shd[,,rep]-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
  sigma.l2.adjshd[rep]=sqrt(sum((sigma.est.adjshd[,,rep]-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
  sigma.l2.cpsf[rep]=sqrt(sum((sigma.est.cpsf[,,rep]-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
  sigma.l2.lncpsf[rep]=sqrt(sum((sigma.est.lncpsf[,,rep]-sigma.true[,,rep])^2))/sqrt(sum(sigma.true[,,rep]^2))
}

mean(sigma.l2.bic)
mean(sigma.l2.shd)
mean(sigma.l2.adjshd)
mean(sigma.l2.cpsf)
mean(sigma.l2.lncpsf)

sd(sigma.l2.bic)
sd(sigma.l2.shd)
sd(sigma.l2.adjshd)
sd(sigma.l2.cpsf)
sd(sigma.l2.lncpsf)

##con
con.l2.bic=numeric(nrep)
con.l2.shd=numeric(nrep)
con.l2.adjshd=numeric(nrep)
con.l2.cpsf=numeric(nrep)
con.l2.lncpsf=numeric(nrep)

for(rep in 1:nrep){
  con.l2.bic[rep]=sqrt(sum((con.est.bic[,,rep]-con.true[,,rep])^2))/sqrt(sum(con.true[,,rep]^2))
  con.l2.shd[rep]=sqrt(sum((con.est.shd[,,rep]-con.true[,,rep])^2))/sqrt(sum(con.true[,,rep]^2))
  con.l2.adjshd[rep]=sqrt(sum((con.est.adjshd[,,rep]-con.true[,,rep])^2))/sqrt(sum(con.true[,,rep]^2))
  con.l2.cpsf[rep]=sqrt(sum((con.est.cpsf[,,rep]-con.true[,,rep])^2))/sqrt(sum(con.true[,,rep]^2))
  con.l2.lncpsf[rep]=sqrt(sum((con.est.lncpsf[,,rep]-con.true[,,rep])^2))/sqrt(sum(con.true[,,rep]^2))
}

mean(con.l2.bic)
mean(con.l2.shd)
mean(con.l2.adjshd)
mean(con.l2.cpsf)
mean(con.l2.lncpsf)

sd(con.l2.bic)
sd(con.l2.shd)
sd(con.l2.adjshd)
sd(con.l2.cpsf)
sd(con.l2.lncpsf)

##entropy
KL.loss.bic=numeric(nrep)
KL.loss.shd=numeric(nrep)
KL.loss.adjshd=numeric(nrep)
KL.loss.cpsf=numeric(nrep)
KL.loss.lncpsf=numeric(nrep)

for(rep in 1:nrep){
  KL.loss.bic[rep]=sum(diag(con.est.bic[,,rep]%*%sigma.true[,,rep]))-log(det(con.est.bic[,,rep]%*%sigma.true[,,rep]))-p
  KL.loss.shd[rep]=sum(diag(con.est.shd[,,rep]%*%sigma.true[,,rep]))-log(det(con.est.shd[,,rep]%*%sigma.true[,,rep]))-p
  KL.loss.adjshd[rep]=sum(diag(con.est.adjshd[,,rep]%*%sigma.true[,,rep]))-log(det(con.est.adjshd[,,rep]%*%sigma.true[,,rep]))-p
  KL.loss.cpsf[rep]=sum(diag(con.est.cpsf[,,rep]%*%sigma.true[,,rep]))-log(det(con.est.cpsf[,,rep]%*%sigma.true[,,rep]))-p
  KL.loss.lncpsf[rep]=sum(diag(con.est.lncpsf[,,rep]%*%sigma.true[,,rep]))-log(det(con.est.lncpsf[,,rep]%*%sigma.true[,,rep]))-p
}
mean(KL.loss.bic) ##229.05
mean(KL.loss.shd) ##196.72
mean(KL.loss.adjshd) ##182.19
mean(KL.loss.cpsf) ##193.53
mean(KL.loss.lncpsf) ##189.96

sd(KL.loss.bic) ##229.05
sd(KL.loss.shd) ##196.72
sd(KL.loss.adjshd) ##182.19
sd(KL.loss.cpsf) ##193.53
sd(KL.loss.lncpsf) ##189.96

####
entropy.loss.bic=numeric(nrep)
entropy.loss.shd=numeric(nrep)
entropy.loss.adjshd=numeric(nrep)
entropy.loss.cpsf=numeric(nrep)
entropy.loss.lncpsf=numeric(nrep)

for(rep in 1:nrep){
  entropy.loss.bic[rep]=sum(diag(sigma.est.bic[,,rep]%*%con.true[,,rep]))-log(det(sigma.est.bic[,,rep]%*%con.true[,,rep]))-p
  entropy.loss.shd[rep]=sum(diag(sigma.est.shd[,,rep]%*%con.true[,,rep]))-log(det(sigma.est.shd[,,rep]%*%con.true[,,rep]))-p
  entropy.loss.adjshd[rep]=sum(diag(sigma.est.adjshd[,,rep]%*%con.true[,,rep]))-log(det(sigma.est.adjshd[,,rep]%*%con.true[,,rep]))-p
  entropy.loss.cpsf[rep]=sum(diag(sigma.est.cpsf[,,rep]%*%con.true[,,rep]))-log(det(sigma.est.cpsf[,,rep]%*%con.true[,,rep]))-p
  entropy.loss.lncpsf[rep]=sum(diag(sigma.est.lncpsf[,,rep]%*%con.true[,,rep]))-log(det(sigma.est.lncpsf[,,rep]%*%con.true[,,rep]))-p
}
mean(entropy.loss.bic)  ##614.71
mean(entropy.loss.shd)  ##681.16
mean(entropy.loss.adjshd)  ##596.59 
mean(entropy.loss.cpsf) ##671.31
mean(entropy.loss.lncpsf)  ##637.52

sd(entropy.loss.bic)  ##614.71
sd(entropy.loss.shd)  ##681.16
sd(entropy.loss.adjshd)  ##596.59 
sd(entropy.loss.cpsf) ##671.31
sd(entropy.loss.lncpsf)  ##637.52
