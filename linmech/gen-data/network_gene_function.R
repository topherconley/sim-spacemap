######################
#####update on 09-18-2011 by jie; two functions: 1. gene.pancr and 2. gene.pancr.label
### updated by Ru on 10/17/2012: now allows both positive and negative beta coefficients by the "flip=TRUE" option
##add bernoulli 
######################


######################
gene.pancr<-function(N, pancr.adj, bmin=0.3, bmax=0.5, SN=runif(nrow(pancr.adj), min=0.5, max=1.5),flip=FALSE)
{
### originally written by Pei. need to order the nodes such that parents nodes always have smaller indices than children nodes 
#### also the DAG can only has one root which is node 1
### pancr.adj[i,j] indicates the existence of edge i--->j
### generate data based on a Gaussian linear mechnism and a DAG 
### para: N -- sample size, pancr.adj: adjacency matrix, bmin -- lower bound for the coefficients, bmax -- upper bound for the coefficients 
### SN: a p-vector , signal to noise ratio, used to determine the variance of the residuals 
### return: N by p data matrix 

  size=nrow(pancr.adj)   ##number of nodes 

  beta.simu=matrix(0, nrow=size, ncol=size)
  beta.simu[pancr.adj==1]=runif(sum(pancr.adj), min=bmin, max=bmax)  ##coefficients for linear mechanism 

  if(flip==TRUE){
    sign.temp=numeric(sum(pancr.adj==1))
    sign.temp=sample(c(-1,1),length(sign.temp),replace=TRUE)
    beta.simu[pancr.adj==1]=beta.simu[pancr.adj==1]*sign.temp
  }
  result=matrix(0, N, size)   ##data matrix 
  result[,1]=rnorm(N, mean=0, sd=1)  ##node one is the root, generate the root node 
  
  esd.v=rep(0, size)    ##sd for the residuals in the linear mechanism 
  esd.v[1]=1
  
  
  for(i in 2:size)
   {
      temp=beta.simu[,i]    ## 1:(i-1): parents coefficients for the ith node  
      if(i==2)
      {
        yhat=result[,1:(i-1)]*temp[1:(i-1)]
      } else {
         yhat=result[,1:(i-1)]%*%as.matrix(temp[1:(i-1)])
       }  
      ysd=sd(yhat)
      esd.v[i]=ysd/SN[i]    ## sd for the residual in the ith linear mechanism, set such that the ith SNR is SN[i]
      error=rnorm(N, mean=0, sd=esd.v[i])
      result[,i]=yhat+error
   }
   return(list(data=result, beta=beta.simu, esd=esd.v))
}


#######################
gene.pancr.label<-function(N, pancr.adj, bmin=0.3, bmax=0.5, SN=runif(nrow(pancr.adj), min=0.5, max=1.5),flip=FALSE)
{
### originally written by Ru on 09-16-2011: to generate a more general Gaussian linear mechanism
### applicable to any DAG 
### pancr.adj[i,j] indicates the existence of edge i--->j
### para: N -- sample size, pancr.adj: adjacency matrix, bmin -- lower bound for the coefficients, bmax -- upper bound for the coefficients 
### SN: a p-vector , signal to noise ratio, used to determine the variance of the residuals 
### return: N by p data matrix 

  size=nrow(pancr.adj)

  beta.simu=matrix(0, nrow=size, ncol=size)
  beta.simu[pancr.adj==1]=runif(sum(pancr.adj), min=bmin, max=bmax)
  if(flip==TRUE){
    sign.temp=numeric(sum(pancr.adj==1))
    sign.temp=sample(c(-1,1),length(sign.temp),replace=TRUE)
    beta.simu[pancr.adj==1]=beta.simu[pancr.adj==1]*sign.temp
  }

  result=matrix(0, N, size)
  
  
  tt=which(colSums(pancr.adj)==0)  ### find out nodes without parents 
  for(j in tt){
    result[,j]=rnorm(N, mean=0, sd=1)  ### generate them from N(0,1)
  }

  esd.v=rep(0, size)
  esd.v[tt]=1

  tt2=(1:size)[-tt]   ##nodes with at least one parents 
  current=tt2

  while(length(current)>0){

      i=current[1]                      ##look at node i  
      temp=beta.simu[,i]
      temp2=which(temp!=0)     ##parents for node i
      check=numeric(length(temp2))

      for(m in 1:length(temp2)){
        check[m]=all(result[,temp2[m]]==0)   ##check whether node i's parents have been generated: check =0 yes, check =1, not yet 
      }

      if(any(check!=0)){  ##if some parents of node i have not been generated 
        current=current[-1]
        current=c(current,i)
      }

      if(all(check==0)){  ## if all parenets of node i have been generated 
      
         yhat=as.matrix(result[,temp2],nrow=N)%*%as.matrix(temp[temp2])
         ysd=sd(as.vector(yhat))
         esd.v[i]=ysd/SN[i]
         error=rnorm(N, mean=0, sd=esd.v[i])
         result[,i]=yhat+error
         current=current[-1]
      } 

   }##end while 

   return(list(data=result, beta=beta.simu, esd=esd.v))
}

#######################
gene.pancr.label.ber<-function(N, pancr.adj, index.x, p.x, bmin=0.3, bmax=0.5, SN=runif(nrow(pancr.adj), min=0.5, max=1.5),flip=FALSE)
{
### originally written by Ru on 09-16-2011: to generate a more general Gaussian linear mechanism
### applicable to any DAG 
### index.x: indices for those parents which will be Ber
### p.x: same length as index.x which specifies Ber prob. 
### pancr.adj[i,j] indicates the existence of edge i--->j
### para: N -- sample size, pancr.adj: adjacency matrix, bmin -- lower bound for the coefficients, bmax -- upper bound for the coefficients 
### SN: a p-vector , signal to noise ratio, used to determine the variance of the residuals 
### return: N by p data matrix 

  size=nrow(pancr.adj)

  beta.simu=matrix(0, nrow=size, ncol=size)
  beta.simu[pancr.adj==1]=runif(sum(pancr.adj), min=bmin, max=bmax)

  if(flip==TRUE){
    sign.temp=numeric(sum(pancr.adj==1))
    sign.temp=sample(c(-1,1),length(sign.temp),replace=TRUE)
    beta.simu[pancr.adj==1]=beta.simu[pancr.adj==1]*sign.temp
  }

  result=matrix(0, N, size)
  
  if(!all(colSums(pancr.adj)[index.x]==0)){
     print("wrong index of bernoulli variable, generate Gaussian variable instead")
  }

  index.npar=which(colSums(pancr.adj)==0)  ### find out nodes without parents 
  index.par=(1:size)[-index.npar]   ##nodes with at least one parents 
  esd.v=rep(0, size)
  esd.v[index.npar]=1

  ###generate gaussian
  temp.gindex=(!(index.npar%in%index.x))
  if(sum(temp.gindex)>0){
    gaussian.npar=index.npar[temp.gindex]
    for(j in gaussian.npar){
      result[,j]=rnorm(N, mean=0, sd=1)  ### generate them from N(0,1)
    }
  }

  ###generate bernoulli
  for(j in 1:length(index.x)){
    result[,index.x[j]]=rbinom(N,1,prob=p.x[j])
  }

  current=index.par

  while(length(current)>0){

      i=current[1]                      ##look at node i  
      temp=beta.simu[,i]
      temp2=which(temp!=0)     ##parents for node i
      check=numeric(length(temp2))

      for(m in 1:length(temp2)){
        check[m]=all(result[,temp2[m]]==0)   ##check whether node i's parents have been generated: check =0 yes, check =1, not yet 
      }

      if(any(check!=0)){  ##if some parents of node i have not been generated 
        current=current[-1]
        current=c(current,i)
      }

      if(all(check==0)){  ## if all parenets of node i have been generated 
      
         yhat=as.matrix(result[,temp2],nrow=N)%*%as.matrix(temp[temp2])
         ysd=sd(as.vector(yhat))
         esd.v[i]=ysd/SN[i]
         error=rnorm(N, mean=0, sd=esd.v[i])
         result[,i]=yhat+error
         current=current[-1]
      } 

   }##end while 

   return(list(data=result, beta=beta.simu, esd=esd.v))
}

