result_skeleton<-function(adj.m, true.ske){
## find total number of skeleton and number of corrected skeleton from the adj of an estimated DAG (compared with the true.ske)
##para: adj.m: adjacency matrix, true.ske -- true skeleton (note: symmetric)

diag(adj.m)=0
tt=adj.m+t(adj.m)
correct.c=sum((tt>0)&(true.ske>0))/2    
total.c=sum(tt>0)/2

return(c(total.c, correct.c))

}

vstructures<-function(adj.matrix){
### 06/27/2011
### function to find vstructures in an adjacency matrix 
##para: adj.matrix: adjacency matrix
##return: 3-column matrix: (par1, child , par2): note, par 1< par 2     p=ncol(adj.matrix)
     res=NULL
  for(i in 1:p){
    parent.node=which(adj.matrix[,i]!=0)
    n.par=length(parent.node)    if(n.par>1){
      for(j in 1:(n.par-1)){
        for(k in (j+1):n.par){
          if(abs(adj.matrix[parent.node[j],parent.node[k]])+abs(adj.matrix[parent.node[k],parent.node[j]])==0){
            res=rbind(res,c(parent.node[j],i,parent.node[k]))   
          }
        }
      }##for loop
    }##end if
  }## end i loop
  if(!is.null(res)){
    colnames(res)=c("par1","child","par2")
  }
  return(res)
}result_vstruct<-function(adj.m,true.v){
##para: adj.: adjacency matrix 
##true.v: 3 column matrix (par1, child, par2), note par1<par2
##result: total number of vstructure and number of correct vstructure. also returned are the estimated v-structure and the correct subset

  total.v=vstructures(adj.m)   ##get the v-structure of adj.m: 3 column matrix 

  if(is.null(total.v)){
     total.v=matrix(0,ncol=3)
     corr.v=matrix(0,ncol=3)
     colnames(total.v)=c("par1","child","par2")
     colnames(corr.v)=c("par1","child","par2")
     
     total.num=0
     corr.num=0
  }else{
     total.num=nrow(total.v)
     corr.v=compare.vstructures(total.v,true.v)
  
  if(is.null(corr.v)){
     corr.v=matrix(0,ncol=3)
     colnames(corr.v)=c("par1","child","par2")
     
     corr.num=0
   }else{
     corr.num=nrow(corr.v)
   }  
  }
   return(list(estimated.vstru=total.v,corr.vstru=corr.v,total.num=total.num,corr.num=corr.num))
}

compare.vstructures<-function(target.vstructures,true.vstructures){
##compare two v-structures: two 3-column matrix
##return: a 3-column matrix, the subset of true-vstructure in the target vstructure

   corr.v=NULL
   if(!is.null(target.vstructures)){
   target.vstructures=matrix(target.vstructures,ncol=3)

     target.l=nrow(target.vstructures)
     for(i in 1:target.l){
       res=apply(true.vstructures,1,function(x) all(x==target.vstructures[i,]))
       
       if(any(res==TRUE)){
         corr.v=rbind(corr.v,true.vstructures[which(res==TRUE),])
       }##end if 
     }##end for
   }##end if 
   
   return(corr.v)
}

moral_graph<-function(adj.matrix){
p=nrow(adj.matrix)
moral.adj.matrix=adj.matrix
n.parents=colSums(adj.matrix)
for(i in 1:p){
  if(n.parents[i]>=2){
    pars=which(adj.matrix[,i]>0)
    for(j in 1:(n.parents[i]-1)){
      for(k in (j+1):n.parents[i]){
          moral.adj.matrix[pars[j],pars[k]]=1
      }#end for
    }#enf for
  }#end if
}#end for
moral.adj.matrix=((moral.adj.matrix+t(moral.adj.matrix))>0)
return(moral.adj.matrix)
}

weight.avg<-function(res.list){
  nrep=length(res.list)
  n1_max=NULL
  n2_max=NULL
  for(rep in 1:nrep){
    tmp1=ncol(res.list[[rep]])
    n1_max=c(n1_max,res.list[[rep]][1,tmp1])
  }
  tmp2=max(n1_max)

  cves=matrix(0,2,tmp2)
  last.cves=0
  for(j in 1:tmp2){
    tmp3=numeric(2)
    tmp.count=0
    for(rep in 1:nrep){
      tmp4=which(res.list[[rep]][1,]==j)
      if(length(tmp4)>0){
        tmp3=tmp3+res.list[[rep]][,tmp4[1]]
        tmp.count=tmp.count+1
      }
    }##end for
    if(tmp.count>0){
      cves[,j]=tmp3/tmp.count
      last.cves=cves[,j]
    }
    else{
      cves[,j]=last.cves
    }

  }##end for

  return(cves)

}

move_adj<-function(p,movement){
   adj.result=NULL
   adj.matrix=matrix(0,p,p)
   step=nrow(movement)

   for(i in 1:step){
     if(movement[i,3]>0){
       adj.matrix[movement[i,1],movement[i,2]]=(2-movement[i,3])>0   ##updated by jie on 1/10/2012 for efficiency 
       adj.matrix[movement[i,2],movement[i,1]]=(2-movement[i,3])<0
     }
     adj.result[[i]]=adj.matrix
   }
   return(adj.result)
}

