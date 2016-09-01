###############updated by jie on 11-03-2011
############### functions to generate DAG's from linear mechanism 


############### part I: simulate DAG topology
#### 
DAG_net<-function(p,Pare.prob){
##para: p -- number of nodes; 
##Pare.prob: p by p matrix: indicating the chance of a lower numbered node being parent of a given higher numbered node
## 
##result: p by p adjacent matrix. Note that for a DAG, only a lower numbered node can be a parent of
##        a higher numbered node; thus only the lower half of the matrix holds information. 
##        i.e., result[i,j]=0 for j>=i

result<-matrix(0,p,p)
  for(i in 2:p){
   prob.c<-as.vector(Pare.prob[i,1:(i-1)])
   prob.c<-cbind(1-prob.c,prob.c)
   pare.c<-apply(prob.c, 1, sample, x=0:1, replace=FALSE, size=1)
   result[i,1:(i-1)]<-pare.c
  }
  
 rownames(result)<-paste("node",1:p)
 colnames(result)<-paste("node",1:p)

 return(result)

}

#Author: Chris Conley
dag_xy_xhub <- function(px, py, nxhub = 20, min_hub_size = 5, mean_hub_size = 12) { 
  
  #construct x hubs with y targets ( x->y )
  index_xhubs <- sample(px, nxhub)
  xy_pool <- expand.grid(x = index_xhubs, y = seq_len(py))
  hub_size_pool <- pmax(rpois(1000, mean_hub_size),min_hub_size)
  #hist(hub_size_pool)
  y_pool <- seq_len(py)
  samp_y_edges <- function(x) { 
    
  }
  
  xy_edges <- list()
  for(i in seq_along(index_xhubs)) { 
    y_hits <- sample(x = y_pool, size = sample(hub_size_pool,1))
    xy_edges[[i]] <- y_hits
    y_still_in_pool <- as.logical(rbinom(n = length(y_hits), size = 1, prob = 0.05))
    y_pool <- setdiff(y_pool, y_hits[y_still_in_pool])
  }
  
  x2y_adj <- matrix(0,px,py)
  for (i in seq_along(index_xhubs)) x2y_adj[index_xhubs[i],xy_edges[[i]]] <- 1
  x2y_adj
}


###############################
#####  generate layer 1:  X->Y network: only relevant X and Y-nodes in Markov blanket of X  are remained in the adjacent matrix
DAG_XY<-function(px=800,py=800,rel=NULL,pc=0.1,pt=0.1/py,pco=0.3/py, xhub = TRUE){
##pare: px -- number of X-nodes; py -- number of Y-nodes; 
##      rel -- 0-1 matrix (px by py) indicating cis, trans between X and Y nodes
##      pc -- chance of cis-regulation between X-Y cis pairs; pt -- chance of trans-regulation
##      pco -- chance of Y-nodes co-regulate another Y-node with a X -- node
##result: weighted adjacent matrix: 1 -- cis-regulate, 2 -- trans-regulate, 3 -- co-regulate; 
## only relevant X and Y-nodes in Markov blanket of X  are remained in the adjacent matrix

  if(is.null(rel)){
   rel<-matrix(0,px,py)
  }


Pare.prob<-matrix(0,px+py,px+py)

##X->Y prob.
if (xhub) { 
  temp <- t(rel)
} else { 
  temp<-pc*t(rel)+pt*(1-t(rel))              
}


Pare.prob[(px+1):(px+py),1:px]<-temp

temp<-matrix(pco,py,py)                ##Y-->Y prob.
Pare.prob[(px+1):(px+py),(px+1):(px+py)]<-temp

adj1<-DAG_net(px+py,Pare.prob)

## get rid of irrelevant Y and X: i.e.,those do not have any children and parents 
temp<-adj1[,1:px]
index.x<- (1:px)[apply(temp,2,sum)>0]        ##find out X with at least one kids

temp<-adj1[-(1:px),1:px]
index.yx<-(1:py)[apply(temp,1,sum)>0]       ##find out Y with at least one parent in X 

temp<-adj1[index.yx+px, -(1:px)]
index.co<-(1:py)[apply(temp,2,sum)>0]        ##find out co-parents of X in Y  

index.y<-sort(union(index.yx,index.co))      ## Ys in the first layer    
index.1<-c(index.x,index.y+px)               ##X and first layer Y

adj.1<-adj1[index.1,index.1]                 ##first layer subnetwork

##find out trans-edges and weighted as 2

px.new<-length(index.x)
adj.1[-(1:px.new),1:px.new]<-adj.1[-(1:px.new),1:px.new]+(1-t(rel[index.x,index.y]))*adj.1[-(1:px.new),1:px.new]


##find out the indices for index.co nodes, and then use different weights for Y--> Y edges: coded as 3
index.co.new<-(1:length(index.y))[is.element(index.y,index.co)]
adj.1[,index.co.new+px.new]<-adj.1[,index.co.new+px.new]*3


##weighted (1,2,3) adjacent matrix for further use 
adj.xy.we<-adj.1                         

return(list("adjacent"=adj.xy.we,"X_index"=index.x, "Y_index"=index.y,"co_Y_index"=index.co))


}


#####  plot layer 1 X->Y netowrk
Plot_XY<-function(adj.xy.we, index.x,  hub.size=2, lay=layout.fruchterman.reingold){
## para: adj.xy.we --weighted ad network from DAG_XY
##       index.x: indices for X-nodes: 1:px.new; index.co: indices for co-parent nodes
##       hub.size: definition of hub-nodes
##       lay: layout 

library(igraph)
temp=graph.adjacency(adjmatrix=t(adj.xy.we), mode="directed",weighted=TRUE)
temp$layout<-lay

V(temp)$color<- "green" #3                               ##X--red (2),Y-green (3), Y-co: blue (4)
V(temp)$color[1:length(index.x)]<- "red" #2
index.co<-(1:nrow(adj.xy.we))[apply(adj.xy.we==3,2,sum)>0]        ##find out co-parents of X in Y  
V(temp)$color[index.co]<-  "blue" #4


V(temp)$size<-numeric(nrow(adj.xy.we))+3         ##hub nodes larger in size
hub.index<-apply(adj.xy.we>0,2,sum)>=hub.size            ###indices for hub nodes 
V(temp)$size[hub.index]<-5

V(temp)$label<-""

E(temp)$color= "red" #2                            ##cis edge: red (2); trans edge: black (1); Y edge: green (3)
E(temp)$color[E(temp)$weight==2]= "black" #1
E(temp)$color[E(temp)$weight==3]= "green" #4

##
plot.igraph(temp,edge.width=0.5, edge.arrow.size=0.2)

x.no<-length(index.x)
y.no<-nrow(adj.xy.we)-x.no
no.cis<-sum(adj.xy.we==1)     ##number of cis-edges      
no.tr<-sum(adj.xy.we==2)      ##number of trans-edges     
no.Y<-sum(adj.xy.we==3)       ##number of Y-edges         
edge.no<-no.cis+no.tr+no.Y

title(paste("X (red):", x.no, ", Y (gre, blue):", y.no, "; cis (red):",no.cis, ", trans (black):", no.tr, ", Y-co  (blue):", no.Y,sep=" "))
}


#########################################################
###### generate layer 2 network: Y->Z
DAG_YZ<-function(pz,result.xy,hub.no=c(5,5,5,10),hub.size.pool, p.o){

##para: pz -- number of Z nodes : i.e., nodes only connected with Y, but not directly connected with X
##result.xy -- X-Y network from DAG_XY
##      hub.no -- number of each types of hubs: cis, trans, co-parent, and Z-hub
##      hub.size.pool: an integer vector where we draw the size of each hub from 
##      p.o  -- chance of non-hub edges between Z-nodes

adj.xy.we<-result.xy$adjacent
index.x<-result.xy$X_index
index.co<-result.xy$co_parent_index-length(index.x)  ##co-parenet-Y-node index

px.new<-length(index.x)
py.new<-nrow(adj.xy.we)-px.new

pyz=pz+py.new 

temp<-adj.xy.we[-(1:px.new),1:px.new]   ##X-->Y
index.cis.new<-(1:py.new)[apply(temp==1,1,sum)>0]   ##cis-Y-node index
index.trans.new<-(1:py.new)[apply(temp==2,1,sum)>0]  ##trans-Y-node index 


##specify size of hubs: drawing from hub.size.pool
hub.y.size<-sample(hub.size.pool,size=hub.no[1]+hub.no[2]+hub.no[3], replace=TRUE)
hub.z.size<-sample(hub.size.pool,size=hub.no[4],replace=TRUE)

##randomly select hub-Y based on the hub numbers 
#set.seed(2)

hub.no[1]=min(hub.no[1], length(index.cis.new))
hub.no[2]=min(hub.no[2], length(index.trans.new))
hub.no[3]=min(hub.no[3], length(index.co))

index.hub.cis<-sample(index.cis.new, hub.no[1],replace=FALSE)
index.hub.tr<-sample(index.trans.new, hub.no[2],replace=FALSE)
index.hub.co<-sample(index.co, hub.no[3],replace=FALSE)

index.hub.y<-sort(union(union(index.hub.cis,index.hub.tr),index.hub.co)) ##Y-hub indices
hub.y.size<-hub.y.size[1:length(index.hub.y)]        ##deal with possibly overlapping hubs in Y


#index.hub.z<-sort(sample((1+py.new):pyz, hub.no[4],replace=FALSE))
hub.no[4]=min(hub.no[4], pz-1)
index.hub.z<-(1+py.new): (1+py.new+hub.no[4]-1)  ##Z hub index: set as the first hub.no[4 ] Z-nodes 

## selected hub-kids from Z nodes
Pare.prob.yz<-matrix(0,pyz,pyz)
 for(i in 1:length(index.hub.y)){  ##Y-hubs: kids
  index.y.c<-index.hub.y[i]
  size.c<-min(hub.y.size[i],pz)  ##number of kids can not exceed pz
 
  kid.c<-sample((1+py.new):pyz,size=size.c,replace=FALSE)
  Pare.prob.yz[kid.c,index.y.c]<-1
 }

for(i in 1:length(index.hub.z)){ ##Z-hubs: kids (note all kids need to have higher indices than parents)
  index.z.c<-index.hub.z[i]
  size.c<-min(hub.z.size[i], pyz-index.z.c)  ##number of kids can not exceed the number of Z-nodes with a higher indices 
   kid.c<-sample((1+index.z.c):pyz,size=size.c,replace=FALSE)
  Pare.prob.yz[kid.c,index.z.c]<-1
 }


##specify chances of other Y1->Z and Z1->Z connections, where Y1, Z1 are not hubs 
Pare.prob.yz[-(1:py.new), -c(index.hub.y,index.hub.z)]<-p.o

##
#set.seed(3)
adj2<-DAG_net(pyz,Pare.prob.yz)

##
adj2.we<-adj2        ##weight the adjacent matrix: Y->Z,coded as 4, Z-->Z coded as 5; 
adj2.we[-(1:py.new),1:py.new]<-4*adj2.we[-(1:py.new),1:py.new]
adj2.we[-(1:py.new),-(1:py.new)]<-5*adj2.we[-(1:py.new),-(1:py.new)]
adj2.we[1:py.new,1:py.new]<-adj.xy.we[-(1:px.new), -(1:px.new)] ##Y->Y network from layer 1: coded as 3

##
result<-list("adjacent"=adj2.we, "cis_hub"=index.hub.cis, "trans_hub"=index.hub.tr, "co_hub"=index.hub.co, "Z_hub"=index.hub.z)

return(result)
}


##### plot layer 2 Y->Z netowrk
Plot_YZ<-function(adj.yz.we, pz, hub.size=5, lay=layout.fruchterman.reingold){
## para: adj.yz.we --weighted ad network from DAG_YZ
##       pz -- number of Z-nodes
##       hub.size: definition of hub-nodes
##       lay: layout 

pyz<-nrow(adj.yz.we)
py<-pyz-pz

temp=graph.adjacency(adjmatrix=t(adj.yz.we), mode="directed",weighted=TRUE)
temp$layout<-lay


V(temp)$color<-4                               ##Z--blue (4),Y-green (3)
V(temp)$color[1:py]<-3

V(temp)$size<-numeric(pyz)+3         ##hub nodes larger in size
hub.index.YZ<-apply(adj.yz.we>0,2,sum)>=hub.size            ###index for hub Y,Z 
V(temp)$size[hub.index.YZ]<-5
                        
E(temp)$color[E(temp)$weight==3]=3              ##Y-Y edge: green; Y-Z edge: black; Z-Z edge: blue
E(temp)$color[E(temp)$weight==4]=1
E(temp)$color[E(temp)$weight==5]=4

V(temp)$label<-""

##
plot(temp,edge.width=0.5, edge.arrow.size=0.2)
yy.edge<-sum(adj.yz.we==3)            ##Y-Y edge
yz.edge<-sum(adj.yz.we==4)            ##Y-Z edge
zz.edge<-sum(adj.yz.we==5)            ##Z-Z edge
pz.conn<-pz-sum(apply(adj.yz.we[-(1:py),]>0,1,sum)==0&apply(adj.yz.we[,-(1:py)]>0,2,sum)==0)  ##Z-nodes connect to the network

title(paste("Y (gre):", py, ";Z-conn :", pz.conn, ";Y-Y (gre):", yy.edge, ";Y-Z (bla):", yz.edge, ";Z-Z (blue):",zz.edge))

}

############################################################
#####  generate whole network: X->Y->Z
DAG_XYZ<-function(px,py, pz, adj.xy.we,adj.yz.we){
##para: px -- number of X nodes; py -- number of Y nodes; pz --number of Z nodes 
##      adj.xy.we -- X-Y adjacent matrix from DAG_XY;  adj.yz.we -- Y-Z adjacent matrix from DAG_YZ; 
##
adj.we<-matrix(0,px+py+pz,px+py+pz)
adj.we[1:(px+py),1:(px+py)]<-adj.xy.we  ##X-->Y network
adj.we[-(1:px),-(1:px)]<-adj.yz.we   ##(Y,Z) network

##summary
index.no.all<-(1:(px+py+pz))[apply(adj.we>0,2,sum)==0&apply(adj.we>0,1,sum)==0]  ##nodes without parents and kids

index.cis.y<-(1:(px+py+pz))[apply(adj.we==1,1,sum)>0]
index.tr.y<-(1:(px+py+pz))[apply(adj.we==2,1,sum)>0]
index.co.y<-(1:(px+py+pz))[apply(adj.we==3,2,sum)>0]


index.conn<-setdiff(1:(px+py+pz),index.no.all)
index.conn.x<-index.conn[index.conn<=px]
index.conn.z<-index.conn[index.conn>(px+py)]
index.conn.y<-(px+1):(px+py)

##
result<-list("adjacent"=adj.we, "cis_Y"=index.cis.y, "trans_Y"=index.tr.y, "co_Y"=index.co.y, "X_conn"=index.conn.x, "Y_conn"=index.conn.y, "Z_conn"=index.conn.z)

return(result)
}

#### plot X->Y->Z network
Plot_XYZ<-function(result.xyz,hub.size=5, lay=layout.fruchterman.reingold){
##para: result.xyz.we: X-Y-Z network from DAG_XYZ
##      hub.size -- definition of hubs 
##      lay: layout

adj.xyz.we<-result.xyz$adjacent
index.conn.x<-result.xyz$X_conn
index.conn.y<-result.xyz$Y_conn
index.conn.z<-result.xyz$Z_conn

px.new<-length(index.conn.x)
py.new<-length(index.conn.y)
pz.new<-length(index.conn.z)
index.conn<-c(index.conn.x,index.conn.y,index.conn.z)


##plot: get rid of non-connection nodes
adj.we.new<-adj.xyz.we[index.conn,index.conn]
temp=graph.adjacency(adjmatrix=t(adj.we.new), mode="directed",weighted=TRUE)
temp$layout<-lay


V(temp)$color<-4                               ##X--red (2),Y-green (3),Z--blue (4)
V(temp)$color[1:px.new]<-2
V(temp)$color[(px.new+1):(px.new+py.new)]<-3

V(temp)$size<-numeric(nrow(adj.we.new))+3         ##hub nodes larger in size
hub.index.all<-apply(adj.we.new>0,2,sum)>=hub.size            ###indices for hub Y,Z 
V(temp)$size[hub.index.all]<-5


E(temp)$color[E(temp)$weight==1]=2              ##cis edge: red; trans edge: black            
E(temp)$color[E(temp)$weight==2]=1             
E(temp)$color[E(temp)$weight==3]=3             ##Y-Y edge: green; Y-Z edge: black; Z-Z edge: blue    
E(temp)$color[E(temp)$weight==4]=1
E(temp)$color[E(temp)$weight==5]=4

V(temp)$label<-""

##
plot(temp,edge.width=0.5, edge.arrow.size=0.2)
cis.edge<-sum(adj.we.new==1)
tr.edge<-sum(adj.we.new==2)
yy.edge<-sum(adj.we.new==3)
yz.edge<-sum(adj.we.new==4)
zz.edge<-sum(adj.we.new==5)

title(paste("X (red):",px.new,";Y (gre):", py.new, ";Z (blue):", pz.new, ";cis (red):", cis.edge, ";trans (bla):", tr.edge, "; yy (gre):", yy.edge, "; yz (bla):", yz.edge, "; zz (blu):", zz.edge))

}



################### convert a DAG to the undirected graph (Moral graph)
Moral<-function(adj.dir){
## para: adj.dir: adjacent matrix  of the DAG
## result:  ajacent matrix of the corresponding moral graph
## this is weighted: 1--original directed edge; 2 -- undirected edge caused by directed edge
## 3 -- undirected edge by marrying parents
  
  p<-nrow(adj.dir)
  result<-matrix(0,p,p)
  result<-(adj.dir>0)
  
  result<-result+t(2*result)
    
   for(i in 1:p){ 
    index.p<-(1:p)[adj.dir[i,]>0]
     if(length(index.p)>1){  ##marrying parents
       result[index.p,index.p]<-3    
      }    
   }

  diag(result)<-0
  rownames(result)<-paste("p",1:p)
  colnames(result)<-paste("c",1:p)
  
  return(result)
  
}


#################################### Get the coefficient matrix of
####################################
DAG_coef<-function(B.m,inter.v){
##B.m: p by p coefficient matrix which specifies the linear mechanism between nodes; when there is no edge from j--> i: B.m[i,j]=0
##inter.v: p by 1 intercept vector for the linear mechanisms
##result:p by (p+1) matrix which specifies each node as a linar combination of the residual vector
##the last colum corresponds to intercept of each linear equation
## return:  DAG nodes with respect to the residual vector, i.e.,  node_values=result[,1:p]%*%residual 

 p<-length(inter.v)  
 result<-matrix(0,p,p+1)
  
 result[1,1]<-1
 result[1,p+1]<-inter.v[1]
 
 for (i in 2:p){
 coef.c<-matrix(result[1:(i-1),], nrow=i-1, ncol=p+1)
 beta.c<-matrix(B.m[i,1:(i-1)], nrow=1,ncol=i-1)
 
 result[i,]<-beta.c%*%coef.c
 result[i,i]<-1
 result[i,p+1]<-result[i,p+1]+inter.v[i]
 
  }
  
 colnames(result)<-c(paste("resi",1:p),"intercept")
 rownames(result)<-c(paste("node",1:p))
 return(result)

}


################################ Generate reisdual variance based on signal to noise ratio (sn.pool)    
Resi_var<-function(adj.we, coef.r, cov.X, sn.pool,sig2=1,prop=FALSE){
##para: adj.we: weighted adjacent matrix 
##      coef.r: coefficient matrix from DAG_coef; 
##      cov.X: covariance matrix of X-nodes
##      sn.pool: the pool for signal noise ratio
##      sig2: varinace of Y/Z nodes which do not have parents
##      prop:  whether signal to noise ratio should be proportional to the number of parents 
##result:residual variance: for X-node, use diag(cov.X), for Y,Z-node without parents, use sig2  
## SNR := signal part variance/residual variance 

px<-nrow(cov.X)
pyz<-nrow(adj.we)-px
index.wp<-(1:(px+pyz))[apply(adj.we>0,1,sum)>0]    ##indices for variables with parents 
sn<-sample(sn.pool,size=length(index.wp))          ##get signal to noise ratio for nodes with parents

if(prop){
pare.no<-apply(adj.we[index.wp,]>0,1,sum) ##number of parents
sn[order(pare.no)]<-sort(sn)
}

##Find out residual variance to match the signal to noise ratio for nodes with parents
resi.var<-numeric(px+pyz)
resi.var[-index.wp]<-sig2          ##for Y-nodes without parents, set the residual (i.e, Y) variance as sig2
resi.var[1:px]<-diag(cov.X)        ## set variance for X nodes 

 for (i in 1:length(index.wp)){
   index.c<-index.wp[i]
   sn.c<-sn[i]
   coef.c<-matrix(coef.r[index.c,1:(index.c-1)])
   resi.cov<-matrix(0, index.c-1, index.c-1)
   diag(resi.cov)<-resi.var[1:(index.c-1)]
   resi.cov[1:px,1:px]<-cov.X
   signal.c<-t(coef.c)%*%resi.cov%*%coef.c
   resi.var[index.c]<-signal.c/sn.c
 }

return(resi.var)


}

################################ Get mean and covariance for the nodes
################################ Also get regression coefficients: Y ~ X+Y
################################# and residual variance in the regressions 
Dist_node<-function(coef.r, cov.X, mean.X, resi.var){
##para: coef.r: coefficient matrix from DAG_coef; 
##      cov.X: covariance matrix of X-nodes; mean.X: mean of X-nodes
##      resi.var: residual variance from Resi_var
##return: mean and covariance of nodes; regression coefficients, and regression residual variance

px<-length(mean.X)
py<-length(resi.var)-px

resi.cov<-matrix(0,px+py,px+py)
diag(resi.cov)<-resi.var
resi.cov[1:px,1:px]<-cov.X
cov.node<-coef.r[,1:(px+py)]%*%resi.cov%*%t(coef.r[,1:(px+py)])

resi.mean<-numeric(px+py)
resi.mean[1:px]<-mean.X
mean.node<- coef.r[,1:(px+py)]%*%matrix(resi.mean)+coef.r[,px+py+1]


##regression coefficient for moral graph
 beta.m<-matrix(0,py,px+py) 
 beta.o<-numeric(py)
 resi.reg.var<-numeric(py)
 cov.inv<-solve(cov.node)

 
 for(i in 1:py){
 print(i)
 y.index<-i+px
 x.index<-(1:(px+py))[-y.index]
 ##x.cov<-cov.node[x.index,x.index]
 xy.cov<-matrix(cov.node[y.index,x.index])
 
 x.cov.inv<-Inv_Sub(cov.inv,y.index)
 beta.c<-x.cov.inv%*%xy.cov            
 beta.m[i,-y.index]<-beta.c                
 beta.o[i]<-mean.node[y.index]-sum(mean.node[x.index]*beta.c)
 
 resi.reg.var[i]<-cov.node[y.index,y.index]-t(xy.cov)%*%beta.c
 }

## 
result<-list("cov_node"=cov.node, "mean_node"=mean.node, "reg_coef"=beta.m, "reg_inter"=beta.o, "reg_resi_var"=resi.reg.var)
return(result)

}







############################# auxiliary functions
### given the inverse of a p by p matrix, get the inverse of its p-1 by p-1 submatrix (deleting ith row and ith col)
### only applicable to symmetric matrices
Inv_Sub<-function(Sig.Inv,i){
##para: Sig.inv: inverse of the p by p matrix
##i: the index for the row and column deleted from the original matrix
##result: the inverse of the corresponding submatrix

p<-nrow(Sig.Inv)
temp<-Sig.Inv
temp[1:(p-1),]<-Sig.Inv[-i,]
temp[p,]<-Sig.Inv[i,]
temp11<-temp[,i]
temp[,1:(p-1)]<-temp[,-i]
temp[,p]<-temp11


temp1<-temp[1:(p-1),1:(p-1)]
d<-temp[p,p]
b<-temp[1:(p-1),p]

result<-temp1-1/d*b%*%t(b)

return(result)

}








########### alternative way to generate node data: this is not as efficient
########### since for each replicate, the whole process needs to be repeated
DAG_var<-function(index.root,value.obs, coef.m, inter.v, sd.v){
##para: index.root: indices for nodes without parents, i.e., root nodes
##value.obs: p by 1 observed value vector; Initially only exisit for the root nodes, 
## and for other nodes, the value needs to be updated
##coef.m: p by p coefficient matrix which specifies the linear mechanism; when there is no edge: coef.m[i,j]=0
##inter.v: p by 1 intercept vector for the linear mechanism
##sd.v: p by 1 residual sd 
##return: result=updated value.obs

result<-value.obs
p<-length(value.obs)

  for (i in 1:p){
     if(!is.element(i,index.root)){
          coef.c<-coef.m[i,]
          obs.c<-sum(result*coef.c)+inter.v[i]+sd.v[i]*rnorm(1)
          #obs.c<-sum(result*coef.c)+inter.v[i]+sd.v[i]*resi[i] 
          result[i]<-obs.c 
     }
   }

return(result)

}















