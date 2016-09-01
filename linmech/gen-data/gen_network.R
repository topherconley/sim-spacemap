library(igraph)
#################### Part I: generate the first layer network: X->Y
source("~/repos/sim-spacemap/linmech/gen-data/gen_dag_functions.R")

########## set parameters 
px.ini<-250                     ## initial number of X nodes
py.ini<-300                     ##initial  number of Y nodes 

xhub <- TRUE

if (xhub) { 
  rel <- dag_xy_xhub(px = px.ini, py = py.ini, nxhub = 10, min_hub_size = 5, mean_hub_size = 12)
} else { 
  rel<-matrix(0,px.ini,py.ini)        ## px by py matrix indicating cis (1), trans (0) relationship
  diag(rel[,1:px.ini])<-1   ##set X, Y nodes with same indices as cis-pair
}

pc<-0.5                    ##chance of cis-regulation between cis-pairs
pt<-0.5/py.ini                  ##chance of trans-regulation between trans-pairs
pco<-0.6/py.ini                 ##chance of being co-parent



set.seed(100)
result.xy<-DAG_XY(px.ini,py.ini,rel,pc,pt,pco, xhub = xhub)
adj.xy.we<-result.xy$adjacent         ## cis:1; trans :2 ; co-expression (Y-Y): 3

index.x<-result.xy$X_index
px=length(index.x)              ##resulting number of X nodes
#hist(colSums(adj.xy.we[,index.x]))
#hist(rowSums(adj.xy.we[index.x,]))
index.y<-result.xy$Y_index
py=length(index.y)     ##resulting number of Y nodes: either children of X nodes, of co-parents with X nodes 

index.co<-result.xy$co_Y_index
dim(adj.xy.we)           
px           ##176  X-nodes 
py           ##203 Ynodes  
length(index.co)      ##among the Y nodes: 51 co-parenet nodes 

sum(adj.xy.we>0)

out.deg=apply(adj.xy.we, 1, sum)
in.deg=apply(adj.xy.we, 2, sum)

par(mfrow=c(2,1))
hist(out.deg)
hist(in.deg)

par(mfrow=c(1,1))


##plot X-->Y layer 
Plot_XY(adj.xy.we, index.x, hub.size=2, lay=layout.fruchterman.reingold)


##################### Part II:generate Y->Z network

set.seed(2)
pz=125
hub.no<-c(4,5,4,5)              ##number of hubs for Y-cis, Y-trans, Y-co, and Z
#hub.size.pool<-sample(5:10, size=100000, replace=TRUE)         ##create pool for hub sizes
hub.size.pool<-rpois(100000, 5)+5
summary(hub.size.pool)
hist(hub.size.pool)
 ##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 ##  5.00    8.00   10.00   10.01   11.00   21.00 

p.o<-0.8/(py+pz)                      ## chnaces for Z-Z connection

result.yz<-DAG_YZ(pz,result.xy,hub.no,hub.size.pool, p.o)
adj.yz.we<-result.yz$adjacent         ##Y->Y: 3; Y->Z: 4; Z->Z: 5 
sum(adj.yz.we==3)   ##Y-->Y
sum(adj.yz.we==4)  ##Y-->Z
sum(adj.yz.we==5)  ##Z-->Z

##plot Y-Z network
Plot_YZ(adj.yz.we, pz, hub.size=5, lay=layout.fruchterman.reingold)

####################### Part III: generate X->Y->Z network
result.xyz<-DAG_XYZ(px,py, pz, adj.xy.we,adj.yz.we)
adj.xyz.we<-result.xyz$adjacent
dim(adj.xyz.we)
all(adj.xyz.we[upper.tri(adj.xyz.we)]==0)  ##all edges are from lower index to higher index: true, i.e., no edge from j (col)-->i (row), if j>i


table(apply(adj.xyz.we>0,2,sum))          ##number of kids

table(apply(adj.xyz.we>0,1,sum))          ##number of parents

sum(adj.xyz.we>0)                         ##total edge


table(apply(adj.xyz.we[,-(1:px)]>0,2,sum))          ##number of kids for Y,Z net

table(apply(adj.xyz.we[-(1:px),]>0,1,sum))          ##number of parents

sum(adj.xyz.we[-(1:px),-(1:px)]>0)       ##total edge


##plot 
set.seed(3)
Plot_XYZ(result.xyz,hub.size=10, lay=layout.fruchterman.reingold)

#############################################
################### use  for DAG builting purpose 
#############################################
adj.result=t(adj.xyz.we>0)

true.dir=adj.result
true.ske=(adj.result+t(adj.result))>0

true.moral <- Moral(true.dir)
sum(true.moral[upper.tri(true.moral)]  == 3)
sum(true.moral[upper.tri(true.moral)]  == 2)
sum(true.moral[upper.tri(true.moral)]  == 1)

dim(true.moral)
library(ggplot2)
qplot(rowSums(true.moral))

hist(colSums(true.moral))
#save(true.dir,true.ske,file="dense_p500_adj_matrix.Rdata")
save(true.dir,true.ske,file="~/tmp/dense_p500_adj_matrix.Rdata")
