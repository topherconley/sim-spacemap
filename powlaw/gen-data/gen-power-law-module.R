
########################################################
#Function to generate a power law module 
gen_powlaw_module <- function(id, R, S, ntophubs, prop_xhubs, 
                              xbg_clusters = 5, prob_xbg_edge = 0.3) { 
  library(igraph)
  #out-degree of new vertices being attached
  out.seq <- pmax(1,rpois(n = R, lambda = 2.2)) #2.2
  #simulate a Barbasi-Albert scale-free preferential attachment network with R nodes
  pag <- sample_pa(n = R, power = 1.3, out.seq =  out.seq, out.pref = TRUE, directed = FALSE)
  #the number of edges
  #ecount(pag)
  #the number of vertices
  #vcount(pag)
  deg_pag <- degree(pag)
  #degree distribution summary
  #summary(deg_pag)
  #boxplot(degree(pag))
  #Inspect the power law fit
  plf <- power.law.fit(x = deg_pag)
  
  #Identify hubs for X and Y
  hub_pag <- hub.score(pag, scale = TRUE, weights=NULL, options = igraph.arpack.default)
  #plot(hub_pag$vector, deg_pag, pch = 19)
  #top hubs split in half among X and Y
  tophubs_index <- sort(hub_pag$vector, index.return = TRUE, decreasing = TRUE)$ix[1:ntophubs]
  xhubs_index <- sample(tophubs_index, size = ceiling(length(tophubs_index)*prop_xhubs))
  yhubs_index <- setdiff(tophubs_index, xhubs_index)
  #hist(degree(pag)[xhubs_index])
  #hist(degree(pag)[yhubs_index])
  
  #Identify the non-hub Y nodes
  ybg_index <- setdiff(1:R, c(xhubs_index, yhubs_index))
  
  #How many x-x hub edges exist? 
  #very few edges  ut of the whole graph  and the total X hub edges, okay to efffectively delete X--X hub edges
  xx_edges <- E(pag)[ xhubs_index %--% xhubs_index]
  #delete these edges from the graph as we are trying to create no geographic (genomic) correlation between true hubs
  cgpag <- delete_edges(pag, xx_edges)
  #ecount(cgpag)
  #new degrees for X hubs
  #degree(cgpag, v = xhubs_index)
  
  #Still remains a power-law with 2 < alpha < 3 in the range expected for biological networks
  noxx_plf <- power.law.fit(x = degree(cgpag))
  #noxx_plf$alpha
  
  #Assert the module has no disconnected components among the original power-law distribution
  stopifnot(count_components(cgpag) == 1)
  
  #add additional X vertices of size @param<S> which will have some edges between them. 
  cgpag <- add_vertices(graph = cgpag, nv = S)
  xbg_index <- (R + 1):(R + S)
  
  #Create X--X subgraphs of size @param<xbg_clusters> from an Erdos-Reny model, which has the potential to include an Xhub. 
  start <- R + 1
  ostart <- start
  for (xhub in xhubs_index) { 
    xid <- c(xhub, start:(start + xbg_clusters - 1))
    tmpg <- sample_gnp(n = length(xid), p  = prob_xbg_edge)
    tmpg <- set_vertex_attr(graph = tmpg, name = "name", value = xid)
    cgpag <- add_edges(cgpag, edges = as.vector(t(get.edgelist(tmpg))))
    start <- start + 5
  }
  #set of X's in X clusters with potential connection to xhubs
  xset1 <- ostart:(start - 1)
  
  #Create X--X subgraphs of size @param<xbg_clusters> from an Erdos-Reny model but do not include an X hub. 
  for(xbs in seq(from = start, to = (R + S), by = xbg_clusters)) {
    xid <- xbs:(xbs + xbg_clusters - 1)
    tmpg <- sample_gnp(n = length(xid), p  = prob_xbg_edge)
    tmpg <- set_vertex_attr(graph = tmpg, name = "name", value = xid)
    cgpag <- add_edges(cgpag, edges = as.vector(t(get.edgelist(tmpg))))
  }
  #set of X's in X clusters with no connection to xhubs
  xset2 <- start:(R + S)
  
  #unique vertex id
  #set_vertex_attr(cgpag, name = "name", index = V(cgpag), value = paste(id, V(cgpag), sep = ":"))
  V(cgpag)$name <- paste(id, V(cgpag), sep = ":")
  #cgpag <- set_vertex_attr(cgpag, name = "module_id", index = V(cgpag), value = id)
  list(mod = cgpag, 
       xhubs = paste(id, xhubs_index, sep = ":"), 
       yhubs = paste(id, yhubs_index, sep = ":"),
       xbg  = paste(id, xbg_index, sep = ":"),
       ybg = paste(id, ybg_index, sep = ":"), 
       xsubhubs = paste(id, xset1, sep = ":"), 
       xsubbg = paste(id, xset2, sep = ":"))
}

