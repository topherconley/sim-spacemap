
dag_xy_xhub <- function(px, py, nxhub = 20, min_hub_size = 5, mean_hub_size = 12) { 
  
  #construct x hubs with y targets ( x->y )
  index_xhubs <- sample(px, nxhub)
  xy_pool <- expand.grid(x = index_xhubs, y = seq_len(py))
  hub_size_pool <- pmax(rpois(1000, mean_hub_size),min_hub_size)
  #hist(hub_size_pool)
  y_pool <- seq_len(py)
  samp_y_edges <- function(x) { 
    sample(x = y_pool, size = sample(hub_size_pool,1))
  }
  xy_edges <- lapply(index_xhubs, samp_xy_edges)
  
  x2y_adj <- matrix(0,px,py)
  for (i in seq_along(index_xhubs)) x2y_adj[index_xhubs[i],xy_edges[[i]]] <- 1
  x2y_adj
}


x2y_adj <- dag_xy_xhub(px = 250, py = 300)
dim(x2y_adj)
