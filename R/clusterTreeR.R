#' clusterTreeR
#' @export
clustertree_ex <- function(x, k = 5L){
  suppressMessages({ require("igraph"); require("FNN") })

  ## Initialize variables
  { n <- nrow(x); dist_x <- dist(x); }

  ## Get the smallest radius as a starting point for r
  r <- min(dist_x)

  ## Start off with empty adjacency graph
  G_r <- igraph::graph_from_adjacency_matrix(diag(n))

  ## Get balls centered at each x of radius r_k
  ## (containing k points inclusive of x_i itself)
  r_k <- apply(FNN::knn.dist(x, k = k - 1, algorithm = "kd_tree"), 1, max)

  ## Initiate cluster tree and distance matrix to simplify indexing
  clustertree <- list()
  l2_dist <- as.matrix(dist_x)
  diag(l2_dist) <- Inf

  ## Vector of sorted radii to iterate through and a counter
  lambda <- sort(dist_x)
  alpha <- 1
  i <- 1

  ## expand eps-Ball from 0 -> Inf
  for (r in lambda){

    ## Wisharts scheme: Only connect points that have at least
    ## k neighbors within distance r
    eps <- mapply(function(i, j) ( l2_dist[i, j] <= r * alpha &&
                                     r_k[i] <= r && # radius of x_i
                                     r_k[j] <= r ), # radius of x_j
                  row(l2_dist), col(l2_dist))

    ## Construct the adjacency graph, recording distinct CCs as the level sets
    adj_matrix <- matrix(as.integer(eps), nrow = n, ncol = n)
    G_r <- igraph::graph_from_adjacency_matrix(adj_matrix)
    CC <- igraph::components(G_r)$membership

    ## Record distinct level sets
    if (i == 1 || any(CC != clustertree[[i-1]]$cluster)){
      clustertree[[i]] <- list(cluster=CC, radius=r * alpha)
      i <- i + 1
    }
    if (length(clustertree) == n - 1) break
  }
  return(clustertree)
}
