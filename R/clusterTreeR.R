#' #' clusterTreeR
#' clustertree_ex <- function(x, k = 5L, alpha = sqrt(2), type = 0){
#'   suppressMessages({ require("igraph"); require("FNN") })
#'
#'   ## Initialize variables
#'   { n <- nrow(x); dist_x <- dist(x); }
#'
#'   ## Start off with empty adjacency graph
#'   G_r <- igraph::graph_from_adjacency_matrix(diag(n))
#'
#'   ## Get balls centered at each x of radius r_k
#'   ## (containing k points inclusive of x_i itself)
#'   r_k <- apply(dbscan::kNNdist(x, k = k - 1), 1, max)
#'
#'   ## Initiate cluster tree and distance matrix to simplify indexing
#'   clustertree <- list()
#'   l2_dist <- as.matrix(dist_x)
#'   diag(l2_dist) <- Inf
#'
#'   ## Vector of sorted radii to iterate through and a counter
#'   lambda <- sort(dist_x)
#'
#'   ## Initialize with every point as a singleton
#'   clustertree[[1]] <- list(cluster=1:n, radius=0, iter = 0)
#'   i <- 2
#'
#'   iter <- 0
#'   ## expand eps-Ball from 0 -> Inf
#'   pb <- txtProgressBar(max = n - 1L, style = 3)
#'
#'   ## Get the from <--> to indices
#'   from <- row(l2_dist)[lower.tri(l2_dist)]
#'   to <- col(l2_dist)[lower.tri(l2_dist)]
#'   mo <- cbind(from, to)
#'   mo <- mo[order(dist_x),]
#'
#'   ## Calculate R based on type
#'
#'   for (r in lambda){
#'
#'     ## Wisharts scheme: Only connect points that have at least
#'     ## k neighbors within distance r
#'     eps <- mapply(function(i, j) ( l2_dist[i, j] <=  alpha * max(c(r_k[i], r_k[j])) &&
#'                                      r_k[i] <= r && # radius of x_i
#'                                      r_k[j] <= r ), # radius of x_j
#'                   row(l2_dist), col(l2_dist))
#'
#'     ## Construct the adjacency graph, recording distinct CCs as the level sets
#'     adj_matrix <- matrix(as.integer(eps), nrow = n, ncol = n)
#'     G_r <- igraph::graph_from_adjacency_matrix(adj_matrix)
#'     CC <- igraph::components(G_r)$membership
#'     CC_admitted <- CC
#'     CC_admitted[which(r_k > r)] <- 0
#'
#'     # get.edgelist(G_r)[unique(c(which(get.edgelist(G_r)[, 1] %in% c(67, 78)), which(get.edgelist(G_r)[, 2] %in% c(67,78)))),]
#'     ## Record distinct level sets
#'     if (any(CC != clustertree[[i-1]]$cluster)){
#'       clustertree[[i]] <- list(cluster=CC, radius=r * alpha, dist_ij = r, iter = iter, CC = CC_admitted,
#'                                n_admitted = sum(r_k <= r))
#'       i <- i + 1
#'     }
#'     setTxtProgressBar(pb, value = i)
#'     if (length(clustertree) == n) break
#'     iter <- iter + 1
#'   }
#'   close(pb)
#'   attr(clustertree, "call") <- match.call()
#'   return(clustertree)
#' }
