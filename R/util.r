## Clusterpath preprocessing
tri2vec <- cmpfun(function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
})

vec2tri <- cmpfun(function(k,n) {
  i <- ceiling(0.5*(2*n-1 - sqrt((2*n-1)^2 - 8*k)))
  j <- k - n*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
})

#' Construct indices matrices
#' 
#' \code{compactify_edges} constructs M1, M2, and ix index matrices. Note that storage conventions are different for ama and admm.
#' 
#' @param w weights vector
#' @param n number of points to cluster
#' @param method 'ama' or 'admm'
#' @export
#' @examples
#' p <- 10
#' n <- 20
#' X <- matrix(rnorm(p*n),p,n)
#' w <- kernel_weights(X,0)
#' w <- knn_weights(w,3,n)
#' method <- 'admm'
#' edge_info <- compactify_edges(w,n,method=method)
#' method <- 'ama'
#' edge_info <- compactify_edges(w,n,method=method)
compactify_edges <- cmpfun(function(w,n,method='ama') {
  if (!is.null(method) && !(method %in% c("ama","admm")))
    stop("method must be 'ama', 'admm', or NULL.")  
  sizes1 <- double(n)
  sizes2 <- double(n)
  
  if (method=='ama') {
    P <- vec2tri(which(w > 0),n)
  } else {
    P <- vec2tri(1:(n*(n-1)/2),n)
  }
  
  M1 <- matrix(0,nrow(P),n)
  M2 <- matrix(0,nrow(P),n)
  
  for (i in 1:n) {
    group1 <- which(P[,1] == i)
    sizes1[i] <- length(group1)
    if (sizes1[i] > 0) {
      M1[1:sizes1[i],i] <- group1
    }
    group2 <- which(P[,2] == i)
    sizes2[i] <- length(group2)
    if (sizes2[i] > 0) {
      M2[1:sizes2[i],i] <- group2
    }
  }
  
  M1 <- M1[1:max(sizes1),,drop=FALSE]
  M2 <- M2[1:max(sizes2),,drop=FALSE]
  
  return(list(ix=P,M1=M1,M2=M2,s1=sizes1,s2=sizes2))
})

#' Create adjacency matrix from V
#' 
#' \code{create_adjacency} creates an n-by-n sparse adjacency matrix from the matrix of centroid differences.
#' 
#' @param V Matrix of centroid differences
#' @param w Weights vector
#' @param ix Key mapping centroid difference indices to node pairs
#' @param n Number of points to cluster
#' @export
#' @examples
#' ## Clusterpaths for Mammal Dentition
#' data(mammals)
#' X <- as.matrix(mammals[,-1])
#' X <- t(scale(X,center=TRUE,scale=FALSE))
#' n <- ncol(X)
#' 
#' ## Pick some weights and a sequence of regularization parameters.
#' k <- 5
#' phi <- 0.5
#' w <- kernel_weights(X,phi)
#' w <- knn_weights(w,k,n)
#' gamma <- seq(0.0,43, length.out=100)
#' 
#' ## Perform clustering
#' nu <- AMA_step_size(w,n)
#' sol <- cvxclust_path_ama(X,w,gamma,nu=nu)
#' 
#' ## Construct adjacency matrix
#' A <- create_adjacency(sol$V[[10]],w,n)
#' G <- graph.adjacency(A, mode = 'upper')
#' plot(G,vertex.label=as.character(mammals[,1]),vertex.label.cex=0.65,vertex.label.font=2)
create_adjacency <- function(V,w,n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0)
  ix <- vec2tri(which(w>0),n)
  i <- ix[connected_ix,1]
  j <- ix[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}

#' Find clusters
#' 
#' \code{find_clusters} uses breadth-first search to identify the connected components of the corresponding
#' adjacency graph of the centroid differences vectors.
#' 
#' @param A adjacency matrix
#' @export
#' @examples
#' ## Clusterpaths for Mammal Dentition
#' data(mammals)
#' X <- as.matrix(mammals[,-1])
#' X <- t(scale(X,center=TRUE,scale=FALSE))
#' n <- ncol(X)
#' 
#' ## Pick some weights and a sequence of regularization parameters.
#' k <- 5
#' phi <- 0.5
#' w <- kernel_weights(X,phi)
#' w <- knn_weights(w,k,n)
#' gamma <- seq(0.0,43, length.out=100)
#' 
#' ## Perform clustering
#' nu <- AMA_step_size(w,n)
#' sol <- cvxclust_path_ama(X,w,gamma,nu=nu)
#' 
#' ## Construct adjacency matrix
#' A <- create_adjacency(sol$V[[10]],w,n)
#' find_clusters(A)
find_clusters <- function(A) {
  G <- graph.adjacency(A, mode = 'upper')
  n <- nrow(A)
  node_seen <- logical(n)
  cluster <- integer(n)
  k <- 1
  for (i in 1:n) {
    if (!node_seen[i]) {
      connected_set <- graph.bfs(G, root=i, unreachable = FALSE)$order
      node_seen[connected_set] <- TRUE
      cluster[connected_set] <- k
      k <- k + 1
    }
  }
  nClusters <- k - 1
  size <- integer(nClusters)
  for (j in 1:nClusters) {
    size[j] <- length(which(cluster == j))
  }
  return(list(cluster=cluster, size=size))
}
