#' Compute step size Anderson-Morely upper bound on the largest eigenvalue of the Laplacian
#' 
#' \code{AMA_step_size} computes a step size based on the better of two bounds derived by Anderson 
#' and Morely.
#' 
#' @param w vector of weights
#' @param p number of points to cluster
#' @export
#' @examples
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
#' AMA_step_size(w,n)
AMA_step_size <- cmpfun(function(w,n,inflate=1.9999) {
  ix <- compactify_edges(w,n)$ix
  nEdges <- nrow(ix)
  nVertex <- max(ix)
  deg <- integer(nVertex)
  for (i in 1:nVertex) {
    deg[i] <- length(which(ix[,1] == i)) + length(which(ix[,2] == i))
  }
  bound <- -Inf
  for (j in 1:nEdges) {
    bound <- max(bound,deg[ix[j,1]] + deg[ix[j,2]])
  }
  return(inflate/min(bound,nVertex))
})