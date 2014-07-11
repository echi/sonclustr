#' "Thin" a weight vector to be positive only for its k-nearest neighbors
#' 
#' \code{knn_weights} takes a weight vector \code{w} and sets the ith 
#' component \code{w[i]} to zero if either of the two corresponding nodes
#' is not among the other's \code{k} nearest neighbors.
#' 
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param k The number of nearest neighbors
#' @param p The number of data points.
#' @author Eric C. Chi, Kenneth Lange
#' @export
#' @return A vector \cite{w} of weights for convex clustering.
knn_weights <- function(w,k,p) {
  i <- 1
  neighbors <- tri2vec(i,(i+1):p,p)
  keep <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(p-1)) {
    group_A <- tri2vec(i,(i+1):p,p)
    group_B <- tri2vec(1:(i-1),i,p)
    neighbors <- c(group_A,group_B)
    knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep <- union(knn,keep)
  }
  i <- p
  neighbors <- tri2vec(1:(i-1),i,p)
  knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep <- union(knn,keep)
  w[-keep] <- 0
  return(w)
}

#' Compute Gaussian Kernel Weights
#' 
#' \code{kernel_weights} computes Gaussian kernel weights given a data matrix \code{X} and a scale parameter \code{phi}. Namely,
#' the lth weight \code{w[l]} is given by
#' \deqn{
#' w[l] = exp(-phi ||X[,i]-X[,j]||^2)
#' }, where the lth pair of nodes is (\code{i},\code{j}).
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param phi The nonnegative parameter that controls the scale of kernel weights
#' @author Eric C. Chi, Kenneth Lange
#' @export
#' @return A vector \cite{w} of weights for convex clustering.
kernel_weights <- function(X,phi) {
  p <- ncol(X)
  k <- 1
  w <- matrix(0,p*(p-1)/2,1)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      w[k] <- exp(-phi*norm(as.matrix(X[,i,drop=FALSE]-X[,j,drop=FALSE]),'f')^2)
      k <- k+1
    }
  }
  return(weights=w)
}