#' Create a random clustering problem
#' 
#' \code{create_clustering_problem} makes a random clustering problem for testing purposes.
#' 
#' @param p Dimension of space of points to be clustered
#' @param n Number of points
#' @param seed Random number seed
#' @param nnn Number of nearest neighbors
#' @param method 'ama' or 'admm'
#' @export
#' @examples
#' p <- 10
#' n <- 20
#' seed <- 12345
#' rnd_problem_admm <- create_clustering_problem_new(p,n,seed)
create_clustering_problem <- function(p,n,seed=12345,nnn=3,method='ama') {
  if (!is.null(method) && !(method %in% c("ama","admm")))
    stop("method must be 'ama', 'admm', or NULL.")    
  set.seed(seed)
  X <- matrix(rnorm(p*n),p,n)
  w <- kernel_weights(X,0)
  w <- knn_weights(w,nnn,n)
  if (method=='ama') {
    w <- w[w>0]
  }
  edge_info <- compactify_edges(w,n,method=method)
  ix <- edge_info$ix
  M1 <- edge_info$M1
  M2 <- edge_info$M2
  s1 <- edge_info$s1
  s2 <- edge_info$s2
  ix <- edge_info$ix
  return(list(X=X,ix=ix-1,M1=M1-1,M2=M2-1,s1=s1,s2=s2,w=w))
}
