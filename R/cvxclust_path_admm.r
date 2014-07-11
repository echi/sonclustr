#' Convex Clustering Path via ADMM
#' 
#' \code{cvxclust_path_admm} estimates the convex clustering path via ADMM.
#' Required inputs include a data matrix \code{X} (rows are features; columns are samples), a vector of weights
#' \code{w}, and a sequence of regularization parameters \code{gamma}.
#' Two penalty norms are currently supported: 1-norm and 2-norm. ADMM admits acceleration by extrapolated steps akin to those in FISTA.
#' This speed-up is employed by default.
#' 
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param gamma A sequence of regularization parameters.
#' @param nu A positive penalty parameter for quadratic deviation term.
#' @param tol_abs The convergence tolerance (absolute).
#' @param tol_rel The convergence tolerance (relative).
#' @param max_iter The maximum number of iterations.
#' @param accelerate If \code{TRUE} (the default), acceleration is turned on.
#' @return \code{U} A list of centroid matrices.
#' @return \code{V} A list of centroid difference matrices.
#' @return \code{Lambda} A list of Lagrange multiplier matrices.
#' @export
#' @author Eric C. Chi, Kenneth Lange
#' @seealso \code{\link{cvxclust_path_ama}} for estimating the clustering path with AMA. \code{\link{kernel_weights}} and \code{\link{knn_weights}} compute useful weights.
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
#' sol <- cvxclust_path_admm(X,w,gamma)
#' 
#' ## Plot the cluster path
#' library(ggplot2)
#' svdX <- svd(X)
#' pc <- svdX$u[,1:2,drop=FALSE]
#' pc.df <- as.data.frame(t(pc)%*%X)
#' nGamma <- sol$nGamma
#' df.paths <- data.frame(x=c(),y=c(), group=c())
#' for (j in 1:nGamma) {
#'   pcs <- t(pc)%*%sol$U[[j]]
#'   x <- pcs[1,]
#'   y <- pcs[2,]
#'   df <- data.frame(x=pcs[1,], y=pcs[2,], group=1:n)
#'   df.paths <- rbind(df.paths,df)
#' }
#' X_data <- as.data.frame(t(X)%*%pc)
#' colnames(X_data) <- c("x","y")
#' X_data$Name <- mammals[,1]
#' data_plot <- ggplot(data=df.paths,aes(x=x,y=y))
#' data_plot <- data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
#' data_plot <- data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name), position=position_jitter(h=0.125,w=0.125))
#' data_plot <- data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
#' data_plot <- data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
#' data_plot + theme_bw()
cvxclust_path_admm <- function(X,w,gamma,nu=1,tol_abs=1e-5,tol_rel=1e-4,max_iter=1e4,accelerate=TRUE) {
  call <- match.call()
  nGamma <- length(gamma)
  n <- ncol(X)
  p <- nrow(X)
  edge_info <- compactify_edges(w,n,method='admm')
  nK <- n*(n-1)/2
  Lambda <- matrix(0,p,nK)
  V <- matrix(0,p,nK)
  list_U <- vector(mode="list",length=nGamma)
  list_V <- vector(mode="list",length=nGamma)
  list_Lambda <- vector(mode="list",length=nGamma)
  ix <- edge_info$ix  
  M1 <- edge_info$M1
  M2 <- edge_info$M2
  s1 <- edge_info$s1
  s2 <- edge_info$s2
  iter_vec <- integer(nGamma)  
#  print("gamma    its | primal res      dual res      max res     ")
#  print("---------------------------------------------------------")    
  for (ig in 1:nGamma) {
    gam <- gamma[ig]
    cc <- cvxclust_admm(X,Lambda,ix-1,M1-1,M2-1,s1,s2,w,gam,nu=nu,max_iter=max_iter,tol_abs=tol_abs,tol_rel=tol_rel,accelerate=accelerate)
    iter_vec[ig] <- cc$iter
    Lambda <- cc$Lambda
    V <- cc$V
    list_U[[ig]] <- cc$U
    list_V[[ig]] <- V
    list_Lambda[[ig]] <- Lambda 
 #   print(sprintf("%5d  %5d | %5f        %5f      %7f", ig, cc$iter, signif(cc$primal[cc$iter],4),
 #                 signif(cc$dual[cc$iter],4),
 #                 signif(max(cc$primal[cc$iter],cc$dual[cc$iter]),4)))    
    #    if (norm(cc$V,'1')==0) {
    #      print('Single cluster')
    #      break
    #    }
  }
  cvxclust_obj <- list(U=list_U,V=list_V,Lambda=list_Lambda,nGamma=ig,iters=iter_vec,call=call)
  class(cvxclust_obj) <- "cvxclustobject"
  return(cvxclust_obj)  
}
