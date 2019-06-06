
#' The Topic SCORE algorithm
#' 
#' This function obtains the word-topic matrix A from the word-document matrix X through the Topic SCORE algorithm.
#' @param K The number of topics.
#' @param X  The p-by-n word-document matrix, with each column being a distribution over a fixed set of vocabulary.
#' This matrix can be of class \code{simple_triplet_matrix} defined in \strong{slam} package, 
#' or any other class that can be transformed to class \code{dgRMatrix} defined in \strong{Matrix} package through
#' \code{as} function in \strong{methods} package.
#' @param K0 The number of greedy search steps in vertex hunting. If the value is missing it will be set to ceiling(1.5*K).
#' @param m  The number of centers in the kmeans step in vertex hunting. If the value is missing it will be set to 10*K.
#' @param Mquantile The percentage of the quantile of the diagonal entries of matrix M, 
#'  which is used to upper truncate the diagonal entries of matirx M.
#'  When it's 0, it will degenerate the case when there is no normalization.
#'  When it's 1, it means there is no truncation.
#'  Default is 0.
#' @param scatterplot  Whether a scatterplot of rows of R will be generated.
#' @param num_restart  The number of random restart in the kmeans step in vertex hunting.
#'  Default is 1.
#' @param seed The random seed. Default value is NULL.
#' @return A list containing \describe{
#'   \item{A_hat}{The estimated p-by-K word-topic matrix.}
#'   \item{R}{The p-by-(K-1) left singular vector ratios matrix.}
#'   \item{V}{The K-by-(K-1) vertices matrix, with each row being a vertex found through the vertex hunting algorithm in the simplex formed by the rows of R.}
#'   \item{Pi}{The p-by-K convex combinations matrix, with each row being the convex combination coefficients of a row of R using V as vertices.}
#'   \item{theta}{The K0-by-(K-1) matrix of K0 potential vertices found in the greedy step of the vertex hunting algorithm.}
#' }
#' @importFrom graphics par plot points
#' @importFrom utils install.packages installed.packages
#' @importFrom RSpectra svds
#' @importFrom methods as
#' @import Matrix 
#' @import slam
#' @importFrom stats kmeans quantile
#' @importFrom quadprog solve.QP
#' @importFrom utils combn
#' @export
#' @examples
#' data("AP")
#' K <- 3
#' tscore_obj <- topic_score(K, AP)
#' 
#' # Visualize the result
#' plot(tscore_obj$R[,1], tscore_obj$R[,2])
#' @author Minzhe Wang
#' @references Ke, Z. T., & Wang, M. (2017). A new SVD approach to optimal topic estimation. arXiv preprint arXiv:1704.07016.
#' @keywords models

topic_score <- function(K, X, K0, m, Mquantile=0, scatterplot=FALSE, num_restart = 1, seed=NULL){
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (missing(K0)){
    K0 <- ceiling(1.5*K)
  }
  if (missing(m)){
    m <- 10*K
  }
  
  if (class(X)=='simple_triplet_matrix'){
    X <- as.matrix(X)
  }
  X <- as(X, 'dgRMatrix')
  
  p <- dim(X)[1]
  n <- dim(X)[2]
  M <- rowMeans(X)
  M_trunk <- pmin(M,quantile(M,Mquantile))
  
  obj <- svds(sqrt(M_trunk^(-1))*X, K)
  Xi <- obj$u
  
  #Step 1
  Xi[,1] <- abs(Xi[,1])
  R <- apply(Xi[,2:K],2,function(x) x/Xi[,1])
  
  #Step 2
  vertices_est_obj <- vertices_est(R,K0,m,num_restart)
  V <- vertices_est_obj$V
  theta <- vertices_est_obj$theta
  
  if (scatterplot){
    par(mar=c(1,1,1,1))
    plot(R[,1],R[,2])
    points(V[,1],V[,2],col=2,lwd=5)
  }
  
  #Step 3
  Pi <- cbind(R, rep(1,p))%*%solve(cbind(V,rep(1,K)))
  Pi <- pmax(Pi,0)
  temp <- rowSums(Pi)
  Pi <- apply(Pi,2,function(x) x/temp)
  
  #Step 4
  A_hat <- sqrt(M_trunk)*Xi[,1]*Pi
  
  #Step 5
  temp <- colSums(A_hat)
  A_hat <- t(apply(A_hat,1,function(x) x/temp))
  
  return(list(A_hat=A_hat, R=R,V=V, Pi=Pi, theta=theta))
}




#' The vertex hunting in the Topic SCORE algorithm
#' 
#' This function conducts the vertex hunting in the Topic SCORE algorithm.
#' More generally this function finds a simplex with K vertices that best approximates
#' the given p data points in a (K-1) dimensional space.
#' @param R The p-by-(K-1) data matrix, with each row being a data point.
#' @param K0 The number of greedy search steps.
#' @param m  The number of centers in the kmeans step.
#' @param num_restart  The number of random start in the kmeans step.
#' @return A list containing \describe{
#'   \item{V}{The K-by-(K-1) vertices matrix, with each row being a vertex in the found simplex.}
#'   \item{theta}{The K0-by-(K-1) matrix of potential K0 vertices found in the greedy step.}
#' }
#' @importFrom graphics par plot points
#' @importFrom utils install.packages installed.packages
#' @import Matrix 
#' @import slam
#' @importFrom quadprog solve.QP
#' @importFrom utils combn
#' @export
#' @examples
#' # Generate 3 vertices
#' V <- rbind(c(-1/2,-1/2), c(1,0), c(0,1))
#' 
#' # Randomly generate the convex combination weights of 1000 points
#' Pi <- matrix(runif(3*1000),3,1000)
#' Pi <- apply(Pi, 2, function(x){x/sum(x)})
#' 
#' R <- t(Pi)%*%V
#' v_est_obj <- vertices_est(R, 1.5*3, 10*3, 1)
#' 
#' # Visualize the result
#' plot(R[,1], R[,2])
#' points(v_est_obj$V[,1], v_est_obj$V[,2], col=2, lwd=5)
#' @author Minzhe Wang
#' @references Ke, Z. T., & Wang, M. (2017). A new SVD approach to optimal topic estimation. arXiv preprint arXiv:1704.07016.


vertices_est <- function(R,K0,m,num_restart){
  K <- dim(R)[2] + 1
  
  #Step 2a
  obj <- kmeans(R,m,iter.max=100,nstart=num_restart)
  theta <- as.matrix(obj$centers)
  theta_original <- theta
  
  #Step 2b'
  inner <- theta%*%t(theta)
  distance <- diag(inner)%*%t(rep(1,m)) + rep(1,m)%*%t(diag(inner)) - 2*inner
  top2 <- which(distance==max(distance),arr.ind=TRUE)[1,]
  theta0 <- as.matrix(theta[top2,])
  theta <- as.matrix(theta[-top2,])
  
  if (K0 > 2){
    for (k0 in 3:K0){
      inner <- theta%*%t(theta)
      distance <- rep(1,k0-1)%*%t(diag(inner))-2*theta0%*%t(theta)
      ave_dist <- colMeans(distance)
      index <- which(ave_dist==max(ave_dist))[1]
      theta0 <- rbind(theta0, theta[index,])
      theta <- as.matrix(theta[-index,])
    }
    theta <- theta0
  }
  
  #Step 2b
  comb <- combn(1:K0, K)
  max_values <- rep(0, dim(comb)[2])
  for (i in 1:dim(comb)[2]){
    for (j in 1:K0){
      max_values[i] <- max(simplex_dist(as.matrix(theta[j,]), as.matrix(theta[comb[,i],])), max_values[i])
    }
  }
  
  min_index <- which(max_values == min(max_values))
  
  return(list(V=theta[comb[,min_index[1]],], theta=theta_original))
}



#' The l_2 distance between a point and a simplex
#' 
#' This function computes the l_2 distance between a point and a simplex, 
#' that is the shortest l_2 distance between the given point and any point in the simplex.
#' @param theta A (K-1) dimensional vector, representing a point.
#' @param V The K-by-(K-1) vertices matrix, with each row being a vertex.
#' @return The l_2 distance.
#' @importFrom graphics par plot points
#' @importFrom utils install.packages installed.packages
#' @import Matrix 
#' @import slam
#' @importFrom quadprog solve.QP
#' @export
#' @examples
#' # Generate 3 vertices
#' V <- rbind(c(-1/2,-1/2), c(1,0), c(0,1))
#' 
#' theta <- c(3,1)
#' simplex_dist(theta, V)
#' @author Minzhe Wang
#' @references Ke, Z. T., & Wang, M. (2017). A new SVD approach to optimal topic estimation. arXiv preprint arXiv:1704.07016.


simplex_dist <- function(theta, V){
  VV <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))%*%V
  M <- VV%*%t(VV)
  d <- VV%*%(theta-V[dim(V)[1],])
  
  A <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))
  b0 <- c(rep(0,dim(V)[1]-1),-1)
  
  obj <- solve.QP(M, d, A, b0)
  return( sqrt( max(sum((theta-V[dim(V)[1],])^2)+ 2*obj$value, 0) ) )
}

