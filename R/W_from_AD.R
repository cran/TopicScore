#' Estimation of W from A and X
#' 
#' This function estimates the topic-document matrix W from the word-topic matrix A and the word-document X through quadratic programming.
#' @param A The p-by-K word-topic matrix.
#' @param X The p-by-n word-document matrix.
#' @return The estimated K-by-n topic-document matrix W_hat.
#' @importFrom graphics par plot points
#' @importFrom utils install.packages installed.packages
#' @import Matrix 
#' @import slam
#' @importFrom quadprog solve.QP
#' @export
#' @examples
#' data("AP")
#' K <- 3
#' tscore_obj <- topic_score(K, AP)
#' W_hat <- W_from_AD(tscore_obj$A_hat, AP)
#' @author Minzhe Wang

W_from_AD <- function(A, X){
  K <-dim(A)[2]
  n <- dim(X)[2]
  
  W_hat <- matrix(0, K, n)
  M <- rbind(diag(K-1), rep(-1,K-1))
  bM <- diag(K)[,K]
  Dmat <- 2*t(A%*%M)%*%(A%*%M)
  Amat <- t(M)
  bvec <- -bM
  
  AM <- A%*%M
  AbM <- A%*%bM
  
  if (class(X)!="matrix"){
    X <- as.matrix(X)
  }
  
  for (i in 1:n){
    dvec <- 2*t(X[,i]-AbM)%*%AM
    qp_sol <- solve.QP(Dmat, dvec, Amat, bvec)$solution
    W_hat[,i] <- c(qp_sol, 1-sum(qp_sol))
  }
  W_hat <- pmax(W_hat,0)
  
  return(W_hat)
}
