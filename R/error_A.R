#' The l_1 distance between two thin matrices up to a column permuation
#' 
#' This function computes l_1 distance between two thin matrices up to a column permuation, 
#' that is to find the smallest sum of absolute value entry-wise difference between two matrices 
#' over all possible permutations over the columns of the first matrix. This can be done
#' either universially or greedily.  
#' @param A The first p-by-K matrix.
#' @param A_hat The second p-by-K matrix.
#' @param type The search type for the best permutation. If it's 'u', the search is 
#'   done universially, that is over all possible permuations of the columns of A.
#'   If it's 'g', the search is done greedily, that is at kth step find the closest 
#'   column in the remaining columns of A to the kth column of A_hat in terms of l_1 distance. Greedy search may
#'   result in sub-optimal solutions, but it can be computed much faster than
#'   the universal way when K is large. The default value is 'u'.
#' @return The l_1 distance.
#' @importFrom graphics par plot points
#' @importFrom utils install.packages installed.packages
#' @importFrom combinat permn
#' @export
#' @examples
#' # The example uses the runif() function in the 'stats' package
#' A <- matrix(runif(10*3),10,3)
#' A_hat <- A + 0.1*matrix(runif(10*3),10,3)
#' error_A(A, A_hat)
#' error_A(A, A_hat, type='g')
#' @author Minzhe Wang

error_A <- function(A, A_hat, type='u'){
  K <- dim(A)[2]
  
  if (type=='u'){
    all_perm <- combinat::permn(1:K)
    error <- Inf
    
    for (i in 1:length(all_perm)){
      error <- min(error, sum(colSums(abs(A[,all_perm[[i]]]-A_hat))))
    }
  }
  else{
    used <- rep(1,K)
    A_perm <- matrix(0,dim(A)[1],dim(A)[2])
    
    for (k in 1:K){
      dis <- colSums(abs(A-A_hat[,k]))*used
      index <- which(dis == min(dis))
      index <- index[1]
      A_perm[,k] <- A[,index]
      used[index] <- Inf
    }
    
    error <- sum(colSums(abs(A_perm-A_hat)))
  }
  return(error)
}

