#'general_CG
#'
#'Conduct general Conjugate Gradient Algorithm to find numeric roots of coupled matrix equation sets AXB = C
#'
#'@param A input matrixes list A with length k^2
#'
#'@param B input matrixes list B with length k^2
#'
#'@param C input matrices list C with length k
#'
#'@param tol tolerance of residuals from 0 per X element, when res < tol, return current X as the numeric root
#'
#'@param max.iter the maxium iteration before ending the for loop, with default 100000
#'
#'@return solve of matrix function
#'
#'@examples
#'A11 <- matrix(runif(12),3,4)
#'A12 <- matrix(runif(6),3,2)
#'A13 <- matrix(runif(15),3,5)
#'A21 <- matrix(runif(16),4,4)
#'A22 <- matrix(runif(8),4,2)
#'A23 <- matrix(runif(20),4,5)
#'A31 <- matrix(runif(8),2,4)
#'A32 <- matrix(runif(4),2,2)
#'A33 <- matrix(runif(10),2,5)
#'B11 <- matrix(runif(12),4,3)
#'B12 <- matrix(runif(9),3,3)
#'B13 <- matrix(runif(6),2,3)
#'B21 <- matrix(runif(8),4,2)
#'B22 <- matrix(runif(6),3,2)
#'B23 <- matrix(runif(4),2,2)
#'B31 <- matrix(runif(12),4,3)
#'B32 <- matrix(runif(9),3,3)
#'B33 <- matrix(runif(6),2,3)
#'C1 <- matrix(c(1:9),3,3)
#'C2 <- matrix(c(-4:3),4,2)
#'C3 <- matrix(c(-1:4),2,3)
#'C <- list(C1,C2,C3)
#'A <- list(A11,A12,A13,A21,A22,A23,A31,A32,A33)
#'B <- list(B11,B12,B13,B21,B22,B23,B31,B32,B33)
#'X <- general_CG(A,B,C)
#'#Verification
#'A11%*%X[[1]][[1]]%*%B11 +
#'A12%*%X[[1]][[2]]%*%B12 + A13%*%X[[1]][[3]]%*%B13
#'C1
#'A21%*%X[[1]][[1]]%*%B21 +
#'A22%*%X[[1]][[2]]%*%B22 + A23%*%X[[1]][[3]]%*%B23
#'C2
#'A31%*%X[[1]][[1]]%*%B31 +
#'A32%*%X[[1]][[2]]%*%B32 + A33%*%X[[1]][[3]]%*%B33
#'C3
#'
#'@export
#'
general_CG <- function(A,B,C,tol= 1e-20,max.iter = 100000){
  if (length(A) != length(B) | length(C)^2 != length(A)){
    stop("Please Check the number of matrixes of A, B and C!")
  }
  dim_X <- dimCalculator(A,B,C)
  X <- vector("list",length = dim(dim_X)[2])
  for (i in 1:length(X)) {
    X[[i]] <- matrix(0,dim_X[1,i],dim_X[2,i])
  }
  zeros_X <- X
  save_sum <- C
  length_A <- length(A)
  length_X <- sqrt(length_A)
  for (i in 1:length_A){
    mod_X <- (i-1)%%length_X+1
    group_X <- (i-1)%/%length_X+1
    test <- crossprod(t(A[[i]]),X[[mod_X]]%*%B[[i]])
    save_sum[[group_X]] <- save_sum[[group_X]] - test
  }
  R1 <- unlist(save_sum)
  prepare <- rep(X,length_X)
  p <- X
  for (i in 1:length_A){
    group_X <- (i-1)%/%length_X+1
    mod_X <- (i-1)%%length_X+1
    prepare[[i]] <-crossprod(A[[i]],save_sum[[group_X]]%*%t(B[[i]]))
    p[[mod_X]] <- p[[mod_X]] + prepare[[i]]
  }

  ## Step 5

  new_X <- zeros_X
  for (iter in 1:max.iter){
    for (i in 1:length_X){
      new_X[[i]] <- X[[i]] + sum(R1^2)/sum(unlist(p)^2)*p[[i]]
    }
    new_save_sum <- C
    for (i in 1:length_A){
      mod_X <- (i-1)%%length_X+1
      group_X <- (i-1)%/%length_X+1
      test <- crossprod(t(A[[i]]),new_X[[mod_X]]%*%B[[i]])
      new_save_sum[[group_X]] <- new_save_sum[[group_X]] - test
    }
    new_R <- unlist(new_save_sum)
    for (i in 1:length_X){
      p[[i]] <- sum(new_R^2)/sum(R1^2)*p[[i]]
    }
    prepare <- rep(zeros_X,length_X)

    for (i in 1:length_A){
      group_X <- (i-1)%/%length_X+1
      mod_X <- (i-1)%%length_X+1
      prepare[[i]] <-crossprod(A[[i]],new_save_sum[[group_X]]%*%t(B[[i]]))
      p[[mod_X]] <- p[[mod_X]] + prepare[[i]]
    }
    if (sum(new_R^2)/length(new_R)<tol){
      break
    }
    X <- new_X
    save_sum <- new_save_sum
    R1 <- new_R
  }
  return(list(X,iter))
}


 # ## Verification
 # A11%*%X[[1]][[1]]%*%B11 + A12%*%X[[1]][[2]]%*%B12 + A13%*%X[[1]][[3]]%*%B13
 # C1
 # A21%*%X[[1]][[1]]%*%B21 + A22%*%X[[1]][[2]]%*%B22 + A23%*%X[[1]][[3]]%*%B23
 # C2
 # A31%*%X[[1]][[1]]%*%B31 + A32%*%X[[1]][[2]]%*%B32 + A33%*%X[[1]][[3]]%*%B33
 # C3

