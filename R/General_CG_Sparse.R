#'general_CG_sparse
#'
#'@param A input matrixes
#'
#'@param B input matrixes
#'
#'
#'@return solve of matrix function
#'
#'@examples
#'A11 <- matrix(runif(12),3,4)
#'A11 <- Matrix(A11,sparse = T)
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
#'general_CG(A,B,C)
#'
#'@export
#'
general_CG_sparse <- function(A,B,C,tol= 1e-20,max.iter = 10000){
  if (length(A) != length(B) | length(C)^2 != length(A)){
    stop("Please Check the number of matrixes of A, B and C!")
  }
  dim_X <- dimCalculator(A,B,C)
  num_X <- sum(dim_X[1,]*dim_X[2,])
  X <- vector("list",length = dim(dim_X)[2])
  for (i in 1:length(X)) {
    X[[i]] <- matrix(0,dim_X[1,i],dim_X[2,i])
    X[[i]] <- Matrix(X[[i]],sparse = T)
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
  R1 <- c()
  for (i in 1:length_X){
    R1 <- c(R1,as.vector(save_sum[[i]]))
  }
  prepare <- rep(X,length_X)
  p <- X
  for (i in 1:length_X){
    p[[i]] <- as.matrix(X[[i]])
  }
  for (i in 1:length_A){
    group_X <- (i-1)%/%length_X+1
    mod_X <- (i-1)%%length_X+1
    prepare[[i]] <-t(A[[i]])%*%save_sum[[group_X]]%*%t(B[[i]])
    p[[mod_X]] <- as.matrix(p[[mod_X]] + prepare[[i]])
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
      test <- crossprod(t(A[[i]]),crossprod(t(new_X[[mod_X]]),B[[i]]))
      new_save_sum[[group_X]] <- new_save_sum[[group_X]] - test
    }
    new_R <- c()
    for (i in 1:length_X){
      new_R <- c(new_R,as.vector(new_save_sum[[i]]))
    }
    for (i in 1:length_X){
      p[[i]] <- sum(new_R^2)/sum(R1^2)*p[[i]]
    }
    prepare <- rep(zeros_X,length_X)

    for (i in 1:length_A){
      group_X <- (i-1)%/%length_X+1
      mod_X <- (i-1)%%length_X+1
      prepare[[i]] <-crossprod(A[[i]],tcrossprod(new_save_sum[[group_X]],B[[i]]))
      p[[mod_X]] <- as.matrix(p[[mod_X]] + prepare[[i]])
    }
    if (sum(new_R^2)/length(new_R)<tol){
      break
    }
    X <- new_X
    save_sum <- new_save_sum
    R1 <- new_R
  }
  return(X)
}

A11 <- Matrix(A11,sparse = T)
A12 <- Matrix(A12,sparse = T)
A13 <- Matrix(A13,sparse = T)
A21 <- Matrix(A21,sparse = T)
A22 <- Matrix(A22,sparse = T)
A23 <- Matrix(A23,sparse = T)
A31 <- Matrix(A31,sparse = T)
A32 <- Matrix(A32,sparse = T)
A33 <- Matrix(A33,sparse = T)
B11 <- Matrix(B11,sparse = T)
B12 <- Matrix(B12,sparse = T)
B13 <- Matrix(B13,sparse = T)
B21 <- Matrix(B21,sparse = T)
B22 <- Matrix(B22,sparse = T)
B23 <- Matrix(B23,sparse = T)
B31 <- Matrix(B31,sparse = T)
B32 <- Matrix(B32,sparse = T)
B33 <- Matrix(B33,sparse = T)
C1 <- Matrix(C1,sparse = T)
C2 <- Matrix(C2,sparse = T)
C3 <- Matrix(C3,sparse = T)
C <- list(C1,C2,C3)
A <- list(A11,A12,A13,A21,A22,A23,A31,A32,A33)
B <- list(B11,B12,B13,B21,B22,B23,B31,B32,B33)

R1
