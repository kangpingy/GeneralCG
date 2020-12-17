#'general_CG
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
#'X1 <- matrix(0,dim(A11)[2],dim(B11)[1])
#'X2 <- matrix(0,dim(A12)[2],dim(B12)[1])
#'X3 <- matrix(0,dim(A13)[2],dim(B13)[1])
#'zeros_X <- list(X1,X2,X3)
#'X <- list(X1,X2,X3)
#'general_CG(A,B,C)
#'
#'@export
#'
#'
general_CG <- function(A,B,C,tol= 1e-20,max.iter = 10000){
  save_sum <- C
  length_A <- length(A)
  length_X <- sqrt(length_A)
  for (i in 1:length_A){
    mod_X <- (i-1)%%length_X+1
    group_X <- (i-1)%/%length_X+1
    test <- A[[i]]%*%X[[mod_X]]%*%B[[i]]
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
      test <- A[[i]]%*%new_X[[mod_X]]%*%B[[i]]
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
  return(X)
}
A11 <- matrix(runif(12),3,4)
A12 <- matrix(runif(6),3,2)
A13 <- matrix(runif(15),3,5)
A21 <- matrix(runif(16),4,4)
A22 <- matrix(runif(8),4,2)
A23 <- matrix(runif(20),4,5)
A31 <- matrix(runif(8),2,4)
A32 <- matrix(runif(4),2,2)
A33 <- matrix(runif(10),2,5)
B11 <- matrix(runif(12),4,3)
B12 <- matrix(runif(9),3,3)
B13 <- matrix(runif(6),2,3)
B21 <- matrix(runif(8),4,2)
B22 <- matrix(runif(6),3,2)
B23 <- matrix(runif(4),2,2)
B31 <- matrix(runif(12),4,3)
B32 <- matrix(runif(9),3,3)
B33 <- matrix(runif(6),2,3)
C1 <- matrix(c(1:9),3,3)
C2 <- matrix(c(-4:3),4,2)
C3 <- matrix(c(-1:4),2,3)
C <- list(C1,C2,C3)
A <- list(A11,A12,A13,A21,A22,A23,A31,A32,A33)
B <- list(B11,B12,B13,B21,B22,B23,B31,B32,B33)

## Step 1
X1 <- matrix(0,dim(A11)[2],dim(B11)[1])
X2 <- matrix(0,dim(A12)[2],dim(B12)[1])
X3 <- matrix(0,dim(A13)[2],dim(B13)[1])
zeros_X <- list(X1,X2,X3)
X <- list(X1,X2,X3)
## Return Result
X <- general_CG(A,B,C)
## Verification
A11%*%X[[1]]%*%B11 + A12%*%X[[2]]%*%B12 + A13%*%X[[3]]%*%B13
C1
A21%*%X[[1]]%*%B21 + A22%*%X[[2]]%*%B22 + A23%*%X[[3]]%*%B23
C2
A31%*%X[[1]]%*%B31 + A32%*%X[[2]]%*%B32 + A33%*%X[[3]]%*%B33
C3

#library("microbenchmark")
#microbenchmark(general_CG(A,B,C,max.iter = 100000),sylvester(A,B,C))
