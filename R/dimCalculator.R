#'dimCalculator
#'
#'Give the list of input A,B and C check if the input is valid in dimension for martix group X, and return the
#'dim of X if it is valid
#'
#'@param A input list consists of n^2 coupled matrixes A
#'
#'@param B input list consists of n^2 coupled matrixes B
#'
#'@param C input list consists of n coupled matrixes B
#'
#'@return return the dimension of X if the input is valid.
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
#'dimCalculator(A,B,C)
#'
#'@export
#'
dimCalculator <- function(A,B,C){
  len <- length(C)
  dim_X <- matrix(0,2,len)
  for (i in 1: len^2){
    mod_X <- (i-1)%%len+1
    group_X <- (i-1)%/%len+1
    if (dim(A[[i]])[1] != dim(C[[group_X]])[1]){
      stop("matrix A do not fit with matrix C")
    }
    if (dim(B[[i]])[2] != dim(C[[group_X]])[2]){
      stop("matrix B do not fit with matrix C")
    }
    if (i <= len){
      dim_X[,mod_X] <- c(dim(A[[i]])[2],dim(B[[i]])[1])
    }else{
      if (!all.equal(dim_X[,mod_X],c(dim(A[[i]])[2],dim(B[[i]])[1]))){
        stop("Matrix pair A1 and B1 do not fit with Matrix pair A",i, " and B",i)
      }
    }
  }
  return(dim_X)
}
