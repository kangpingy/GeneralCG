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
