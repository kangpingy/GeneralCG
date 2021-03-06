% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Algr2_1.R
\name{general_CG}
\alias{general_CG}
\title{general_CG}
\usage{
general_CG(A, B, C, tol = 1e-20, max.iter = 1e+05)
}
\arguments{
\item{A}{input matrixes list A with length k^2}

\item{B}{input matrixes list B with length k^2}

\item{C}{input matrices list C with length k}

\item{tol}{tolerance of residuals from 0 per X element, when res < tol, return current X as the numeric root}

\item{max.iter}{the maxium iteration before ending the for loop, with default 100000}
}
\value{
solve of matrix function
}
\description{
Conduct general Conjugate Gradient Algorithm to find numeric roots of coupled matrix equation sets AXB = C
}
\examples{
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
X <- general_CG(A,B,C)
#Verification
A11\%*\%X[[1]][[1]]\%*\%B11 +
A12\%*\%X[[1]][[2]]\%*\%B12 + A13\%*\%X[[1]][[3]]\%*\%B13
C1
A21\%*\%X[[1]][[1]]\%*\%B21 +
A22\%*\%X[[1]][[2]]\%*\%B22 + A23\%*\%X[[1]][[3]]\%*\%B23
C2
A31\%*\%X[[1]][[1]]\%*\%B31 +
A32\%*\%X[[1]][[2]]\%*\%B32 + A33\%*\%X[[1]][[3]]\%*\%B33
C3

}
