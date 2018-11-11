
library("limSolve") 

b = c(1,2,3)

A = matrix(nrow = 3, ncol = 3)
A[,1] = c(0.25,0.50,3/4)
A[,2] = c(0.1,2/10,3/10)
A[,3] = c(1,2,3)

E = rep(1, length(b))
f = 1


lsei(A, b, E, f)
