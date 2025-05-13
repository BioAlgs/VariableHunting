
n = 20 
x = rnorm(n, 0, 1/sqrt(20))
y = rnorm(n, 0, 1/sqrt(20))

sum(x*y)



X = matrix(c(1,0,1,1,0,0,1,1,0,1,1,1), 4,3)

y = X %*% c(1,3,0)


Xo = X[,c(1,3)]

solve(t(Xo) %*% Xo) %*% t(Xo) %*% y



X = matrix(c(1,0,1,1,1,0,0,1,1), 3,3)
y = X %*% c(1/2,1/2,0)
Xo = X[,c(1,3)]

solve(t(Xo) %*% Xo) %*% t(Xo) %*% y

