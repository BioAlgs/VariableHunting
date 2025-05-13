##
xmat = matrix(c(-1,1,1 , -1,-1,1) , 3,2);

amat = matrix(c(0.01, 2, 0.01, 3),2,2);
ymat = xmat %*% amat;

solve(t(xmat)%*%xmat)%*%t(xmat)%*%ymat


solve(t(xmat[,i])%*%xmat[,i])%*%t(xmat[,i])%*%ymat


xmat = matrix(c(1,-1,1,1,-1,1,1,1,-1,-1,1,1),4,3);
amat = matrix(c(1, 0, 3, 2,0,2, 3,0,1),3,3)
ymat = xmat %*% amat;

solve(t(xmat)%*%xmat)%*%t(xmat)%*%ymat
i= 3

beta1 = solve(t(xmat[,-i])%*%xmat[,-i])%*%t(xmat[,-i])%*%ymat
ymat - xmat[,-i] %*% beta1


k=2
p = m = 10
n = 100
mu = 0.05
dmat = mu * rep(1,m)%*%t(rep(1,m)) + (1-mu) * diag(rep(1,m))
amat = matrix(0,n,p)
for(i in 1:n){
	amat[i,sample(p,k)] = 1
}
ymat = amat %*% dmat

cal_ang = function(x,y){
	sum(x * y)/sqrt(sum(x^2))/sqrt(sum(y^2))
}

cal_ang(dmat[,1], dmat[,2])

n=200
p = 1000
amat = matrix(rbinom(n*p, 1,0.2), n,p)
ymat = amat[,1:m] %*% dmat

cmat = matrix(0, n, n)
for(i in 1:n){
	for(j in 1:n){
		cmat[i,j] = cal_ang(ymat[i,],ymat[j,])
	}
}
diag(cmat) = 0

dhat = matrix(0,10,10)
for(d in 1:10){
	dhat[,d] = svd(ymat[which(amat[,d]!=0),])$v[,1]
}

yind = which(cmat>0.5, arr.ind = T)




svd(ymat[c(1,12,21,35,37,49,57,70,84,97),])$v[,1]
svd(ymat[c(1,12,21,35,37,49,57,70,84,97),])$u[,1]

svd(ymat[which(amat[,2]!=0),])$v[,1]

i =10
a1 = amat[yind[i,1],]
a2 = amat[yind[i,2],]

library(igraph)
g<-graph.empty(n=200, directed=TRUE)
edges <- t(yind)
g<-graph(edges, n=max(edges), directed=TRUE)

neighbors(g, V(g)[2])

plot(g, layout = layout.fruchterman.reingold,vertex.label=V(g)$number, edge.arrow.size=0.5)

a = cliques(g, min=4,max=4)

n=100
a =numeric(1000)
for(i in 1:1000){
	e1 = rnorm(n,1, )
	e2 = rnorm(n)
	a[i] = sum(e1*e2)
}


f1 = function(x){
	# return((x-0.5)^3)
	return(sin(2*pi*x)*(1-x))
}

curve(f1, 0,1)

