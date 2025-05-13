

pvec = numeric(20)
count = 1
for(t in seq(.1,.3, length.out=20)){
	X = matrix(rbinom(n*p ,1, t),n,p)

	XT = X[,1:q]

	eMat = matrix(0, n,n)
	for(i in 1: n){
		for(j in 1:n){
			eMat[i,j] = sum(XT[i,] * XT[j,])
		}
	}

	pvec[count] =length(which(eMat>0))/n^2
	count= count+1
}

pdf('1.pdf')

plot(seq(.1,.3, length.out=20), pvec)

dev.off()

lm(pvec ~ seq(.1,.3, length.out=20))