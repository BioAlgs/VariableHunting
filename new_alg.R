# k=2
# p = m = 10
# n = 100
# mu = 0.05
# dmat = mu * rep(1,m)%*%t(rep(1,m)) + (1-mu) * diag(rep(1,m))
# amat = matrix(0,n,p)
# for(i in 1:n){
# 	amat[i,sample(p,k)] = 1
# }
# ymat = amat %*% dmat + matrix(rnorm(n*m, 0, 0.1), n,m)

cal_ang = function(x,y){
	sum(x * y)/sqrt(sum(x^2))/sqrt(sum(y^2))
}
varHunting = function(ymat, xmat, thred = 0.2){
# ymat =Y
# xmat = XO
	n = nrow(ymat)
	m = ncol(ymat)
	cmat = matrix(0, n, n)
	for(i in 1:n){
		for(j in 1:n){
			cmat[i,j] = cal_ang(ymat[i,],ymat[j,])
		}
	}
	
	diag(cmat) = 0
	yind = which(cmat>thred, arr.ind = T)
	library(igraph)
	g<-graph.empty(n, directed=TRUE)
	edges <- t(yind)
	g<-graph(edges, n=max(edges), directed=TRUE)

###estimate the sparsity
	pSparsity = (min(0.5, nrow(yind)/(n^2)) + 0.15)/3.5
	nNeighbor = max(5, qbinom(0.0001, n, pSparsity))

	
###Select the candidate clique
	niter = 5000
	cliqueMat = vector("list",0)
	count= 1
	for(t in 1:niter){
		idx = sample(1:nrow(yind),1)
		vex = numeric(3)
		vex[1:2] = yind[idx,]
		if(length(intersect(neighbors(g, V(g)[vex[1]]), neighbors(g, V(g)[vex[2]])))==0){
			next
		}
		vex[3] = sample(intersect(neighbors(g, V(g)[vex[1]]), neighbors(g, V(g)[vex[2]])),1)

		n1 = intersect(neighbors(g, vex[1]), neighbors(g, vex[2]))
		n2 = intersect(n1, neighbors(g, vex[3]))
		
		cn = length(n2)
		if(cn >= nNeighbor){
			cliqueMat[[count]] =  union(n2, vex)
			count = count + 1
			# cat(vex,'\n')
		}
	}



###Searching for maximum clique
	maxCliqueMat = vector("list",0)
	idMaxClique = 1
	id = numeric(length(cliqueMat))
	for(s in 1:length(cliqueMat)){
		if(s==1){
			maxCliqueMat[[idMaxClique]] = cliqueMat[[s]]
			idMaxClique =  idMaxClique+1
			# id[s] = which.max(apply(amat[tmp,],2,sum))
		}else{
			tmp = cliqueMat[[s]]
			# id[s] = which.max(apply(amat[tmp,],2,sum))
			idOverlap = 0
			for(i in 1:length(maxCliqueMat)){
				if(length(intersect(tmp, maxCliqueMat[[i]]))>0.85* length(tmp)){
					# maxCliqueMat[[i]] = intersect(tmp, maxCliqueMat[[i]])
					idOverlap = idOverlap+1
				}
			}
			if(idOverlap == 0){
				maxCliqueMat[[idMaxClique]] = tmp
				idMaxClique = idMaxClique+1
				# cat(tmp,'\n')
			}
		}
	}


###Variable Hunting

	cat('Variable hunting ......\n')
	varId = numeric(length(maxCliqueMat))
	Xhat=NULL
	for(i in 1:length(maxCliqueMat)){
		if(i==1){
			tmp = rep(0, n)
			tmp[maxCliqueMat[[1]]] = 1
			Xhat = tmp
		}else{
			tmp = rep(0, n)
			tmp[maxCliqueMat[[i]]] = 1
			Xhat = cbind(Xhat,tmp)
		}
		pval = numeric(ncol(xmat))
		for(j in 1:ncol(xmat)){
			test = chisq.test(table(tmp,xmat[,j]))
			pval[j] = test$p.value
		}
		qval = p.adjust(pval,method='bonferroni')
		if(min(qval)<0.01){
			varId[i] = which.min(pval)
		}else{
			varId[i] = NA
		}
	}

return(varId)

}






which.max(cor(Xhat[,1],X))
chisq.test(table(Xhat[,1],X[,16]))

XT[maxCliqueMat[[1]],]
XT[maxCliqueMat[[2]],]
XT[maxCliqueMat[[3]],]
XT[maxCliqueMat[[4]],]
XT[maxCliqueMat[[5]],]
XT[maxCliqueMat[[6]],]
XT[maxCliqueMat[[7]],]
XT[maxCliqueMat[[8]],]


XT[maxCliqueMat[[which(idCluster==18)[1]]],]
XT[maxCliqueMat[[which(idCluster==18)[2]]],]

idCluster = numeric(length(maxCliqueMat))
for(i in 1:length(maxCliqueMat)){
	tmp = maxCliqueMat[[i]]
	idCluster[i] = which.max(apply(XT[tmp,],2,sum))
	cat(which.max(apply(XT[tmp,],2,sum)),'\n')
}
unique(idCluster)
length(unique(idCluster))


XT[cliqueMat[[1]],]


idCluster = numeric(length(cliqueMat))
for(i in 1:length(cliqueMat)){
	tmp = cliqueMat[[i]]
	idCluster[i] = which.max(apply(XT[tmp,],2,sum))
	# cat(which.max(apply(XT[tmp,],2,sum)),'\n')
}
unique(idCluster)
length(unique(idCluster))

apply(XT[,-unique(idCluster)],2,sum)

amat[maxCliqueMat[[7]],]
amat[maxCliqueMat[[1]],]


n * 2 / (10*m)




ymat1 = ymat[which(XT[,18]!=0),]
	n = nrow(ymat1)
	m = ncol(ymat1)
	cmat = matrix(0, n, n)
	for(i in 1:n){
		for(j in 1:n){
			cmat[i,j] = cal_ang(ymat1[i,],ymat1[j,])
		}
	}




