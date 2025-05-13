

##Simulation Setting 1
n = 200
m = 50
q = 20
p = 500

nMissVec= c(0,4,8,12,16)
probVec = c(0.1,0.2,0.3)

FPMat = array(0, dim =c(length(nMissVec), length(probVec), 3, 30))
FNMat = array(0, dim =c(length(nMissVec), length(probVec), 3, 30))


for(tt in 1:30){
	for(ii in 1:length(nMissVec)){
		for(jj in 1:length(probVec)){


			X = matrix(rbinom(n*p ,1, probVec[jj]),n,p)
			beta = matrix(rnorm(q*m,0,1),q,m)


			XT = X[,1:q]
			sigma = 0.1
			Y = XT%*%beta + matrix(sigma * rnorm(n*m,0,1),n,m)

			nMiss = nMissVec[ii]
			idxOT = sample(1:q, size=q-nMiss, replace=F)
			XO = cbind(XT[,idxOT], X[,(q+1):p])
			betaO = rbind(beta[idxOT,], matrix(0, p-q, m))


			###Using VH
			fit1 = varHunting(Y, XO)
			FNMat[ii,jj,1,tt] = ((q-nMiss) - length(which(unique(fit1) %in% (1:(q-nMiss))))) * m
			FPMat[ii,jj,1,tt] = length(which(unique(fit1) %in% ((q-nMiss+1):(p-nMiss)))) * m



			###Using Lasso
			library(glmnet)
			# install.packages('glmnet', repo='https://cran.rstudio.com')

			coefList = vector("list",m)
			coefMat = matrix(0, p-nMiss ,m)
			for(k in 1:m){
				fit = cv.glmnet(XO, Y[,k] , family="gaussian", nfolds=5)
				tmp = coef(fit, s=fit$lambda.min)[-1,1]
				coefList[[k]] = which(tmp != 0)
				coefMat[,k] = tmp
			}

			FNMat[ii,jj,2,tt] =  length(which(coefMat[1:(q-nMiss), ] == 0 ) )
			FPMat[ii,jj,2,tt] =  length(which(coefMat[(q-nMiss+1):(n-nMiss), ] != 0 )) 


			###Using remMap
			library(remMap)

			lamL1.v=exp(seq(log(10),log(20), length=3)) 
			lamL2.v=seq(0,5, length=3)
			try1 =remMap.CV(X=XO, Y=Y,lamL1.v, lamL2.v, C.m=NULL, fold=5, seed=1)
			pick=which.min(as.vector(try1$rss.cv))	
			lamL1.pick=try1$l.index[1,pick]
			lamL2.pick=try1$l.index[2,pick]
			result=remMap(XO, Y,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
			FPMat[ii,jj,3,tt]=sum(result$phi!=0 & betaO==0)
			FNMat[ii,jj,3,tt]=sum(result$phi==0 & betaO!=0)

	cat(ii,jj,tt,FPMat[ii,jj,,tt],'\n')
	cat(ii,jj,tt,FNMat[ii,jj,,tt],'\n')

		}		
	}
}


##Simulation Setting 2
n = 200
m = 50
q = 20
p = 500

library(bindata)

## Construct a binary correlation matrix

XT = matrix(0,n,q)
XR = NULL
for(i in 1:q){
	rho <- 0.5
	mcov = matrix(0.5, p/q,p/q)
	diag(m) = 1
	tmpX = rmvbin(n, margprob = rep(0.1, p/q) , bincorr = mcov) 
	tmpIdx = sample(1:(p/q),1)
	XT[,i] = tmpX[,tmpIdx]
	XR = cbind(XR, tmpX[, -tmpIdx])
}

X = cbind(XT, XR)
beta = matrix(rnorm(q*m,0,1),q,m)

sigma = .1
Y = XT%*%beta + matrix(sigma * rnorm(n*m,0,1),n,m)



library(glmnet)
# install.packages('glmnet', repo='https://cran.rstudio.com')

nMiss = 2
idxOT = sample(1:q, size=q-nMiss, replace=F)

XO = cbind(XT[,idxOT], X[,(q+1):p])

coefList = vector("list",m)
for(j in 1:m){
	fit = cv.glmnet(XO, Y[,j] , family="gaussian")
	coefList[[j]] = which(coef(fit, s=fit$lambda.min) != 0)
}


unique(unlist(coefList)) 


