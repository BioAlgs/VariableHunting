

##Simulation Setting 2

m = 50

q = 20
p = 500


# load('/data/xx/project/variablehunting/SNP.Rdata')
dat = read.table('/data/xx/project/variablehunting/SNP_info.csv', sep=',',header=T)

X = matrix(0, ncol(dat),2000)
for(i in 1:2000){
	idx = sample(1:nrow(dat),1)
	tmp = dat[idx,]
	while(sum(is.na(tmp))>0 ||  sum(tmp<0,na.rm=T)<.1*length(tmp) || sum(tmp==1)<.1*length(tmp) ){
		if(sum(tmp<0,na.rm=T)>.5*length(tmp) && sum(tmp>0, na.rm=T)>.5*length(tmp)) {
			next
		}
		idx = sample(1:nrow(dat),1)
		tmp = dat[idx,]
	}
	aa = numeric(length(tmp))
	if(sum(tmp<0,na.rm=T) > sum(tmp>0,na.rm=T)){
		aa[which(tmp>0)]=1
	}else{
		aa[which(tmp<0)]=1
	}
	X[,i] = aa
	cat(i,sum(aa),'\n')
}
X = X[, sample(which(apply(X,2,sum)<.4*n),500)]
save(X, file="singlecellsim.RData")
library(R.matlab)
writeMat("singlecellsim.mat", X=X)

n = dim(X)[1]

nMissVec= c(0,4,8,12,16)
# probVec = c(0.1,0.2,0.3)
errVec = c(1, 2, 5)/sqrt(m)

FPMat = array(0, dim =c(length(nMissVec), length(errVec), 3, 3))
FNMat = array(0, dim =c(length(nMissVec), length(errVec), 3, 3))


for(tt in 1:3){
	for(ii in 1:length(nMissVec)){
		for(jj in 1:length(errVec)){


			# X = matrix(rbinom(n*p ,1, 0.2),n,p)
			beta = matrix(rnorm(q*m,0,1),q,m)


			XT = X[,1:q]
			sigma = errVec[jj]
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


apply(FPMat,c(1,2,3),mean)
apply(FNMat,c(1,2,3),mean)
