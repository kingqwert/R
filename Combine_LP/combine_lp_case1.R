library(Matrix)
library(MASS)
library(ROCR)
library(tidyverse)
library(qpcR)
library(ggplot2)
library(clusterGeneration)
library(rms)


# to prepare AUC, mean vectors, Sigma matrix from logits
res.logits = function(x,y){
	fit = glm(y~x,family="binomial")
	#fit = lda(y ~ x)

	AUC = prediction(fitted(fit), y) %>%
  	performance(measure = "auc")%>%
  	.@y.values %>% as.numeric

  	return(list(
		beta=coef(fit)[-1], # coefs
		AUC.k=AUC,
  		mu1 = apply(x[y==1,],2,mean,na.rm=T),
  		mu0 = apply(x[y==0,],2,mean,na.rm=T),
  		mu=apply(x[y==1,],2,mean,na.rm=T)-apply(x[y==0,],2,mean,na.rm=T),
  		sd1 = apply(x[y==1,],2,sd,na.rm=T),
  		sd0 = apply(x[y==0,],2,sd,na.rm=T)
  		))
}


make.studies = function(N1,N0,m1,m0,Sigma1,Sigma0,cutoff1){	
 	# simulate design matrix
	X.g1 = mvrnorm(N1,mu=m1,Sigma=Sigma1)
	X.g0 = mvrnorm(N0,mu=m0,Sigma=Sigma0)
	
	
	tmp = data.frame(
		y=c(rep(0,N0),rep(1,N1)),
		x=rbind(X.g0,X.g1))
	
	tmp$x.4 = ifelse(tmp$x.4>cutoff1,1,0)
	#tmp$x.5 = ifelse(tmp$x.5>cutoff2,1,0)
	return(tmp)
}

#######################################
# Prepare for calculating AUC using our method and conventional methods
#######################################
AUC_total_our = matrix(NA,ncol=5,nrow=1000)
AUC_total_mean = matrix(NA,ncol=5,nrow=1000)
AUC_total_study1 =matrix(NA,ncol=5,nrow=1000)
AUC_total_study2 =matrix(NA,ncol=5,nrow=1000)
AUC_total_all = matrix(NA,ncol=5,nrow=1000)
WhichAUC = matrix(NA,ncol=5,nrow=1000)

# simulation parameters
cors = c(-0.3,-0.1,0.1,0.3,0.6) 
N1=N0=500

for(c in 1:5){
	for(sim in 1:1000){
		studies =list()
		data.list =list()
		Ck=list()
		
		mu.g1 = c(rep(1,5),rep(0.5,5))
		mu.g0 = rep(0,10)
		
		# correlation matrix
		cor=0.5
		cor_studies = cors[c]
		R = matrix(rep(cor,10^2), ncol=10) 
		# R12
		R[c(1:3),c(4:6)]=cor_studies
		R[c(4:6),c(1:3)]=cor_studies
		# R13
		R[c(1:3),c(7:10)]=cor_studies
		R[c(7:10),c(1:3)]=cor_studies
		# R23
		R[c(4:6),c(7:10)]=cor_studies
		R[c(7:10),c(4:6)]=cor_studies
		diag(R) = 1 

		
		# prepare Sigma
		# Sigma1 = diag(rep(2,10)) %*% R %*% diag(rep(2,10))
		Sigma1 = diag(c(1:10)) %*% R %*% diag(c(1:10))

		
		# simulate design matrix
		X.g1 = mvrnorm(N1,mu=mu.g1,Sigma=Sigma1)
		X.g0 = mvrnorm(N0,mu=mu.g0,Sigma=Sigma1)
		
		test=data.frame(
			y=c(rep(0,N0),rep(1,N1)),
			x=rbind(X.g0,X.g1))
		
		test$x.4 = ifelse(test$x.4>1,1,0)
		#test$x.5 = ifelse(test$x.5>1,1,0)
		test$x.9 = ifelse(test$x.9>1,1,0)
		#test$x.10 = ifelse(test$x.10>1,1,0)

		Sigma_test = var(test[,-1])
		
	
		# number of variables in each study
		K=2
		
		for(i in 1:K){
			if(i==1){
				data.list[[i]]=make.studies(N1=500,N0=500,m1=mu.g1[1:5],m0=mu.g0[1:5], Sigma1=Sigma1[c(1:5),c(1:5)], Sigma0=Sigma1[c(1:5),c(1:5)],cutoff1=1)
			}else if(i==2){
				data.list[[i]]=make.studies(N1=500,N0=500,m1=mu.g1[6:10],m0=mu.g0[6:10], Sigma1=Sigma1[c(6:10),c(6:10)], Sigma0=Sigma1[c(6:10),c(6:10)],cutoff1=1)
			}
			
			studies[[i]] = res.logits(as.matrix(data.list[[i]][,-1]),data.list[[i]]$y)
		}


		
		# calculate w_opt
		B0hat = cbind(
			c(studies[[1]]$beta,rep(0,5)),
			c(rep(0,5),studies[[2]]$beta)
		)
		muhat = cbind(
			c(studies[[1]]$mu,rep(NA,5)),
			c(rep(NA,5),studies[[2]]$mu)
		)

		delta_hat = apply(muhat,1,mean,na.rm=T)
		
		
		w_opt2 = c(solve(t(B0hat)%*%Sigma_test%*%B0hat)%*%t(B0hat)%*%delta_hat)
		
		final_opt2=list()
		for(i in 1:K){
			betas = studies[[i]]$beta*w_opt2[i]
			if(i==1){
				final_opt2[[i]] = c(betas,rep(0,5))	
			}else if(i==2){
				final_opt2[[i]] = c(rep(0,5),betas)	
			}
		}
	

		###############
		# Comparison
		###############
		# let's compare
		#our model
		score=list()
		for(i in 1:K){
			score[[i]] = as.matrix(test[,-1]) %*%final_opt2[[i]]
		}
		
		pred.our = Reduce('+',score) 
		# phat_our = 1/(1+exp(-pred.our))
		# val.prob(phat_our, test$y, cex=.5)
		AUC_total_our[sim,c]=prediction(pred.our, test$y) %>%
		  	performance(measure = "auc")%>%
		  	.@y.values %>% as.numeric
		
		
		#comparison groyup 1: each models
		score.each=list()
		phat_each=list()
		AUC.each = numeric(K)
		for(i in 1:K){
			if(i==1){
				beta = c(studies[[i]]$beta,rep(0,5))	
			}else if(i==2){
				beta = c(rep(0,5),studies[[i]]$beta)	
			}
			
			score.each[[i]] =as.matrix(test[,-1]) %*%beta
			
			# phat_each[[i]] = 1/(1+exp(-score.each[[i]]))
			# val.prob(phat_each[[i]], test$y, cex=.5)
			
			AUC.each[i]=prediction(score.each[[i]], test$y) %>%
		  		performance(measure = "auc")%>%
		  		.@y.values %>% as.numeric
		}
	
		AUC_total_study1[sim,c]=AUC.each[1]
		AUC_total_study2[sim,c]=AUC.each[2]

		# comparison group 2: just taking mean of f_k
		B0hat2 =  cbind(
			c(studies[[1]]$beta,rep(NA,5)),
			c(rep(NA,5),studies[[2]]$beta)
		)
		mean.coef = apply(B0hat2,1,mean,na.rm=T)
		mean.score = as.matrix(test[,-1]) %*%mean.coef

		
		AUC_total_mean[sim,c]=prediction(mean.score, test$y) %>%
		  		performance(measure = "auc")%>%
		  		.@y.values %>% as.numeric

		WhichAUC_tmp = which.max(sapply(1:2,function(x){studies[[x]]$AUC.k}))
		if(WhichAUC_tmp==1){
			WhichAUC[sim,c] = AUC.each[1]
		}else if(WhichAUC_tmp==2){
			WhichAUC[sim,c] = AUC.each[2]
		}

	}
}



AUC_our = AUC_total_our %>% as.data.frame %>% gather(key=correlation,value=AUC,1:5)
AUC_mean = AUC_total_mean %>% as.data.frame %>% gather(key=correlation,value=AUC,1:5)
AUC_study1 = AUC_total_study1%>% as.data.frame %>% gather(key=correlation,value=AUC,1:5)
AUC_study2 = AUC_total_study2 %>% as.data.frame %>% gather(key=correlation,value=AUC,1:5)
AUC_prior = WhichAUC %>% as.data.frame %>% gather(key=correlation,value=AUC,1:5)

df = cbind(
	method=factor(c(rep("Our",1000*5),rep("Mean",1000*5),rep("Study 1",1000*5),rep("Study 2",1000*5),rep("AUC prior",1000*5)),levels=c("Our","Mean","Study 1","Study 2","AUC prior")),
	rbind(AUC_our,AUC_mean,AUC_study1,AUC_study2,AUC_prior)
	)
df$correlation=factor(df$correlation,levels=c("V1","V2","V3","V4","V5"),labels=c("-0.3","-0.1","0.1","0.3","0.6"))
g <- ggplot(df, aes(x = correlation, y = AUC, fill=method))
g <- g + geom_boxplot()+theme_bw()+ggtitle("")
g

