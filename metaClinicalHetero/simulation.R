library(MASS)
library(Matrix)
library(nlme)
library(mvmeta)
library(ggplot2)


#Set basic statics 
#latent data generating process
set.seed(12345)

#mean
bmu = 0.5

#Make sample studies
MakeSample = function(K, bvar){
	avg = rnorm(K, mean = bmu, sd = sqrt(bvar))
	var = 0.25*rchisq(K,1)
	se = ifelse(var<0.009,sqrt(0.009),
		ifelse(var>0.6,sqrt(0.6),sqrt(var)))

	return(list(avg,se))
}

# revised reml
re.reml = function(y,sigma.i,a=1000,b=0.1,par){
	sigmaREREML =par[1]
	beta = par[2]

	0.5*(sum(log(sigmaREREML^2+sigma.i)+(y-beta)^2/(sigmaREREML^2+sigma.i))+
		log(sum((sigmaREREML^2+sigma.i)^(-1))))-b*exp(-length(y)^2/a)*log(sigmaREREML^2)
}


#Set the matrix to put the results
res.beta = list()
res.psi = list()
res.I = list()

reml.beta = matrix(NA, ncol=6, nrow =1000)
ml.beta = matrix(NA, ncol=6, nrow =1000)
dl.beta = matrix(NA, ncol=6, nrow =1000)
propose.beta = matrix(NA, ncol=6, nrow =1000)

reml.psi = matrix(NA, ncol=6, nrow =1000)
ml.psi = matrix(NA, ncol=6, nrow =1000)
dl.psi = matrix(NA, ncol=6, nrow =1000)
propose.psi = matrix(NA, ncol=6, nrow =1000)

reml.I = matrix(NA, ncol=6, nrow =1000)
ml.I = matrix(NA, ncol=6, nrow =1000)
dl.I = matrix(NA, ncol=6, nrow =1000)
propose.I = matrix(NA, ncol=6, nrow =1000)

# main iterations
for(k in 1:4){
	for (j in 1:6){
		for (i in 1:1000){
			coefs<-MakeSample(K = c(5, 10, 30, 50)[k], bvar = c(0, 0.01, 0.02, 0.05, 0.1, 0.2)[j], 1000, 100)
			
			if(is.na(coefs[1])==FALSE){
				
				reml = mvmeta(formula=coefs[[1]]~1,S=coefs[[2]],method="reml")
				ml = mvmeta(formula=coefs[[1]]~1,S=coefs[[2]],method="ml")
				dl = mvmeta(formula=coefs[[1]]~1,S=coefs[[2]],method="mm")
				
				reml.psi[i,j] = reml$Psi
				ml.psi[i,j] = ml$Psi
				dl.psi[i,j] = dl$Psi
			
				reml.beta[i,j] = coef(reml)
				ml.beta[i,j] = coef(ml)
				dl.beta[i,j] = coef(dl)
				
				reml.I[i,j] = reml$Psi/(reml$Psi+sum(coefs[[2]]^(-1)*(c(5, 10, 30, 50)[k]-1))/(sum((coefs[[2]])^(-1))^2-sum((coefs[[2]]^(-1))^2)))
				ml.I[i,j] = ml$Psi/(ml$Psi+sum(coefs[[2]]^(-1)*(c(5, 10, 30, 50)[k]-1))/(sum((coefs[[2]])^(-1))^2-sum((coefs[[2]]^(-1))^2)))
				dl.I[i,j] = dl$Psi/(dl$Psi+sum(coefs[[2]]^(-1)*(c(5, 10, 30, 50)[k]-1))/(sum((coefs[[2]])^(-1))^2-sum((coefs[[2]]^(-1))^2)))
				
				propose = try(optim(par=c(0.1,0.1),fn=function(r){re.reml(y=coefs[[1]],sigma.i=coefs[[2]],par=r)},method="BFGS"),silent=T)

				if(inherits(propose, "try-error") | propose$convergence ==1 | propose$par[1] > 1 | propose$par[2] > 3){
					propose.psi[i,j] = NA
					propose.beta[i,j] = NA
					propose.I[i,j] = NA
				}

				propose.psi[i,j] = propose$par[1]^2
				propose.beta[i,j] = propose$par[2]
				propose.I[i,j] = propose$par[1]^2/(propose$par[1]^2+sum(coefs[[2]]^(-1)*(c(5, 10, 30, 50)[k]-1))/(sum((coefs[[2]])^(-1))^2-sum((coefs[[2]]^(-1))^2)))
			}else{
		
				i<-i
			}
				
		}
	}
	res.beta[[k]] = list(reml.beta, ml.beta, dl.beta, propose.beta)
	res.psi[[k]] = list(reml.psi, ml.psi, dl.psi, propose.psi)
	res.I[[k]] = list(reml.I, ml.I, dl.I, propose.I)

}



