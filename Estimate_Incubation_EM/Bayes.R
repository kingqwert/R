library(dplyr)
library(ggplot2)
library(EnvStats)
library(MASS)
library(matrixcalc)
library(Matrix)
library(foreach)
library(doParallel)
library(brms)
library(rstan)
totalCores = detectCores()
rstan_options(auto_write = TRUE)


res = list()
res_tmp = list()
fits_tmp = list()
fits = list()
result = list()
fit_model =stan_model(file='bayes_version2.stan')


N = *** 
sample_n = *** # 
sample_m = *** # 
sample_l = *** # 
meanlog = *** # 
sdlog = *** # 
threshold = *** # 

for(k in 1:1000){
		show(c(i,j,k))
		y_true=rlnorm3(N,meanlog=meanlog,sdlog=sdlog,threshold=threshold)
		
		if((mean(y_true)+3*sd(y_true))<max(y_true)){
			while((mean(y_true)+3*sd(y_true))<max(y_true)){
				y_true=rlnorm3(N,meanlog=meanlog,sdlog=sdlog,threshold=threshold)
			}
		}

		y = y_true[1:sample_n]
		y_LU = cbind( 
			runif(sample_m, min = threshold,max = y_true[(sample_n+1):(sample_n+sample_m)]),
			runif(sample_m, min = y_true[(sample_n+1):(sample_n+sample_m)],max = y_true[(sample_n+1):(sample_n+sample_m)]+1) 
			)
		# y_C = y_true[(sample_n+sample_m+1):N]-runif(sample_l,0,sdlog)
		y_C = runif(sample_l, min = threshold, max = y_true[(sample_n+sample_m+1):N])
		
		y_tmp = c(y,y_LU[,1],y_C)
		
		gamma0 = min(min(y),min(y_LU),min(y_C))-0.1
		mu0 = (sum(log(y-gamma0))+sum(log(y_LU[,1]-gamma0))+sum(log(y_C-gamma0)))/length(y_true)
		sigma0 = sqrt((sum((log(y-gamma0)-mu0)^2)+sum((log(y_LU[,1]-gamma0)-mu0)^2)+sum((log(y_C-gamma0)-mu0)^2))/length(y_true))
		init_values <- list(
		  list(mu = 1.0, sigma = 0.5, shift = gamma0-1),  # chain 1
		  list(mu = 1.1, sigma = 0.6, shift = gamma0-2),  # chain 2
		  list(mu = mu0, sigma = sigma0, shift = gamma0-0.1)  # chain 3
		)
		d = data.frame(
			rbind(
				cbind(y,y,rep('none',length(y))),
				cbind(y_LU,rep('interval',dim(y_LU)[1])),
				cbind(y_C,y_true[(sample_n+sample_m+1):N],rep('right',length(y_C)))
			)
		)
		names(d) = c('y','y1','cen')
		d$y = as.numeric(d$y)
		d$y1 = as.numeric(d$y1)
		cens = ifelse(d$cen=='none',0,ifelse(d$cen == 'right',-1,1))
		data = list(
			N = dim(d)[1], 
			y=d$y, 
			censoring_limits = d$y, 
			upper_limits = d$y1, 
			lower_limits = d$y, 
			cens = cens, 
			min_y = gamma0,
			theta = mu0,
			theta1 = sigma0,
			tau = 10,
			tau1 = 10)
		fit = sampling(fit_model,data=data,iter=5000,chains = 3,cores=1, init=init_values, open_progress=FALSE,show_messages=FALSE,control=list(max_treedepth=16))
		fits_tmp[[k]] = fit
		res_tmp[[k]] = cbind(summary(fit)$summary[,1],summary(fit)$summary[,4],summary(fit)$summary[,8])[1:3,]

	}

	# res[[j]] = res_tmp
	# fits[[j]] = fits
	res = res_tmp
	fits = fits_tmp
# }
result = list(res,fits)







