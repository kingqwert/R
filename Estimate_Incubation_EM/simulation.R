library(dplyr)
library(ggplot2)
library(EnvStats)
library(MASS)
library(matrixcalc)
library(Matrix)

psi = function(x, mu, sigma, gamma){
	return((log(x-gamma)-mu)/sigma)
}

phi = function(x){
	return(exp(-x^2/2)/sqrt(2*pi))
}

Cphi = function(x){
	return(pnorm(x))
}

ll = function(y,y_LU,y_C,theta){
	mu = theta[1]
	sigma = theta[2]
	gamma = theta[3]
	
	n = length(y)
	m = dim(y_LU)[1]
	l = length(y_C)
	
	y_L = y_LU[,1]
	y_U = y_LU[,2]

	return(
		-n*log(sigma) - sum(log(y-gamma))-0.5*sum((psi(y,mu,sigma,gamma))^2)+
		sum(log(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma))))+
		sum(log(1-Cphi(psi(y_C,mu,sigma,gamma))))
	)
}

mu_t1 = function(y,y_LU,y_C,theta_t0){
	mu = theta_t0[1]
	sigma = theta_t0[2]
	gamma = theta_t0[3]

	y_L = y_LU[,1]
	y_U = y_LU[,2]

	n = length(y)
	m = dim(y_LU)[1]
	l = length(y_C)

	return(
		(sum(log(y-gamma))+
		sum(sigma*(phi(psi(y_L,mu,sigma,gamma))-phi(psi(y_U,mu,sigma,gamma)))/(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma))))+m*mu+
		sum(sigma*phi(psi(y_C,mu,sigma,gamma))/(1-Cphi(psi(y_C,mu,sigma,gamma))))+l*mu)/(n+m+l)
	)
}


sigma_t1_1 = function(y,y_LU,y_C,theta_t0){
	mu = theta_t0[1]
	sigma = theta_t0[2]
	gamma = theta_t0[3]

	y_L = y_LU[,1]
	y_U = y_LU[,2]

	n = length(y)
	m = dim(y_LU)[1]
	l = length(y_C)

	return(
		sqrt((sum((log(y-gamma)-mu)^2)+
		m*sigma^2+sum(sigma^2*(psi(y_L,mu,sigma,gamma)*phi(psi(y_L,mu,sigma,gamma))-psi(y_U,mu,sigma,gamma)*phi(psi(y_U,mu,sigma,gamma)))/(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma))))+
		l*sigma^2+sum(sigma^2*(psi(y_C,mu,sigma,gamma)*phi(psi(y_C,mu,sigma,gamma)))/(1-Cphi(psi(y_C,mu,sigma,gamma)))))/(n+m+l))
	)
}


opt_gamma_t1_1 = function(gamma,y,y_LU,y_C,mu,sigma){
	y_L = y_LU[,1]
	y_U = y_LU[,2]

	n = length(y)
	m = dim(y_LU)[1]
	l = length(y_C)

	return(
		abs(sum(1/(y-gamma)+(log(y-gamma)-mu)/sigma^2/(y-gamma))+
		sum(exp((sigma^2-2*mu)/2)*((Cphi(psi(y_U,mu,sigma,gamma)+sigma)-Cphi(psi(y_L,mu,sigma,gamma)+sigma))/(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma)))))+
		sum(exp((sigma^2-2*mu)/2)*((phi(psi(y_L,mu,sigma,gamma)+sigma)-phi(psi(y_U,mu,sigma,gamma)+sigma))/sigma/(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma)))+(mu-sigma^2)*(Cphi(psi(y_U,mu,sigma,gamma)+sigma)-Cphi(psi(y_L,mu,sigma,gamma)+sigma))/sigma^2/(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma)))))-
		mu/sigma^2*sum(exp((sigma^2-2*mu)/2)*((Cphi(psi(y_U,mu,sigma,gamma)+sigma)-Cphi(psi(y_L,mu,sigma,gamma)+sigma))/(Cphi(psi(y_U,mu,sigma,gamma))-Cphi(psi(y_L,mu,sigma,gamma)))))+
		sum(exp((sigma^2-2*mu)/2)*(1-Cphi(psi(y_C,mu,sigma,gamma)+sigma))/(1-Cphi(psi(y_C,mu,sigma,gamma))))+
		sum(exp((sigma^2-2*mu)/2)*((phi(psi(y_C,mu,sigma,gamma)+sigma))/sigma/(1-Cphi(psi(y_C,mu,sigma,gamma)))+(mu-sigma^2)*(1-Cphi(psi(y_C,mu,sigma,gamma)+sigma))/sigma^2/(1-Cphi(psi(y_C,mu,sigma,gamma)))))-
		mu/sigma^2*sum(exp((sigma^2-2*mu)/2)*(1-Cphi(psi(y_C,mu,sigma,gamma)+sigma))/(1-Cphi(psi(y_C,mu,sigma,gamma))))
	))
}


gamma_t1 = function(y,y_LU,y_C,mu,sigma,interval){
	return(
		optimize(opt_gamma_t1_1,y=y,y_LU=y_LU,y_C=y_C,mu=mu,sigma=sigma,interval=interval)$minimum
	)
}


epsiolon_accel_tmp = function(vec){
	return(unlist(vec)/c(t(vec)%*%vec))
}

epsiolon_accel = function(theta_tm1,theta_t,theta_tp1){
	tmp1 = epsiolon_accel_tmp(theta_tm1-theta_t)+epsiolon_accel_tmp(theta_tp1-theta_t)
	return(theta_t+epsiolon_accel_tmp(tmp1))
}


EM_imcubation_accel = function(y,y_LU,y_C,theta0,epsilon,interval,max_iter){
	mu1 = mu_t1(y,y_LU,y_C,theta_t0=theta0)
	theta0[1] = mu1
	sigma1 = sigma_t1_1(y,y_LU,y_C,theta_t0=theta0)
	gamma1 = gamma_t1(y,y_LU,y_C,mu=mu1,sigma=sigma1,interval)

	theta1 = c(mu1, sigma1, gamma1)
	delta = 1000
	ll_value = ll(y,y_LU,y_C,theta1)
	thetas = theta1
	iter = 1

	while((delta > sqrt(epsilon)) & (iter<max_iter)){
		iter = iter+1
		mu1 = mu_t1(y,y_LU,y_C,theta_t0=theta1)
		theta1[1] = mu1
		sigma1 = sigma_t1_1(y,y_LU,y_C,theta_t0=theta1)
		gamma1 = gamma_t1(y,y_LU,y_C,mu=mu1,sigma=sigma1,interval)
		
		theta1 = c(mu1, sigma1, gamma1)
		thetas = rbind(thetas,theta1)
		at_t = dim(thetas)[1]
		if(at_t==3){
			theta_dot = epsiolon_accel(theta_tm1 = thetas[1,],theta_t=thetas[2,],theta_tp1=thetas[3,])
			theta_dots = theta_dot
		}else if(at_t>3){
			theta_dot = epsiolon_accel(theta_tm1 = thetas[at_t-2,],theta_t=thetas[at_t-1,],theta_tp1=thetas[at_t,])
			theta_dots = rbind(theta_dots,theta_dot)
			delta = t(c(theta_dots[dim(theta_dots)[1],]-theta_dots[dim(theta_dots)[1]-1,]))%*%c(theta_dots[dim(theta_dots)[1],]-theta_dots[dim(theta_dots)[1]-1,])
		}
		# show(c(delta,iter))

		ll_value = c(ll_value,ll(y,y_LU,y_C,theta1))
	}
	return(list(theta1,ll_value,theta_dots))
}

Boot_Cov = function(y,y_LU,y_C,L){
	n = length(y)
	m = dim(y_LU)[1]
	l = length(y_C)
	N = n+m+l
	theta_boot = list()

	for(i in 1:L){
		id = sort(sample(1:N,replace=T))
		y_boot = y[id[id<=n]]
		
		if(length(id[id<n+m & id>n]-n)==1){
			y_LU_boot = matrix(y_LU[id[id<n+m & id>n]-n,],ncol=2)
		}else{
			y_LU_boot = y_LU[id[id<n+m & id>n]-n,]	
		}
		y_C_boot = y_C[id[id>n+m]-n-m]
		
		gamma0 = min(min(y_boot),min(y_LU_boot),min(y_C_boot))-0.0001
		mu0 = (sum(log(y_boot-gamma0))+sum(log(y_LU_boot[,1]-gamma0))+sum(log(y_C_boot-gamma0)))/N
		sigma0 = sqrt((sum((log(y_boot-gamma0)-mu0)^2)+sum((log(y_LU_boot[,1]-gamma0)-mu0)^2)+sum((log(y_C_boot-gamma0)-mu0)^2))/N)
			
		theta0 = c(mu0, sigma0, gamma0)
			
		theta_boot[[i]] = unlist(EM_imcubation_accel(y_boot,y_LU_boot,y_C_boot,theta0,epsilon=0.001,interval=c(threshold*2,gamma0),max_iter=10000)[1])
	}
	boot = do.call(rbind,theta_boot)
	return(
		cbind(
			apply(boot,2,function(x){quantile(x,0.025)}),
			apply(boot,2,function(x){quantile(x,0.975)})
		)
	)

}


res = list()
res_tmp = list()
LL3 = function(X, m, s, t)(1/((X-t)*s*(2*pi)^0.5))*exp(((-(log(X-t)-m)^2)/(2*s^2)))


N = ***
sample_n = ***
sample_m = ***
sample_l = ***
meanlog = ***
sdlog = ***
threshold = ***

for(k in 1:1000){
	# show(k)
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
	y_C = runif(sample_l, min = threshold, max = y_true[(sample_n+sample_m+1):N])
	
	gamma0 = min(min(y),min(y_LU),min(y_C))-0.1
	mu0 = (sum(log(y-gamma0))+sum(log(y_LU[,1]-gamma0))+sum(log(y_C-gamma0)))/length(y_true)
	sigma0 = sqrt((sum((log(y-gamma0)-mu0)^2)+sum((log(y_LU[,1]-gamma0)-mu0)^2)+sum((log(y_C-gamma0)-mu0)^2))/length(y_true))
	
	theta0 = c(mu0, sigma0, gamma0)
	
	theta = try(EM_imcubation_accel(y,y_LU,y_C,theta0,epsilon=0.001,interval=c(threshold*2,gamma0),max_iter=10000)[1])
	if(class(theta) == "try-error"){
		y_true=rlnorm3(N,meanlog=meanlog,sdlog=sdlog,threshold=threshold)
	
		y = y_true[1:sample_n]

		y_LU = cbind( 
			runif(sample_m, min = threshold,max = y_true[(sample_n+1):(sample_n+sample_m)]),
			runif(sample_m, min = y_true[(sample_n+1):(sample_n+sample_m)],max = y_true[(sample_n+1):(sample_n+sample_m)]+1) 
			)
		# y_C = y_true[(sample_n+sample_m+1):N]-runif(sample_l,0,sdlog)
		y_C = runif(sample_l, min = threshold, max = y_true[(sample_n+sample_m+1):N])
		
		gamma0 = min(min(y),min(y_LU),min(y_C))-0.1
		mu0 = (sum(log(y-gamma0))+sum(log(y_LU[,1]-gamma0))+sum(log(y_C-gamma0)))/length(y_true)
		sigma0 = sqrt((sum((log(y-gamma0)-mu0)^2)+sum((log(y_LU[,1]-gamma0)-mu0)^2)+sum((log(y_C-gamma0)-mu0)^2))/length(y_true))
		
		theta0 = c(mu0, sigma0, gamma0)
		
		theta = try(EM_imcubation_accel(y,y_LU,y_C,theta0,epsilon=0.001,interval=c(threshold*2,gamma0),max_iter=10000)[1])
	}
	cov = Boot_Cov(y,y_LU,y_C,L=1000)
	
	set = c(N,sample_n,sample_m,sample_l,meanlog,sdlog,threshold)
	
	# m: Location Parameter
	# s: Scale Parameter
	# t: Threshold Parameter
	# com1 = try(fitdistr(x=y, densfun=LL3, start=list(m=mu0, s=sigma0, t=min(y)-0.1),method="L-BFGS-B",upper=c(Inf,Inf,min(y)-0.001),lower=c(-Inf,0.0001,-Inf)))
	com1 = try(optim(par=theta0,ll2,y=y,method="L-BFGS-B",upper=c(Inf,Inf,min(y)-0.01),lower=c(min(y)+0.01,0.0001,-Inf),hessian=T))

	if(class(com1) == "try-error"){
		com1 = try(optim(par=theta0,ll2,y=y,method="L-BFGS-B",upper=c(Inf,Inf,min(y)-0.1),lower=c(min(y)+0.1,0.0001,-Inf),hessian=T))
		
		if(class(com1) != "try-error"){
			if(all(eigen(com1$hessian)$values>0)){
				res_tmp[[k]] = list(
						set = set,
						theta = unlist(theta),
						cov = cov,
						c1_theta = com1$par,
						c1_vcov = sqrt(diag(solve(com1$hessian)))
				)
			}else{
				sd = sqrt(diag(solve(Matrix::nearPD(com1$hessian)$mat)))
				res_tmp[[k]] = list(
					set = set,
					theta = unlist(theta),
					cov = cov,
					c1_theta = com1$par,
					c1_vcov = sd
				)
			}
		}else if(class(com1) == "try-error"){
			res_tmp[[k]] = list(
				set = set,
				theta = unlist(theta),
				cov = cov,
				c1_theta = NA,
				c1_vcov = NA)			
		}
	}else{
		if(all(eigen(com1$hessian)$values>0)){
			res_tmp[[k]] = list(
					set = set,
					theta = unlist(theta),
					cov = cov,
					c1_theta = com1$par,
					c1_vcov = sqrt(diag(solve(com1$hessian)))
			)
		}else{
			sd = sqrt(diag(solve(Matrix::nearPD(com1$hessian)$mat)))
			res_tmp[[k]] = list(
				set = set,
				theta = unlist(theta),
				cov = cov,
				c1_theta = com1$par,
				c1_vcov = sd
			)
		}
	}

}


