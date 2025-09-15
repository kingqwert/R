library(outbreaks)
library(dplyr)
library(EpiCurve)
library(rstan)
library(loo)
library(tidyverse)

fit_model =stan_model(file="./adaptive_shifted_knots.stan")
model <- stan_model(file = "./comparison.stan")
fit_model2 =stan_model(file="./shifted_lognormal_mix.stan")


calc_quantile = function(alpha,mu,sigma){
	return(
		exp(mu+sigma*qnorm(alpha))
		)
}


######################
# MERS in Korea
######################
data(mers_korea_2015)
head(mers_korea_2015$linelist[, c("id", "dt_onset")])
head(mers_korea_2015$contacts)

id_list= paste0("SK_",c(
	14,11,16,18,12,20,6,27,3,10,13,26,33,17,21,32,15,32,22,19,28,2,34,7,5, # contact with SK_1 in hosp room or hosp visit (maybe contact at 5/16~17)
	113,116,41,112,100,114,136,47,110,122,120,98,102,111,134,97,139,142,96,103,104,123,131,124,132,135,125,138,140,137,35,146,40,48,39,118,126, # contact with SK_14 in emergency room (maybe contact at 5/27)
	30,82,54,127,38,23,24,31,45,83,85,86,143,128,87,149,84,106,129,130 # contact with SK_16 in hospital room (maybe contact at 5/28)
	))

data1 = mers_korea_2015$contacts %>%
	filter(from %in% c('SK_1','SK_14','SK_16')) %>%
	filter(to %in% id_list) %>%
	left_join(mers_korea_2015[[1]][,c('id','dt_onset')],by = c('to' = 'id')) %>%
	group_by(dt_onset,from) %>%
	summarize(num = n())

EpiCurve(data1,
date = "dt_onset",
freq = "num",
cutvar = "from", colors=c("black",'gray80','gray60'),
period = "3 day",
ylabel="Number of cases",
xlabel=sprintf("From %s to %s", min(data1$dt_onset), max(data1$dt_onset)),
title = "Epidemic Curve",
note = "Daily epidemic curve")


data2 = mers_korea_2015$contacts %>%
	filter(from %in% c('SK_1','SK_14','SK_16')) %>%
	left_join(mers_korea_2015[[1]][,c('id','dt_onset')],by = c('to' = 'id')) %>%
	select(to, dt_onset) %>%
	mutate(days = as.integer(difftime(dt_onset, "2015/05/10",units = "days")))%>%
	filter(to %in% id_list)

x = data2[,3]


# our method
dens <- density(x, n=1000)

tau_grid <- dens$x[dens$x <= max(x)] 
tau_prior_density <- dens$y[dens$x <= max(x)]
tau_prior_density <- tau_prior_density / sum(tau_prior_density)

G <- length(tau_grid)
N = length(x)

stan_data <- list(
  N = N, 
  x = x, 
  K = 10,
  G = G, 
  tau_grid = tau_grid,
  tau_prior_density = tau_prior_density
)


# Run  
fit_our =sampling(
  fit_model,
  data=stan_data,
  open_progress=FALSE,
  show_messages=FALSE,
  chains = 2,
  cores = 1,
  warmup=1000,
  iter=2000,
  seed = 123)



print(fit_our, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
calc_quantile(0.95,1.82,0.49)


###########
# comparison fixed K
# change K from 1-5 with minimum waic
###########
waic_result_K = numeric()
fit_comp1_K = list()
for(i in 1:5){
  stan_data <- list(N = N, x = x, K = i)
  fit_comp1_K[[i]] <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)

  log_lik_K <- extract_log_lik(fit_comp1_K[[i]], parameter_name = "log_lik", merge_chains = FALSE)
  # WAICの計算
  waic_result_K[i] <- waic(log_lik_K)$estimates[3,1]
}

fit_comp1 = fit_comp1_K[[which.min(waic_result_K)]]
print(fit_comp1, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
calc_quantile(0.95,2.47,0.31)

###########
# comparison fixed K = 3
# K is correctly specified
###########
stan_data <- list(N = N, x = x, K = 3)
fit_comp1_K3 <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)
print(fit_comp1_K3, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
calc_quantile(0.95,2.65,0.25)

###########
# comparison fixed K = 2
###########
stan_data <- list(N = N, x = x, K = 2)
fit_comp1_K2 <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)
print(fit_comp1_K2, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
calc_quantile(0.5,2.47,0.31)

###########
# comparison fixed K = 4
###########
stan_data <- list(N = N, x = x, K = 4)
fit_comp1_K4 <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)
print(fit_comp1_K4, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
calc_quantile(0.95,2.71,0.22)

############
# comaprison flexible K, but not adaptive knots
############
K = 10
data_list <- list(
  N = N,
  x = x,
  K = K  
)

fit_comp2 =sampling(
  fit_model2,
  data=data_list,
  open_progress=FALSE,
  show_messages=FALSE,
  chains = 2,
  cores = 1,
  warmup=1000,
  iter=2000,
  seed = 123)
print(fit_comp2, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
calc_quantile(0.95,1.95,0.43)



