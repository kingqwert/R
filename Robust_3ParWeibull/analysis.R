library(ExtDist)
library(dplyr)
library(ggplot2)
library(stats)
library(weibullness)
library(ForestFit)
library(openxlsx)
library(EpiCurve)

############################################################
## R code: 3-parameter Weibull with gamma-divergence (MM algorithm)
## Complete example including data generation and iteration loop
############################################################
likelihood = function(alpha, beta, eta, y, gamma){
  return(
    -1/gamma*log(mean(d3pweibull(y=y,alpha=alpha,beta=beta,eta=eta)^gamma))+
    1/(1+gamma)*log(alpha^(gamma/beta)*beta^gamma*gamma(gamma*(beta-1)/beta+1)/(1+gamma)^(gamma*(beta-1)/beta+1))
    )
}

r3pweibull <- function(n, alpha, beta, eta){
  shape = beta
  scale = (1/alpha)^(1/beta) 
  return(rWeibull(n=n,params=list(shape = shape, scale = scale))+eta)
}

#-----------------------------------------------------------
# 2) PDF of 3-parameter Weibull distribution
#    Assumes y > eta, alpha>0, beta>0.
#-----------------------------------------------------------
d3pweibull <- function(y, alpha, beta, eta){
  z <- y - eta
  out <- alpha * beta * z^(beta - 1) * exp(-alpha * z^beta)
  return(out)
}

#-----------------------------------------------------------
# 3) One iteration step of the MM algorithm:
#    - Compute weights w_i(t)
#    - Update alpha^(t+1) from eq.(6)
#    - Solve eq.(7) & eq.(8) for beta^(t+1), eta^(t+1)
#-----------------------------------------------------------
w_cal = function(alpha_t, beta_t, eta_t, y, gamma){
  # Compute pdf_i(t)
  fval <- d3pweibull(y, alpha_t, beta_t, eta_t)
  
  numerator <- fval^gamma
  denom     <- sum(numerator)
  w <- numerator / denom
  return(w)
}

alpha_update = function(alpha_t, beta_t, eta_t, y, gamma){
  w = w_cal(alpha_t, beta_t, eta_t, y, gamma)
  numerator_a <- 1 - gamma / (beta_t * (1 + gamma))
  sum_den <- sum(w*(y - eta_t)^beta_t)  
  alpha_next <- numerator_a / sum_den
  return(alpha_next)
}

beta_update <- function(beta_next,alpha_t, beta_t, eta_t, y, gamma){  
  w = w_cal(alpha_t, beta_t, eta_t, y, gamma)

  return(
    abs(
      -sum(w*((1/beta_next) + (1-alpha_t*(y-eta_t)^beta_next)*log(y-eta_t)))+
    1/(1+gamma)*(-gamma/beta_next^2*(log(alpha_t*(1+gamma))-digamma(1+gamma*(beta_next-1)/beta_next))+gamma/beta_next)
    )
  )
}
beta_update2 <- function(beta_next,alpha_t, beta_t, eta_t, y, gamma){  
  w = w_cal(alpha_t, beta_t, eta_t, y, gamma)

  return(
    abs(-sum(w*((1/beta_next) + log(y-eta_t) - alpha_t*(y-eta_t)^beta_next*log(y-eta_t)))+
    1/(1+gamma)*(-gamma/beta_next^2*log(alpha_t)+gamma/beta_next-gamma/beta_next^2*log(1+gamma)+
    gamma/beta_next^2*digamma(1+gamma*(beta_next-1)/beta_next))
    )
  )
}
eta_update <- function(eta_next,alpha_t, beta_t, eta_t, y, gamma){  
  w = w_cal(alpha_t, beta_t, eta_t, y, gamma)

  # eq(8): Summation part
  return(
    abs(sum(w*((beta_t-1)/(y-eta_next)-alpha_t*beta_t*(y-eta_next)^(beta_t-1)))
    ))
}


mm_update <- function(alpha_t, beta_t, eta_t, y, gamma){
  y_min = min(y)-0.00001
  beta_min = gamma/(1+gamma)+0.0001
  #beta_min = 1 
  #beta_min = 0
  alpha_n = alpha_update(alpha_t, beta_t, eta_t, y, gamma)

  eta_n_tmp = optimize(eta_update, interval=c(-3,y_min), alpha_t=alpha_n,beta_t=beta_t,eta_t=eta_t,y=y,gamma=gamma)
  eta_n = eta_n_tmp$minimum

  beta_n_tmp = optimize(beta_update, interval=c(beta_min,10), alpha_t=alpha_n,beta_t=beta_t,eta_t=eta_n,y=y,gamma=gamma)
  beta_n = beta_n_tmp$minimum
  
  return(list(alpha=alpha_n, beta=beta_n, eta=eta_n))
}


#########################
# calculate H_n(ganmma)
#########################
g_y <- function(y, alpha, beta, eta) {
  return(beta/(y - eta) - alpha*beta)
}

fprime_y <- function(y, alpha, beta, eta) {
  return( g_y(y, alpha, beta, eta) * d3pweibull(y, alpha, beta, eta) )
}

f2prime_y <- function(y, alpha, beta, eta) {
  return( g_y(y, alpha, beta, eta)^2 * d3pweibull(y, alpha, beta, eta) )
}
calc_Cgamma <- function(alpha, beta, gamma) {
  inside <- alpha^(gamma/beta) *
            beta^gamma *
            (1+gamma)^(-(1+gamma - gamma/beta)) *
            gamma(1+gamma - gamma/beta) 
  val <- inside^( gamma/(1+gamma) )
  return(val)
}
Hn <- function(y, alpha, beta, eta, gamma) {
  Cg <- calc_Cgamma(alpha, beta, gamma)
  f_val   <- d3pweibull(y, alpha, beta, eta)         # f(y_i)
  f1_val  <- fprime_y(y, alpha, beta, eta)    # f'(y_i)
  f2_val  <- f2prime_y(y, alpha, beta, eta)   # f''(y_i)
    
  out <- sum(f1_val^2 * f_val^(gamma - 2) *
      (2*(gamma-1)/Cg + (f_val^gamma/Cg)) + 2*f_val^(gamma-1)*f2_val/Cg)
    
  return(out)
}




###########################
# Calculate Covariance matrix
###########################
## first deriv
partial_f_gamma_alpha <- function(x, alpha, beta, eta, gamma) {
  val <- - gamma * ( alpha*(x - eta)^beta - 1 ) / alpha * 
    d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

partial_f_gamma_beta <- function(x, alpha, beta, eta, gamma) {
  val <- - (gamma / beta) * ( beta * log(x - eta) * ( alpha*(x - eta)^beta - 1 ) - 1 ) *
    d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

partial_f_gamma_eta <- function(x, alpha, beta, eta, gamma) {
  val <- gamma * ( beta*( alpha*(x - eta)^beta - 1 ) + 1 ) / (x - eta) *
    d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

## second deriv
partial2_f_gamma_alpha2 <- function(x, alpha, beta, eta, gamma) {
  val <- gamma * ( gamma*( alpha*(x - eta)^beta - 1 )^2 - 1 ) / alpha^2 *
    d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

partial2_f_gamma_beta2 <- function(x, alpha, beta, eta, gamma) {
  f0 <- d3pweibull(x, alpha, beta, eta)^gamma  # 共通
  term1 <- beta^2 * log(x-eta)^2 * ( gamma*( alpha*(x-eta)^beta - 1 )^2 - alpha*(x-eta)^beta )
  term2 <- - 2 * gamma * beta * log(x-eta) * ( alpha*(x-eta)^beta - 1 )
  term3 <- (gamma - 1)
  
  bracket <- term1 + term2 + term3
  val <- (gamma / beta^2) * bracket * f0
  return(val)
}

partial2_f_gamma_eta2 <- function(x, alpha, beta, eta, gamma) {
  big1 <- beta*( alpha*(x - eta)^beta - 1 ) + 1
  big2 <- alpha*beta*(x - eta)^beta + 1
  
  val <- gamma * ( gamma*big1^2 - (beta - 1)*big2 ) / (x - eta)^2 *
    d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

partial2_f_gamma_alpha_beta <- function(x, alpha, beta, eta, gamma) {
  val <- (gamma / (alpha*beta)) * (
    - gamma * alpha*(x - eta)^beta
    + beta*log(x-eta)*(
      gamma*( alpha*(x - eta)^beta - 1 )^2 - alpha*(x - eta)^beta
    )
    + gamma
  ) * d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

partial2_f_gamma_alpha_eta <- function(x, alpha, beta, eta, gamma) {
  big1 <- alpha*(x - eta)^beta - 1
  big2 <- beta*( alpha*(x - eta)^beta - 1 ) + 1
  
  val <- - (gamma / (alpha*(x - eta))) * (
    gamma * big1 * big2
    - alpha*beta*(x - eta)^beta
  ) * d3pweibull(x, alpha, beta, eta)^gamma
  return(val)
}

partial2_f_gamma_beta_eta <- function(x, alpha, beta, eta, gamma) {
  f0 <- d3pweibull(x, alpha, beta, eta)^gamma

  termA <- beta - alpha*beta*(x - eta)^beta
  termB <- gamma * ( beta - alpha*beta*(x - eta)^beta - 1 )
  
  inside1 <- alpha*(x - eta)^beta - 1
  inside2 <- -beta + alpha*beta*(x - eta)^beta + 1
  partC1 <- gamma * ( inside1 * inside2 ) - alpha*beta*(x - eta)^beta
  
  termC <- beta * log(x-eta) * partC1
  
  bracket <- termA + termB + termC
  
  val <- - (gamma / (beta*(x - eta))) * bracket * f0
  return(val)
}


# grad l_gamma
grad_l_gamma <- function(alpha, beta, eta, gamma, y) {
  fvals_gamma <- d3pweibull(y, alpha, beta, eta)^gamma
  sum_fg <- sum(fvals_gamma)
  
  num_alpha <- sum( (alpha*(y-eta)^beta - 1) * fvals_gamma )
  dldalpha <- - num_alpha / sum_fg + gamma / (alpha * beta * (1+gamma))
  
  num_beta <- sum(( beta * log(y-eta) * ( alpha*(y-eta)^beta - 1 ) - 1 ) * fvals_gamma)
  part1_beta <- - (1/beta) * (num_beta / sum_fg)
  
  part2_beta <- (gamma/(1+gamma)) * (
    -1/(beta^2)*log(alpha) +
      1/beta -
      1/(beta^2)*log(1+gamma) +
      1/(beta^2)*digamma(1 + gamma - gamma/beta)
  )
  dldbeta <- part1_beta + part2_beta
  
  num_eta <- sum( 1/(y-eta) * ( beta*( alpha*(y-eta)^beta - 1 ) + 1 ) * fvals_gamma )
  dldeta <- - ( num_eta / sum_fg )
  
  out <- c(dldalpha, dldbeta, dldeta)
  names(out) <- c("dldalpha","dldbeta","dldeta")
  return(out)
  
}


##Hessian l_gamma
hess_l_gamma <- function(alpha, beta, eta, gamma, y) {
  fvals     <- d3pweibull(y, alpha, beta, eta)
  fvals_g   <- fvals^gamma
  sum_fg    <- sum(fvals_g)
  
  sum_dg_alpha <- sum(partial_f_gamma_alpha(y, alpha, beta, eta, gamma))
  sum_dg_beta  <- sum(partial_f_gamma_beta(y, alpha, beta, eta, gamma))
  sum_dg_eta   <- sum(partial_f_gamma_eta(y, alpha, beta, eta, gamma))
  
  sum_d2_alpha2    <- sum(partial2_f_gamma_alpha2(y, alpha, beta, eta, gamma))
  sum_d2_beta2     <- sum(partial2_f_gamma_beta2(y, alpha, beta, eta, gamma))
  sum_d2_eta2      <- sum(partial2_f_gamma_eta2(y, alpha, beta, eta, gamma))
  sum_d2_alpha_beta<- sum(partial2_f_gamma_alpha_beta(y, alpha, beta, eta, gamma))
  sum_d2_alpha_eta <- sum(partial2_f_gamma_alpha_eta(y, alpha, beta, eta, gamma))
  sum_d2_beta_eta  <- sum(partial2_f_gamma_beta_eta(y, alpha, beta, eta, gamma))
  

  coeff <- -1/gamma
  denom <- sum_fg^2
  
  H_alpha_alpha_core <- ( sum_fg * sum_d2_alpha2 - sum_dg_alpha^2 ) / denom
  H_alpha_beta_core  <- ( sum_fg * sum_d2_alpha_beta - sum_dg_alpha*sum_dg_beta ) / denom
  H_alpha_eta_core   <- ( sum_fg * sum_d2_alpha_eta  - sum_dg_alpha*sum_dg_eta  ) / denom
  H_beta_beta_core   <- ( sum_fg * sum_d2_beta2      - sum_dg_beta^2 ) / denom
  H_beta_eta_core    <- ( sum_fg * sum_d2_beta_eta   - sum_dg_beta*sum_dg_eta  ) / denom
  H_eta_eta_core     <- ( sum_fg * sum_d2_eta2       - sum_dg_eta^2 ) / denom
  
 
  H_aa <- coeff * H_alpha_alpha_core - gamma/( alpha^2 * beta * (1+gamma) )
  
  polyg   <- digamma( gamma + 1 - gamma/beta )
  triplyg <- trigamma( gamma + 1 - gamma/beta )
  
  bracket <- beta*( 2*log(alpha) - beta + 2*log(1+gamma) )-2*beta*polyg + gamma*triplyg
  
  H_bb <- coeff * H_beta_beta_core +
    ( gamma/(beta^4*(1+gamma)) ) * bracket
  H_ee <- coeff * H_eta_eta_core
  H_ab <- coeff * H_alpha_beta_core - gamma/( alpha*beta^2*(1+gamma) )
  H_ae <- coeff * H_alpha_eta_core
  H_be <- coeff * H_beta_eta_core
  
  
  H <- matrix(0, nrow=3, ncol=3)
  rownames(H) <- colnames(H) <- c("alpha","beta","eta")
  
  H["alpha","alpha"] <- H_aa
  H["beta","beta"]   <- H_bb
  H["eta","eta"]     <- H_ee
  
  H["alpha","beta"] <- H_ab
  H["beta","alpha"] <- H_ab
  
  H["alpha","eta"] <- H_ae
  H["eta","alpha"] <- H_ae
  
  H["beta","eta"] <- H_be
  H["eta","beta"] <- H_be
  
  return(H)
}


# calculate covariance matrix 
Cal_cov = function(alpha, beta, eta, gamma, y){
  H = hess_l_gamma(alpha, beta, eta, gamma, y)
  U = grad_l_gamma(alpha, beta, eta, gamma, y)%*%t(grad_l_gamma(alpha, beta, eta, gamma, y))
  out = solve(H) %*% U %*% solve(H)
  return(out)
}





######################
# Calculate the average
######################
mean_SHW <- function(alpha, beta, eta) { 
  return(
    (1/alpha)^(1/beta)*gamma(1+1/beta)+eta
    )
}
mean_SHWv <- function(vec){ 
  alpha = vec[1]
  beta = vec[2]
  eta = vec[3]
  return(
    (1/alpha)^(1/beta)*gamma(1+1/beta)+eta
    )
}
######################
# Calculate q% percentile
######################
q_SHW <- function(u,alpha, beta, eta) { 
  return(
    (-log(1-u)/alpha)^(1/beta)+eta
    )
}
q_SHWv <- function(u,vec) { 
  alpha = vec[1]
  beta = vec[2]
  eta = vec[3]
  return(
    (-log(1-u)/alpha)^(1/beta)+eta
    )
}






#-----------------------------------------------------------
# 4) The main iterative function:
#    We'll generate data from known alpha0,beta0,eta0
#    Then run the iteration until convergence or maxIter
#-----------------------------------------------------------
weibull_gdiv_estimate <- function(y, gamma=0.5, 
                                  alpha_init=1, beta_init=1.5, eta_init=0, 
                                  tol=1e-6, maxIter=100){
  
  alpha_curr <- alpha_init
  beta_curr  <- beta_init
  eta_curr   <- eta_init
  
  for(t in 1:maxIter){
    prev <- c(alpha_curr, beta_curr, eta_curr)
    l_old = likelihood(alpha_curr,beta_curr,eta_curr,y,gamma)

    # Perform one iteration
    out <- mm_update(alpha_curr, beta_curr, eta_curr, y, gamma)
    
    alpha_curr <- out$alpha
    beta_curr  <- out$beta
    eta_curr   <- out$eta
    
    # convergence check
    newv <- c(alpha_curr, beta_curr, eta_curr)
    l_new = likelihood(alpha_curr,beta_curr,eta_curr,y,gamma)
    # show(l_new)
    # show(c(alpha=alpha_curr, beta=beta_curr, eta=eta_curr, iter=t))
    diffparam <- abs(l_old-l_new)
    
    if(diffparam < tol){
      #cat(sprintf("Converged at iteration %d\n", t))

      break
    }
    #if(t==maxIter) cat("Warning: reached maxIter without full convergence.\n")
  }
  
  out <- c(alpha_curr, beta_curr, eta_curr,t)
  names(out) <- c("alpha","beta","eta","iter")
  return(out)
  # return(list(alpha=alpha_curr, beta=beta_curr, eta=eta_curr, iter=t))
}

#################
# start analysis for Wang and Teunis (2020)
#################
data = read.xlsx('~/data/application_Wang_and_Teunis(2020).xlsx')
y = data$num

#comparison groups
com_ml = weibull.mle(y)
com_moment = fitWeibull(y,location=TRUE, method = "moment",start=c(1,1,0))
com_mps = fitWeibull(y,location=TRUE, method = "mps",start=c(1,1,0))

mean_SHWv(
	c(1/as.numeric(com_ml[2])^as.numeric(com_ml[1]),as.numeric(com_ml[1]),as.numeric(com_ml[3]))
)
mean_SHWv(
    c(1/as.numeric(com_moment$estimate[2])^as.numeric(com_moment$estimate[1]),as.numeric(com_moment$estimate[1]),as.numeric(com_moment$estimate[3]))
)
mean_SHWv(
    c(1/as.numeric(com_mps$estimate[2])^as.numeric(com_mps$estimate[1]),as.numeric(com_mps$estimate[1]),as.numeric(com_mps$estimate[3]))
)

q_SHWv(
	0.95,c(1/as.numeric(com_ml[2])^as.numeric(com_ml[1]),as.numeric(com_ml[1]),as.numeric(com_ml[3]))
)
q_SHWv(
    0.95,c(1/as.numeric(com_moment$estimate[2])^as.numeric(com_moment$estimate[1]),as.numeric(com_moment$estimate[1]),as.numeric(com_moment$estimate[3]))
)
q_SHWv(
    0.95,c(1/as.numeric(com_mps$estimate[2])^as.numeric(com_mps$estimate[1]),as.numeric(com_mps$estimate[1]),as.numeric(com_mps$estimate[3]))
)


# our method 
Hns = list()
alpha_init = ifelse(1/as.numeric(com_ml[1])^as.numeric(com_ml[2])>10,1,1/as.numeric(com_ml[1])^as.numeric(com_ml[2]))
beta_init =ifelse(as.numeric(com_ml[1])>10,1,as.numeric(com_ml[1]))
eta_init = ifelse(abs(as.numeric(com_ml[3]))>10,0,as.numeric(com_ml[3]))

for(gamma_val in c(0.1,0.5,1,2)){
  res_tmp = try(weibull_gdiv_estimate(
    y, gamma=gamma_val, 
    alpha_init=alpha_init, beta_init=beta_init, eta_init=eta_init, 
    tol=1e-10, maxIter=10000
    ))
  if(class(res_tmp) == "try-error"){
    Hns[[which(gamma_val == c(0.1,0.5,1,2))]] = NA
  }else{
    Hns[[which(gamma_val == c(0.1,0.5,1,2))]] = Hn(y, alpha = res_tmp[1], beta=res_tmp[2], eta=res_tmp[3], gamma_val) 
  }
}
res = weibull_gdiv_estimate(
    y, gamma=c(0.1,0.5,1,2,3)[which.min(Hns)], 
    alpha_init=alpha_init, beta_init=beta_init, eta_init=eta_init, 
    tol=1e-10, maxIter=1000
)
Cal_cov(alpha = res[1], beta=res[2], eta=res[3], gamma=c(0.1,0.5,1,2,3)[which.min(Hns)], y)
c(res[3]-1.96*sqrt(2.2656257),res[3]+1.96*sqrt(2.2656257))

mean_SHWv(as.numeric(res[1:3]))
q_SHWv(0.95,as.numeric(res[1:3]))


