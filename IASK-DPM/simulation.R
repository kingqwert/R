library(rstan)
library(loo)

fit_model =stan_model(file="./adaptive_shifted_knots.stan")
model <- stan_model(file = ".//comparison.stan")
fit_model2 =stan_model(file="./shifted_lognormal_mix.stan")


# データ生成
set.seed(123)
N <- 1000
true_pi <- c(0.4, 0.35, 0.25)

# COVID-19
# Lauer SA, Grantz KH, Bi Q, et al.
# The Incubation Period of Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: Estimation and Application.
# Annals of Internal Medicine. 2020;172(9):577-582.
# 多数の事例から、対数正規分布を仮定した場合の潜伏期間の中央値が 5.1 日 (95%CI: 4.5–5.8日) と推定。
# true_mu <- 1.6
# true_sigma <- 0.5

# Measles
# CDC: Measles (Rubeola). Epidemiology and Prevention of Vaccine-Preventable Diseases (The Pink Book).
# 「潜伏期間は平均10～12日（範囲7～18日）程度」
# Lessler J, et al. “Incubation periods of acute respiratory viral infections: a systematic review.” Lancet Infect Dis. 2009;9:291–300.
# はしかの潜伏期間を含む、さまざまな呼吸器系ウイルス感染症の潜伏期間のメタ解析。
# Fine PEM, et al. “Transmission and control of measles.” Lancet. 1993;341(8852):449–452.
# はしかの流行特性や潜伏期間などの古典的レビュー論文。
# true_mu <- 2.4
# true_sigma <- 0.3


# O-157
# Mead PS, Griffin PM. Escherichia coli O157:H7. Lancet. 1998;352(9135):1207-1212
true_mu <- 1.25
true_sigma <- 0.57


res_x = list()
fit_our = list()
fit_comp1 = list()
fit_comp1_K3 = list()
fit_comp1_K2 = list()
fit_comp1_K4 = list()
fit_comp2 = list()
true_taus = list()

for(j in 1:500){
  
  # tauを前のmodeから生成（tau1 = 0）
  # modeはdelta + exp(mu - sigma^2)
  tau_1 <- 0
  tau_2 = exp(true_mu-true_sigma^2)
  tau_3 = tau_2 + exp(true_mu-true_sigma^2)
    
  # 最終的なデータ生成
  true_tau <- c(tau_1, tau_2, tau_3)
  
  # データを混合比率に基づいて生成
  z <- sample(1:3, N, replace = TRUE, prob = true_pi)
  x <- sapply(z, function(k) rlnorm(1, true_mu, true_sigma) + true_tau[k])
  res_x[[j]] = x
  #hist(x,breaks=100)
  
  
  #########
  # our method
  #########
  dens <- density(x, n=1000)
  
  # グリッドと事前密度（Stanで使うため）
  tau_grid <- dens$x[dens$x <= max(x)]  # 必ずmax(x)以下で制限
  tau_prior_density <- dens$y[dens$x <= max(x)]
  tau_prior_density <- tau_prior_density / sum(tau_prior_density)
  
  G <- length(tau_grid)
  
  # Stanデータ（完全版）
  stan_data <- list(
    N = N, 
    x = x, 
    K = 10,  # 混合成分数
    G = G, 
    tau_grid = tau_grid,
    tau_prior_density = tau_prior_density
  )
  
  ###########
  # comparison fixed K
  # Kを1-5まで動かして、waicでminimumをとる
  ###########
  # waic_result_K = numeric()
  # fit_comp1_K = list()
  # for(i in 1:5){
  #   stan_data <- list(N = N, x = x, K = i)
  #   fit_comp1_K[[i]] <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)
  
  #   log_lik_K <- extract_log_lik(fit_comp1_K[[i]], parameter_name = "log_lik", merge_chains = FALSE)
  #   # WAICの計算
  #   waic_result_K[i] <- waic(log_lik_K)$estimates[3,1]
  # }
  
  # fit_comp1[[j]] = fit_comp1_K[[which.min(waic_result_K)]]
  # print(fit_comp1, pars=c("pi","mu","sigma","tau"), probs=c(0.025,0.5,0.975))
  
  ###########
  # comparison fixed K = 3
  # K is correctly specified
  ###########
  # stan_data <- list(N = N, x = x, K = 3)
  # fit_comp1_K3[[j]] <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)

  # ###########
  # # comparison fixed K = 2
  # # K is misspecified
  # ###########
  # stan_data <- list(N = N, x = x, K = 2)
  # fit_comp1_K2[[j]] <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)

  # ###########
  # # comparison fixed K = 4
  # # K is misspecified
  # ###########
  # stan_data <- list(N = N, x = x, K = 4)
  # fit_comp1_K4[[j]] <- sampling(model, data = stan_data, iter = 2000, warmup=1000, chains = 2, seed = 123,open_progress=FALSE,show_messages=FALSE)


  ############
  # comaprison flexible K, but not adaptive knots
  ############
  K = 10
  data_list <- list(
    N = N,
    x = x,
    K = K
  )
  
  fit_comp2[[j]] =sampling(
    fit_model2,
    data=data_list,
    open_progress=FALSE,
    show_messages=FALSE,
    chains = 2,
    cores = 1,
    warmup=1000,
    iter=2000,
    seed = 123)
}

