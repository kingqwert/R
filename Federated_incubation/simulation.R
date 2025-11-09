library(truncnorm)
library(metafor)
library(dplyr)

## ---------- simulate one dataset -----------------------------------
sim_data <- function(J=40,
                      mu0=2, tau_mu=0.4,
                      Avgn_j=40,
                      a0=0.5, b0=0.3,
                      pR=0.15, pI=0.15) {

  mu_j  <- rnorm(J, mu0, tau_mu)
  sig_j <- rgamma(J, shape = a0, rate = 1/b0)

  data  <- vector("list", J)
  for (j in seq_len(J)) {
    n  <- rpois(1, Avgn_j) + 5
    y  <- rnorm(n, mu_j[j], sig_j[j])

    ## --- right-censoring threshold (facility-constant) ---
    qR <- qnorm(1 - pR)
    A  <- mu_j[j] + sig_j[j] * qR

    ## --- who is right-censored? ---
    isR <- (y > A)                 # logical index length n
    nR  <- sum(isR)

    ## start all as fully observed
    obs <- data.frame(type = rep("O", n), y = y, u = y)

    ## --- apply right-censoring (safe even if nR == 0) ---
    if (nR > 0) {
      obs[isR, ] <- data.frame(type = rep("R", nR),
                               y    = rep(A, nR),
                               u    = rep(A, nR))
    }

    ## --- interval-censoring indices chosen from NON-right-censored cases ---
    avail <- which(!isR)
    nI_desired <- round(pI * n)
    nI <- max(0, min(length(avail), nI_desired))

    if (nI > 0) {
      idx_I <- sample(avail, nI, replace = FALSE)
      ## construct [L, U] with L < U
      L_val <- y[idx_I] - runif(nI, 0, 1)
      U_val <- y[idx_I] + runif(nI, 0, 1)
      swap <- which(L_val > U_val)
      if (length(swap) > 0) {
        tmp <- L_val[swap]; L_val[swap] <- U_val[swap]; U_val[swap] <- tmp
      }
      obs[idx_I, ] <- data.frame(type = rep("I", nI),
                                 y    = L_val,
                                 u    = U_val)
    }

    data[[j]] <- obs
  }
  list(data = data, mu = mu_j, sigma = sig_j)
}




## ---------- Two-step meta (DerSimonian–Laird) -----------------------
two_step_meta_mean <- function(d){
  ybar <- sapply(d, function(df) mean(df[df$type=="O",]$y))
  s2   <- sapply(d, function(df) var(df[df$type=="O",]$y))
  n    <- sapply(d, function(df) nrow(df[df$type=="O",]))
  vi   <- s2/n
  res  <- rma.uni(yi = ybar, vi = vi, method="DL")
  list(mu_hat = as.numeric(res$b),
       tau2   = res$tau2)
}


meta_variance <- function(data,
                          sd_col    = "sd",    # 標準偏差列
                          n_col     = "n",     # 標本サイズ列
                          vi_method = c("exact", "approx")) {    
  df <- data %>% 
    mutate(
      var_obs = .data[[sd_col]]^2,
      k       = .data[[n_col]] - 1,
      ## バイアス補正付き効果量
      yi      = log(var_obs) - (digamma(k / 2) - log(k / 2)),
      ## サンプリング分散
      vi      = if (vi_method == "exact") trigamma(k / 2) else 2 / k
    )
  
  ## ランダム効果モデル（REML + Hartung–Knapp 補正）
  res <- rma(yi, vi, data = df, method = "REML", test = "knha")
  
  ## 母分散推定値（対数→元スケールに戻す）
  sigma2_hat  <- exp(res$b)
  sigma2_hat  <- as.numeric(exp(res$b + 0.5*res$tau2))
  # sigma2_ci   <- exp(c(res$ci.lb, res$ci.ub))
  
  # list(model = res, pooled_sigma2 = sigma2_hat, ci_sigma2 = sigma2_ci, data = df)
  list(model = res, pooled_sigma2 = sigma2_hat, data = df)
}




fac_loglik_plugin <- function(df, mu, lambda){
  sig <- 1 / sqrt(lambda)
  out <- 0
  ## fully observed
  idxO <- which(df$type=="O")
  if(length(idxO)){
    y <- df$y[idxO]
    out <- out + sum(dnorm(y, mu, sig, log=TRUE))
  }
  ## right-censored
  idxR <- which(df$type=="R")
  if(length(idxR)){
    a <- df$y[idxR]
    za <- (a - mu)/sig
    out <- out + sum(pnorm(za, lower.tail=FALSE, log.p=TRUE))
  }
  ## interval
  idxI <- which(df$type=="I")
  if(length(idxI)){
    L <- df$y[idxI]; U <- df$u[idxI]
    zL <- (L - mu)/sig; zU <- (U - mu)/sig
    ## numerically stable log{Phi(zU)-Phi(zL)}
    lpU <- pnorm(zU, log.p=TRUE)
    lpL <- pnorm(zL, log.p=TRUE)
    out <- out + sum(lpU + log1p(-exp(lpL - lpU)))
    #out <- out + sum(lpU-lpL)
  }
  return(out)
}


log1mexp <- function(x) {
  out <- numeric(length(x))
  cut <- log(2)
  i1 <- (x <= cut)
  out[i1]  <- log(-expm1(-x[i1]))      # 小さいとき
  out[!i1] <- log1p(-exp(-x[!i1]))     # 大きいとき
  out
}


logspace_sub <- function(a, b) {
  out_logabs <- numeric(length(a))
  out_sign   <- integer(length(a))
  idx <- (a >= b)
  # a >= b
  if (any(idx)) {
    d <- a[idx] - b[idx]               # >= 0
    out_logabs[idx] <- a[idx] + log1mexp(d)
    out_sign[idx]   <-  1L
  }
  # a < b
  if (any(!idx)) {
    d <- b[!idx] - a[!idx]             # >= 0
    out_logabs[!idx] <- b[!idx] + log1mexp(d)
    out_sign[!idx]   <- -1L
  }
  list(logabs = out_logabs, sign = out_sign)
}



#######################################################################
## Closed-form KLs
#######################################################################
kl_norm_norm <- function(m, v, M, V){
  0.5 * ( (v + (m - M)^2) / V + log(V / v) - 1 )
}
kl_gamma_gamma <- function(aq, bq, ap, bp){
  ## Ga(shape=a, rate=b)
  (aq - ap) * digamma(aq) - lgamma(aq) + lgamma(ap) +
    ap * log(bq / bp) + aq * (bp - bq) / bq
}


## Facility data-term plug-in ------------------------------------
elbo_facility_dataterm <- function(df, m_j, lam_bar){
  sig <- 1 / sqrt(lam_bar)
  out <- 0
  if(any(df$type=="O")){
    y  <- df$y[df$type=="O"]
    out <- out + sum(dnorm(y, mean = m_j, sd = sig, log = TRUE))
  }
  if(any(df$type=="R")){
    a  <- df$y[df$type=="R"]
    za <- (a - m_j)/sig
    out <- out + sum(pnorm(za, lower.tail = FALSE, log.p = TRUE))
  }
  if(any(df$type=="I")){
    L  <- df$y[df$type=="I"]; U <- df$u[df$type=="I"]
    zL <- (L - m_j)/sig; zU <- (U - m_j)/sig
    ## log( Phi(zU) - Phi(zL) ) in numerically stable form
    logdiff <- pnorm(zU, log.p = TRUE) + log1p(-exp(pnorm(zL, log.p = TRUE) - pnorm(zU, log.p = TRUE)))
    out <- out + sum(logdiff)
  }
  out
}

## Facility ELBO (approx) ----------------------------------------
elbo_facility <- function(df, m_j, v_j, a_j, b_j, M, V, A, B,
                          alpha0, beta0){
  lam_bar <- a_j / b_j
  n <- nrow(df)

  ## data term (plug-in)
  lt <- elbo_facility_dataterm(df, m_j, lam_bar)

  ## prior terms
  ## p(mu_j | mu0, eta): mu0 mean=M, eta mean=A/B -> plug-in
  eta_bar <- A / B
  lt <- lt + dnorm(m_j, mean = M, sd = 1/sqrt(eta_bar), log = TRUE)  # plug-in
  ## p(lambda_j): Ga(alpha0,beta0) expected log under q? -> use E_q[log p]
  lt <- lt + (alpha0-1)*digamma(a_j) - (alpha0-1)*log(b_j) - beta0 * (a_j/b_j) +
            alpha0*log(beta0) - lgamma(alpha0)

  ## -E_q[log q(mu_j)] -E_q[log q(lambda_j)]
  lt <- lt - ( -0.5*log(2*pi*v_j) -0.5 )  # entropy of N(m_j,v_j) = 0.5*(1+log 2πv)
  lt <- lt - ( a_j - log(b_j) + lgamma(a_j) + (1-a_j)*digamma(a_j) ) # entropy of Ga

  lt
}

## Global ELBO -------------------------------------------
elbo_global <- function(M,V,A,B,m0,s0,a0,b0){
  ## -KL q(mu0)||p + -KL q(eta)||p
  -kl_norm_norm(M,V,m0,s0^2) - kl_gamma_gamma(A,B,a0,b0)
}

## Stochastic ELBO (mini-batch) --------------------------
elbo_stoch <- function(S_idx, d, m_j, v_j, a_j, b_j,M,V,A,B, hyper){
  al0 <- hyper$alpha0
  bt0 <- hyper$beta0
  J <- length(d)
  Bsz <- length(S_idx) 
  kappa <- J/Bsz
  contribs <- vapply(S_idx, function(j)
    elbo_facility(d[[j]], m_j[j], v_j[j], a_j[j], b_j[j],
                  M,V,A,B, al0, bt0),
    numeric(1))
  kappa * sum(contribs) + elbo_global(M,V,A,B,hyper$m0,hyper$s0,hyper$a0,hyper$b0)
}




## ---------- Central CAVI -----------------------------------
central_cavi <- function(d, maxit = 200, tol_elbo= 1e-8, tol_par = 1e-6,
                         alpha0=1, beta0=1,
                         m0 = 0, s0 = 5,
                         a0 = 1, b0 = 1,
                         verbose = FALSE) {

  J <- length(d)
  
  ## init
  m_j <- sapply(d, function(df) mean(df$y))
  v_j <- rep(1, J)
  a_j <- rep(alpha0+1, J) 
  b_j <- rep(beta0+1, J)
  M <- m0
  V <- s0^2  
  A <- a0+1 
  B <- b0+1
  
  elbo_old <- -100000
  elbo_hist <- numeric(maxit)
  
  for(it in 1:maxit) {
    mu_old <- m_j
    ## local
    bar_eta <- A/B
    for (j in seq_len(J)) {
      df <- d[[j]]
      lam_j <- a_j[j]/b_j[j]
      mu <- m_j[j] 
      sig <- sqrt(1/lam_j)
      ## moments
      # EY <- ifelse(df$type=="O", df$y,
      #         ifelse(df$type=="R",
      #                mu + sig*dnorm((df$y-mu)/sig)/(1-pnorm((df$y-mu)/sig)),
      #                {L=df$y; U=df$u;
      #                 a=(L-mu)/sig; b=(U-mu)/sig; p=pnorm(b)-pnorm(a)
      #                 mu + sig*(dnorm(a)-dnorm(b))/p
      #                }))

      # EY2 <- ifelse(df$type=="O", df$y^2,
      #         ifelse(df$type=="R",
      #                { z=(df$y-mu)/sig; m1=sig*dnorm(z)/(1-pnorm(z))
      #                  mu^2+sig^2+ (mu+df$y)*m1 },
      #                {L=df$y;U=df$u;a=(L-mu)/sig;b=(U-mu)/sig;p=pnorm(b)-pnorm(a)
      #                 m1 = sig*(dnorm(a)-dnorm(b))/p
      #                 mu^2+sig^2 + (mu+L)*sig*dnorm(a)/p - (mu+U)*sig*dnorm(b)/p
      #                }))
      EY  <- numeric(nrow(df)); EY2 <- numeric(nrow(df))
      idxO <- df$type=="O"; idxR <- df$type=="R"; idxI <- df$type=="I"
      
      ## O
      if(any(idxO)){ y <- df$y[idxO]; EY[idxO] <- y; EY2[idxO] <- y^2 }
      
      ## R: Y>a
      if(any(idxR)){
        a <- df$y[idxR]; z <- (a-mu)/sig
        logphi  <- dnorm(z, log=TRUE)
        logtail <- pnorm(z, lower.tail=FALSE, log.p=TRUE)  # 安定
        r <- exp(logphi - logtail)                          # Mills ratio
        EZ  <- r
        EZ2 <- 1 + z*r
        EY[idxR]  <- mu + sig*EZ
        EY2[idxR] <- mu^2 + 2*mu*sig*EZ + sig^2*EZ2
      }
      
      ## I: l<Y<u
    ## I: interval l<Y<u
    if (any(idxI)) {
      L <- df$y[idxI]; U <- df$u[idxI]
      a <- (L - mu)/sig; b <- (U - mu)/sig
    
      lpU <- pnorm(b, log.p = TRUE)
      lpL <- pnorm(a, log.p = TRUE)
      # log{Phi(b) - Phi(a)} を安定に
      big <- (lpU >= lpL)
      lden <- numeric(length(lpU))
      lden[big]    <- lpU[big]    + log1mexp(lpU[big] - lpL[big])
      lden[!big]   <- lpL[!big]   + log1mexp(lpL[!big] - lpU[!big])
    
      lphi_a <- dnorm(a, log = TRUE)
      lphi_b <- dnorm(b, log = TRUE)
    
      # r1 = (phi(a) - phi(b)) / p
      d1 <- logspace_sub(lphi_a, lphi_b)
      r1 <- d1$sign * exp(d1$logabs - lden)
    
      # r2 = (a*phi(a) - b*phi(b)) / p  （符号は a, b の符号に依存）
      logA <- log(abs(a)) + lphi_a
      logB <- log(abs(b)) + lphi_b
      # logA <- ifelse(!is.finite(a) | a <= 0, -Inf, log(abs(a)) + lphi_a)
      # logB <- ifelse(!is.finite(b) | b <= 0, -Inf, log(abs(b)) + lphi_b)
      sA   <- sign(a); sB <- sign(b)
    
      # 2項の和差を安全に：符号が同じなら差、異なれば和
      same <- (sA == sB)
      r2 <- numeric(length(a))
      if (any(same)) {
        tmp <- logspace_sub(logA[same], logB[same])
        r2[same] <- sA[same] * exp(tmp$logabs - lden[same])
      }
      if (any(!same)) {
        # 符号が異なるときは logsumexp
        m  <- pmax(logA[!same], logB[!same])
        lse <- m + log( exp(logA[!same] - m) + exp(logB[!same] - m) )
        # sA=+1,sB=-1 → 正； sA=-1,sB=+1 → 負
        ssum <- ifelse(sA[!same] > 0, 1, -1)
        r2[!same] <- ssum * exp(lse - lden[!same])
      }
    
      EZ  <- r1
      EZ2 <- 1 + r2
    
      EY[idxI]  <- mu + sig * EZ
      EY2[idxI] <- mu^2 + 2*mu*sig*EZ + sig^2 * EZ2
    
      # きわめて狭い区間で p ≈ 0 の保険（任意）
      tiny <- (lden < -40)                  # p < ~4e-18
      if (any(tiny)) {
        # 近似：中心点でO扱い、などにフォールバック
        ymid <- 0.5 * (L[tiny] + U[tiny])
        EY[tiny]  <- ymid
        EY2[tiny] <- ymid^2
      }
    }
      S1 <- sum(EY)
      S2 <- sum(EY2)
      n <- nrow(df)

      v_j[j] <- 1/(n*lam_j + bar_eta)
      m_j[j] <- v_j[j]*(lam_j*S1 + bar_eta*M)
      a_j[j] <- alpha0 + n/2
      # b_j[j] <- beta0 + 0.5*(S2 - 2*m_j[j]*S1 + n*(m_j[j]^2+v_j[j]))
      ss_term <- S2 - 2*m_j[j]*S1 + n*(m_j[j]^2+v_j[j])
      b_j[j] <- beta0 + 0.5 * pmax(ss_term, 1e-9) # pmaxで負の値を回避
    }

    ## global
    lam_bar <- a_j/b_j
    V_new  <- 1/(J*bar_eta + 1/s0^2)
    M_new  <- V_new * (bar_eta * sum(m_j) + m0/s0^2)  
    A_new  <- a0 + J/2
    B_new  <- b0 + 0.5*(sum(m_j^2+v_j) - 2*M_new*sum(m_j) + J*(M_new^2+V_new))

    M <- M_new; V <- V_new; A <- A_new; B <- B_new

    ## ---- ELBO (deterministic plug-in) -------------------------------
    elbo <- 0
    for (j in seq_len(J)) {
      lam_bar <- a_j[j] / b_j[j]
      elbo <- elbo + fac_loglik_plugin(d[[j]], m_j[j], lam_bar)
      elbo <- elbo - kl_norm_norm(m_j[j], v_j[j],M, 1 / (A / B))   
      elbo <- elbo - kl_gamma_gamma(a_j[j], b_j[j], alpha0, beta0)
      #show(c(elbo, m_j[j], a_j[j],b_j[j],lam_bar))
    }
    elbo <- elbo + elbo_global(M,V,A,B,m0,s0,a0,b0)

    elbo_hist[it] <- elbo
    if (verbose) message(sprintf("central it=%d ELBO=%.6f", it, elbo))

    relchg <- abs(elbo - elbo_old) / (abs(elbo_old) + 1e-12)
    dpar   <- max(abs(c(M_new - M, V_new - V, A_new - A, B_new - B)))  # (here zero; shown for symmetry)

    if (relchg < tol_elbo & dpar < tol_par) {
      elbo_hist <- elbo_hist[seq_len(it)]
      break
    }
    elbo_old <- elbo  
  }
  list(M=M, V=V, A=A, B=B, m_j=m_j, v_j=v_j,a_j=a_j,b_j=b_j,elbo=tail(elbo_hist,1), elbo_hist=elbo_hist)
}



## ---------- Fed-SVI --------------------
fed_svi <- function(d, rounds = 80, Bfrac = 0.2,
                    rho_k = 0.7, t0 = 10,
                    tol_elbo = 1e-5, tol_par = 1e-4, window = 3,
                    alpha0=1,beta0=1,m0=0,s0=5,a0=1,b0=1){

  J <- length(d) 
  B <- ceiling(Bfrac * J)

  hyper <- list(alpha0=alpha0,beta0=beta0,m0=m0,s0=s0,a0=a0,b0=b0)
  
  ## init global
  M <- m0
  V <- s0^2
  A <- a0 + 1
  Bpar <- b0 + 1

  ## init locals
  a_j <- rep(alpha0+2, J)
  b_j <- rep(beta0+2, J)
  m_j <- sapply(d, function(df) mean(df$y))
  v_j <- rep(1, J)
  
  elbo_hist <- numeric(rounds)
  elbo_sm <- NA

  conv_count <- 0
  
  for(t in 0:(rounds-1)){
    rho <- (t0 + t)^(-rho_k)
    S_t <- sample.int(J, B)
    
    ## local updates
    bar_eta <- A/Bpar
    for(j in S_t){
      df <- d[[j]]
      lam_j <- a_j[j]/b_j[j]
      mu <- m_j[j]
      sig <- sqrt(1/lam_j)
      
      # EY <- ifelse(df$type=="O", df$y,
      #         ifelse(df$type=="R",
      #                mu + sig*dnorm((df$y-mu)/sig)/(1-pnorm((df$y-mu)/sig)),
      #                {L=df$y; U=df$u;
      #                 a=(L-mu)/sig; b=(U-mu)/sig; p=pnorm(b)-pnorm(a)
      #                 mu + sig*(dnorm(a)-dnorm(b))/p
      #                }))

      # EY2 <- ifelse(df$type=="O", df$y^2,
      #         ifelse(df$type=="R",
      #                { z=(df$y-mu)/sig; m1=sig*dnorm(z)/(1-pnorm(z))
      #                  mu^2+sig^2+ (mu+df$y)*m1 },
      #                {L=df$y;U=df$u;a=(L-mu)/sig;b=(U-mu)/sig;p=pnorm(b)-pnorm(a)
      #                 m1 = sig*(dnorm(a)-dnorm(b))/p
      #                 mu^2+sig^2 + (mu+L)*sig*dnorm(a)/p - (mu+U)*sig*dnorm(b)/p
      #                }))
      EY  <- numeric(nrow(df)); EY2 <- numeric(nrow(df))
      idxO <- df$type=="O"; idxR <- df$type=="R"; idxI <- df$type=="I"
      
      ## O
      if(any(idxO)){ y <- df$y[idxO]; EY[idxO] <- y; EY2[idxO] <- y^2 }
      
      ## R: Y>a
      if(any(idxR)){
        a <- df$y[idxR]; z <- (a-mu)/sig
        logphi  <- dnorm(z, log=TRUE)
        logtail <- pnorm(z, lower.tail=FALSE, log.p=TRUE)  # 安定
        r <- exp(logphi - logtail)                          # Mills ratio
        EZ  <- r
        EZ2 <- 1 + z*r
        EY[idxR]  <- mu + sig*EZ
        EY2[idxR] <- mu^2 + 2*mu*sig*EZ + sig^2*EZ2
      }
      
    ## I: interval l<Y<u
    if (any(idxI)) {
      L <- df$y[idxI]; U <- df$u[idxI]
      a <- (L - mu)/sig; b <- (U - mu)/sig
    
      lpU <- pnorm(b, log.p = TRUE)
      lpL <- pnorm(a, log.p = TRUE)
      # log{Phi(b) - Phi(a)} を安定に
      big <- (lpU >= lpL)
      lden <- numeric(length(lpU))
      lden[big]    <- lpU[big]    + log1mexp(lpU[big] - lpL[big])
      lden[!big]   <- lpL[!big]   + log1mexp(lpL[!big] - lpU[!big])
    
      lphi_a <- dnorm(a, log = TRUE)
      lphi_b <- dnorm(b, log = TRUE)
    
      # r1 = (phi(a) - phi(b)) / p
      d1 <- logspace_sub(lphi_a, lphi_b)
      r1 <- d1$sign * exp(d1$logabs - lden)
    
      # r2 = (a*phi(a) - b*phi(b)) / p  （符号は a, b の符号に依存）
      logA <- log(abs(a)) + lphi_a
      logB <- log(abs(b)) + lphi_b
      sA   <- sign(a); sB <- sign(b)
    
      # 2項の和差を安全に：符号が同じなら差、異なれば和
      same <- (sA == sB)
      r2 <- numeric(length(a))
      if (any(same)) {
        tmp <- logspace_sub(logA[same], logB[same])
        r2[same] <- sA[same] * exp(tmp$logabs - lden[same])
      }
      if (any(!same)) {
        # 符号が異なるときは logsumexp
        m  <- pmax(logA[!same], logB[!same])
        lse <- m + log( exp(logA[!same] - m) + exp(logB[!same] - m) )
        # sA=+1,sB=-1 → 正； sA=-1,sB=+1 → 負
        ssum <- ifelse(sA[!same] > 0, 1, -1)
        r2[!same] <- ssum * exp(lse - lden[!same])
      }
    
      EZ  <- r1
      EZ2 <- 1 + r2
    
      EY[idxI]  <- mu + sig * EZ
      EY2[idxI] <- mu^2 + 2*mu*sig*EZ + sig^2 * EZ2
    
      # きわめて狭い区間で p ≈ 0 の保険（任意）
      tiny <- (lden < -40)                  # p < ~4e-18
      if (any(tiny)) {
        # 近似：中心点でO扱い、などにフォールバック
        ymid <- 0.5 * (L[tiny] + U[tiny])
        EY[tiny]  <- ymid
        EY2[tiny] <- ymid^2
      }
    }
      S1 <- sum(EY)
      S2 <- sum(EY2)
      n <- nrow(df)
     
      v_j[j] <- 1/(n*lam_j + bar_eta)
      m_j[j] <- v_j[j]*(lam_j*S1 + bar_eta*M)
      a_j[j] <- alpha0 + n/2
      # b_j[j] <- beta0 + 0.5*(S2 - 2*m_j[j]*S1 + n*(m_j[j]^2+v_j[j]))
      
      ss_term <- S2 - 2*m_j[j]*S1 + n*(m_j[j]^2+v_j[j])
      b_j[j] <- beta0 + 0.5 * pmax(ss_term, 1e-9) 
    }
    
    ## secure-aggregate
    cfac  <- J/length(S_t)
    Sm    <- cfac * sum(m_j[S_t])
    Smv   <- cfac * sum(m_j[S_t]^2 + v_j[S_t])

    ## ---------- construct instantaneous natural params -------------
    bar_eta <- A / Bpar                   # E_q[eta] BEFORE update
    th1_hat <- bar_eta * Sm + m0 / s0^2
    th2_hat <- -0.5 * (J * bar_eta + 1 / s0^2)
    th3_hat <- a0 + 0.5 * J - 1           # constant; could skip update
    th4_hat <- -( b0 + 0.5 * (Smv - 2 * M * Sm + J * (M^2 + V)) )

    ## current natural params
    th1 <- M / V
    th2 <- -0.5 / V
    th3 <- A - 1
    th4 <- -Bpar

    ## ---------- Robbins–Monro natural-gradient step ---------------
    th1_new <- (1 - rho) * th1 + rho * th1_hat
    th2_new <- (1 - rho) * th2 + rho * th2_hat
    th3_new <- (1 - rho) * th3 + rho * th3_hat  # or simply th3_hat
    th4_new <- (1 - rho) * th4 + rho * th4_hat

    ## ---------- back-transform to moment params -------------------
    V_prev <- V; M_prev <- M; A_prev <- A; B_prev <- Bpar  # save for convergence
    
    V <- -0.5 / th2_new
    M <- th1_new * V
    A <- th3_new + 1
    Bpar <- -th4_new

    ## ---- stochastic ELBO & stopping ----
    elbo_hat <- elbo_stoch(S_t, d, m_j, v_j, a_j, b_j, M, V, A, Bpar, hyper)
    elbo_hist[t+1] <- elbo_hat

    if(is.na(elbo_sm)) {
      elbo_sm <- elbo_hat
      prev_elbo_sm <- elbo_sm
    } else {
      elbo_sm_old <- elbo_sm
      elbo_sm <- 0.9 * elbo_sm + 0.1 * elbo_hat  # EWMA
      relchg <- abs(elbo_sm - elbo_sm_old) / (abs(elbo_sm_old) + 1e-12)

      dpar <- max(abs(c(M - M_prev,
                        V - V_prev,
                        A - A_prev,
                        Bpar - B_prev)))

      if(relchg < tol_elbo && dpar < tol_par) {
        conv_count <- conv_count + 1L
      } else {
        conv_count <- 0L
      }
      if(conv_count >= window) break
    }

    prev_elbo_sm <- elbo_sm
  }
  list(M=M, V=V, A=A, B=Bpar, m_j=m_j, v_j=v_j,a_j=a_j,b_j=b_j,elbo=tail(elbo_hist,1), elbo_hist=elbo_hist)
}


## ---------- Fed-SVI on *fully observed only*  ------------------------
fed_svi_obs_only <- function(d, ...) {
  d_obs <- lapply(d, function(df) df[df$type=="O", , drop=FALSE])
  fed_svi(d_obs, ...)                # same hyper-parameters & rounds
}


qfun <- function(mu0, tau){
  q80  <- qnorm(0.025)   
  q95  <- qnorm(0.975)   
  return(c(p80 = mu0 + tau*q80,
              p95 = mu0 + tau*q95))
  }


Qhat_X <- function(fit, d, probs = c(0.025, 0.975) ){

  # between-facility variance tau^2
  tau2 <- fit$B / fit$A
  
  # facility-size weights
  n <- vapply(d, nrow, integer(1))

  # E[sigma_j^2] = E[1/lambda_j] = b_j/(a_j-1) (requires a_j>1)
  a <- fit$a_j;  b <- fit$b_j
  # s2j <- b / (a - 1)
  s2j <- b/a
  ok  <- is.finite(s2j) & (a > 1) & (n > 0)

  if (!any(ok)) {
    stop("No valid facility for E[sigma^2]: need some a_j>1 with n_j>0.")
  }
  w <- n[ok] / sum(n[ok])
  s2bar <- sum(w * s2j[ok])

  # quantiles on log scale, then exponentiate
  z <- qnorm(probs)
  # q <- exp(fit$M + z * sqrt(tau2 + s2bar))
  q <- fit$M + z * sqrt(tau2 + s2bar)
  names(q) <- paste0("p", formatC(probs*100, width = 0, format = "f"))
  q
}



eval_run <- function(J,mu0,tau_mu,Avgn_j,a0,b0,pR,pI){
  tryCatch({
    dat <- sim_data(J=J,mu0=mu0, tau_mu=tau_mu,
                  Avgn_j=Avgn_j,
                  a0=a0, b0=b0,
                  pR=pR, pI=pI)                     # 1 Monte-Carlo replicate
  
    dat_example <- data.frame(
      study = paste0("Study_", 1:length(dat$data)),
      sd    = sapply(dat$data,function(df) sd(df[df$type=="O",]$y)),
      n     = sapply(dat$data,function(df) nrow(df[df$type=="O",]))
    )
  
    oracle <- central_cavi(dat$data)
    fed    <- fed_svi(dat$data)
    fed_O  <- fed_svi_obs_only(dat$data)
  
    ## true hyper parameters from generated mu_j ----------------------
    mu_true  <- mu0
    tau_true <- tau_mu
    q_true   <- qfun(mu_true, sqrt(tau_mu^2 + a0*(1+a0)*b0^2)) 
  
    # meta-analysis
    meta_mean   <- two_step_meta_mean(dat$data)
    meta_var <- meta_variance(dat_example, sd_col = "sd", n_col = "n", vi_method = "exact")
  
    # calculate the expected alpha0, beta0
    oracle_var_mean = 1/sqrt(estimate_gamma_hyper(oracle$a_j,oracle$b_j)$alpha0/estimate_gamma_hyper(oracle$a_j,oracle$b_j)$beta0)
    fed_var_mean = 1/sqrt(estimate_gamma_hyper(fed$a_j,fed$b_j)$alpha0/estimate_gamma_hyper(fed$a_j,fed$b_j)$beta0)
    fed_O_var_mean = 1/sqrt(estimate_gamma_hyper(fed_O$a_j,fed_O$b_j)$alpha0/estimate_gamma_hyper(fed_O$a_j,fed_O$b_j)$beta0)
  
    c(                               
      ## point estimates of mu_0
      mu0_oracle = oracle$M,
      mu0_meta   = meta_mean$mu_hat,
      mu0_fed    = fed$M,
      mu0_fedO   = fed_O$M,
      mu0_true   = mu_true,
  
      ## tau (between-site SD)
      tau_oracle = 1/sqrt(oracle$A/oracle$B),
      tau_meta   = sqrt(meta_mean$tau2),
      tau_fed    = 1/sqrt(fed$A/fed$B),
      tau_fedO   = 1/sqrt(fed_O$A/fed_O$B),
      tau_true   = tau_true,
  
      ## 80 % / 95 % quantiles
      oracle025975 = Qhat_X(oracle, dat$data, probs = c(0.025, 0.975)),
      meta025975 = qfun(meta_mean$mu_hat, sqrt(meta_mean$tau2 + meta_var$pooled_sigma2)),
      fed025975 = Qhat_X(fed, dat$data, probs = c(0.025, 0.975)),
      fedO025975 = Qhat_X(fed_O, dat$data, probs = c(0.025, 0.975)),
      q_true = q_true      # truth
    )
  },error = function(e) {
      # Handle the error
      warning(paste("Error:", e$message))
      return(    
        c(                               
        mu0_oracle = NA,
        mu0_meta   = NA,
        mu0_fed    = NA,
        mu0_fedO   = NA,
        mu0_true   = NA,
        tau_oracle = NA,
        tau_meta   = NA,
        tau_fed    = NA,
        tau_fedO   = NA,
        tau_true   = NA,
        oracle025975 = NA,
        meta025975 = NA,
        fed025975 = NA,
        fedO025975 = NA,
        q_true = NA
      )) # Return NA on error
  })
}
