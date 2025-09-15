data {
  int<lower=1> N;
  vector<lower=0>[N] x;
  int<lower=1> K;
}
parameters {
  real<lower=0, upper=min(x)> tau1;       // 最も左端のtau
  positive_ordered[K-1] delta_tau;        // 各tauの差分（順序制約付き）
  real mu;                                // 共通のmu
  real<lower=0> sigma;                    // 共通のsigma
  vector<lower=0, upper=1>[K-1] v;        // Stick-breaking用
}
transformed parameters {
  vector[K] tau;
  simplex[K] pi;
  vector[K] log_pi;

  tau[1] = tau1;
  for (k in 2:K)
    tau[k] = tau[k-1] + delta_tau[k-1];  // tauは明確に昇順

  // Stick-breaking
  pi[1] = v[1];
  real stick = 1 - v[1];
  for (k in 2:(K-1)) {
    pi[k] = v[k] * stick;
    stick *= (1 - v[k]);
  }
  pi[K] = stick;
  log_pi = log(pi);
}
model {
  // 明確な事前分布（識別性を高める）
  tau1 ~ uniform(0, min(x));
  delta_tau ~ exponential(1);
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  v ~ beta(1, 1);

  for (n in 1:N) {
    vector[K] lps;
    for (k in 1:K) {
      if (x[n] > tau[k])
        lps[k] = log_pi[k] + lognormal_lpdf(x[n] - tau[k] | mu, sigma);
      else
        lps[k] = negative_infinity();
    }
    target += log_sum_exp(lps);
  }
}
