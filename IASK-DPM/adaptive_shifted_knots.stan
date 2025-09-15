data {
  int<lower=1> N;
  vector<lower=0>[N] x;
  int<lower=1> K;
  int<lower=1> G;                // グリッド数
  vector[G] tau_grid;            // tauのグリッド（事前密度用）
  vector[G] tau_prior_density;   // tauの事前密度（密度推定に基づく）
}
parameters {
  ordered[K] tau;                // 推定対象のtau
  real mu;
  real<lower=0> sigma;
  vector<lower=0, upper=1>[K-1] v;
}
transformed parameters {
  simplex[K] pi;
  pi[1] = v[1];
  real stick = 1 - v[1];
  for (k in 2:(K-1)) {
    pi[k] = v[k] * stick;
    stick *= (1 - v[k]);
  }
  pi[K] = stick;
}
model {
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  v ~ beta(1, 1);

  // Adaptiveな事前分布（密度推定ベースの離散近似事前分布）
  for (k in 1:K) {
    vector[G] log_p;
    for (g in 1:G) {
      log_p[g] = log(tau_prior_density[g]) - square(tau[k] - tau_grid[g]) / (2 * 0.1^2);
    }
    target += log_sum_exp(log_p); // 密度が高い位置にtauを集中させる
  }

  for (n in 1:N) {
    vector[K] lps;
    for (k in 1:K) {
      if (x[n] > tau[k])
        lps[k] = log(pi[k]) + lognormal_lpdf(x[n] - tau[k] | mu, sigma);
      else
        lps[k] = negative_infinity();
    }
    target += log_sum_exp(lps);
  }
}
