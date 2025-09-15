data {
  int<lower=1> N;        // 観測数
  vector<lower=0>[N] x;  // 観測データ（非負）
  int<lower=1> K;        // 混合成分数（外部から指定可能）
}
parameters {
  simplex[K] pi;      // 各成分の混合重み（和が1の確率ベクトル）
  real mu;               // 全成分共通の lognormal の平均パラメータ
  real<lower=0> sigma;   // 全成分共通の lognormal の標準偏差
  ordered[K] tau;        // 各成分のシフトパラメータ（識別性向上のため順序付け）
}
model {
  // --- 事前分布 ---
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  tau ~ normal(0, 10);
  
  // --- 尤度 ---
  // 各観測値 x[n] に対して、各混合成分の寄与を計算
  for (n in 1:N) {
    vector[K] lps;
    for (k in 1:K) {
      if (x[n] > tau[k])
        lps[k] = log(pi[k]) + lognormal_lpdf(x[n] - tau[k] | mu, sigma);
      else
        lps[k] = negative_infinity();  // サポート外なら尤度は無視
    }
    target += log_sum_exp(lps);
  }
}
generated quantities {
  // 各観測点の対数尤度を保存するための変数
  vector[N] log_lik;
  for (n in 1:N) {
    vector[K] lps;
    for (k in 1:K) {
      if (x[n] > tau[k])
        lps[k] = log(pi[k]) + lognormal_lpdf(x[n] - tau[k] | mu, sigma);
      else
        lps[k] = negative_infinity();
    }
    log_lik[n] = log_sum_exp(lps);
  }
}
