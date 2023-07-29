functions {
  real shifted_lognormal_lpdf(real y, real mu, real sigma, real shift) {
    return lognormal_lpdf(y - shift | mu, sigma);
  }
}

data {
  int<lower=0> N; // number of observations
  real y[N];      // observed data
  vector[N] censoring_limits; // censoring limits for right censoring
  vector[N] upper_limits; // upper limits for right censoring in interval censoring
  vector[N] lower_limits; // lower limits for left censoring in interval censoring
  int<lower=-1,upper=2> cens[N]; // censoring status
  real min_y;     // minimum observed value
  real theta;     // mean of prior for mu
  real theta1;    // mean of prior for sigma
  real tau;       // variance of prior for mu
  // real tau1;      // variance of prior for sigma
}

parameters {
  real mu;       // location parameter
  real<lower=0.001> sigma;    // scale parameter
  real<upper=min_y> shift; // shift parameter
}

model {
  mu ~ normal(theta, tau);     // prior for mu
  sigma ~ exponential(theta1);  // prior for sigma
  shift ~ uniform(-10, min_y);  // prior for shift

  /* 
  for (n in 1:N) {
    if (cens[n] == 0) {
        if (shift > y[n])
          target += negative_infinity();
        else
          target += shifted_lognormal_lpdf(y[n] | mu, sigma, shift);
    } else if (cens[n] == -1) {
        if (shift > censoring_limits[n]) 
          target += negative_infinity();
        else
          target += lognormal_lccdf(censoring_limits[n] - shift | mu, sigma);    
    } else if (cens[n] == 1) {
        if (shift > lower_limits[n]) 
          target += negative_infinity();
        else
          target += lognormal_lcdf(upper_limits[n] - shift | mu, sigma) - lognormal_lcdf(lower_limits[n] - shift | mu, sigma);
    }
  }
  */

  for (n in 1:N) {
    if (cens[n] == 0) {
        target += shifted_lognormal_lpdf(y[n] | mu, sigma, shift);
    } else if (cens[n] == -1) {
        target += lognormal_lccdf(censoring_limits[n] - shift | mu, sigma);    
    } else{
        target += lognormal_lcdf(upper_limits[n] - shift | mu, sigma) - lognormal_lcdf(lower_limits[n] - shift | mu, sigma);
    }
  }
}
