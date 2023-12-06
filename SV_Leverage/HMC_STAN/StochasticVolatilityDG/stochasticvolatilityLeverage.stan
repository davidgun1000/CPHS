data {
  int<lower=0> T;   // # time points (equally spaced)
  matrix[T,1] y;      // mean corrected return at time t
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;// persistence of volatility
  real<lower=-1,upper=1> rho; // correlation parameter
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h;                 // log volatility at time t
}
model {
  (phi+1)/2 ~ beta(100, 1.5);
  atanh(rho) ~ normal(0,10000);
  sigma ~ cauchy(0, 1);
  mu ~ normal(0, 10000);
  h[1] ~ normal(mu, sigma / sqrt(1 - phi * phi));
  for (t in 2:T)
    h[t] ~ normal(mu + phi * (h[t - 1] -  mu) + sigma*rho*exp(-h[t-1]/2)*y[t-1], sigma*sqrt(1-rho^2));
  for (t in 1:T)
    y[t] ~ normal(0, exp(h[t] / 2));
}

