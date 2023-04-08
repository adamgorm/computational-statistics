data {
  int<lower = 1> N;
  real x[N];
  real w[N];  
  int<lower = 1> nu;
}
parameters {
  real mu;
  real<lower = 0> sigma2;
}
model {
  w ~ chi_square(nu);
  for (i in 1:N)
    x[i] ~ normal(mu, sqrt(nu * sigma2 / w[i]));
}
