data {
  int<lower = 1> N;
  real x[N];
  int<lower = 1> nu;
}
parameters {
  real mu;
  real<lower = 0> sigma;
}
model {
  x ~ student_t(nu, mu, sigma);
}
