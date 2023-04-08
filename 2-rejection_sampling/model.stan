data {
  int<lower = 0> N;
  int<lower = 0> z[N];
  vector[N] x;
}
parameters {
  real y;
}
model {
  z ~ poisson_log(y * x);
}
