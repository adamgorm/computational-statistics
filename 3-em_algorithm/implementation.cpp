#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

struct estimates {
  double mu;
  double sigma2;
};

estimates
make_estimates(double mu, double sigma2)
{
  estimates est;
  est.mu = mu;
  est.sigma2 = sigma2;
  return est;
}

double
square(double x)
{
  return x * x;
}

double
two_norm(estimates est)
{
  return square(est.mu) + square(est.sigma2);
}

double
diff_two_norm(estimates new_est, estimates old_est)
{
  return square(new_est.mu - old_est.mu)
    + square(new_est.sigma2 - old_est.sigma2);
}

estimates
update_estimates(estimates old_estimates,
                 const NumericVector& x, double nu)
{
  estimates new_estimates = make_estimates(0, 0);
  double mu_divide_by = 0;  
  unsigned long long n_x = x.length();
  double beta_inverse[n_x];
  
  for (unsigned long long i = 0; i < n_x; i++) {
    beta_inverse[i] = 1 / (1 + square((x[i] - old_estimates.mu))
                        / (nu * old_estimates.sigma2));
    new_estimates.mu += x[i] * beta_inverse[i];
    mu_divide_by += beta_inverse[i];
  }
  new_estimates.mu /= mu_divide_by;
  for (unsigned long long i = 0; i < n_x; i++)
    new_estimates.sigma2 += square((x[i] - new_estimates.mu))
      * (1 + nu) * beta_inverse[i];
  new_estimates.sigma2 /= (n_x * nu);
  return new_estimates;
}

estimates
EM(double epsilon, const NumericVector& x, double nu,
   double mu_guess, double sigma2_guess)
{
  double epsilon_squared = square(epsilon);
  estimates old_est;
  estimates new_est = make_estimates(mu_guess, sigma2_guess);
  do {
    old_est = new_est;
    new_est = update_estimates(old_est, x, nu);
  } while (diff_two_norm(new_est, old_est)
           > epsilon_squared * square(two_norm(old_est) + epsilon));
  return new_est;
}

// [[Rcpp::export]]
NumericVector
EM_cpp(double epsilon, const NumericVector& x, double nu,
       double mu_guess, double sigma2_guess)
{
  estimates est = EM(epsilon, x, nu, mu_guess, sigma2_guess);
  NumericVector est_vec = {est.mu, est.sigma2};
  return est_vec;
}
