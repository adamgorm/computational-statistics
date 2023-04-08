#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

double
square(double x)
{
  return x * x;
}


/* Plug-in */

double sqrt2pi = sqrt(2 * PI);

double
d4_std_gaussian_density_single_cpp(double x)
{
  double x_sq = square(x);
  return exp(-x_sq/2) * ((x_sq - 6) * x_sq + 3) / sqrt2pi;
}

NumericVector
d4_std_gaussian_density_cpp(const NumericVector& x)
{
  unsigned long long n = x.length();
  NumericVector y(n);
  for (unsigned long long i = 0; i < n; i++)
    y[i] = d4_std_gaussian_density_single_cpp(x[i]);
  return y;
}

double d4_std_gaussian_density_0 = d4_std_gaussian_density_single_cpp(0);

// [[Rcpp::export]]
double
d2f0_two_norm_squared_estimate_cpp(const NumericVector& x, double r)
{
  unsigned long long n = x.length();
  double result = 0;
  double rsqrt2 = r * sqrt(2);
  for (unsigned long long i = 0; i < n - 1; i++)
    for (unsigned long long j = i + 1; i < j && j < n; j++)
      result += d4_std_gaussian_density_single_cpp((x[j] - x[i]) / rsqrt2);
  return (2 * result + n * d4_std_gaussian_density_0) /
    (square(n) * pow(rsqrt2, 5));
}


/* Cross validation */

double
ep_kernel_single_cpp(double x)
{
  if (-1 < x && x < 1)
    return (1 - square(x)) * 0.75;
  else
    return 0;
}

// [[Rcpp::export]]
NumericVector
ep_kernel_cpp(const NumericVector& x)
{
  unsigned long long n = x.length();
  NumericVector y(n);
  for (unsigned long long i = 0; i < n; i++)
    y[i] = ep_kernel_single_cpp(x[i]);
  return y;
}


/* Density estimate */

double
ep_density_single_cpp(double x, const NumericVector& x_obs, double h)
{
  double result = 0;
  unsigned long long n = x_obs.length();
  for (unsigned long long i = 0; i < n; i++)
    result += ep_kernel_single_cpp((x - x_obs[i]) / h);
  return result / (h * n);
}

// [[Rcpp::export]]
NumericVector
ep_density_cpp(const NumericVector& x_obs, double h,
               const NumericVector& grid)
{
  NumericVector density_estimates(grid.length());
  for (int i = 0; i < grid.length(); i++)
    density_estimates[i] = ep_density_single_cpp(grid[i], x_obs, h);
  return density_estimates;
}

/* Binning */

// [[Rcpp::export]]
NumericVector
bin_weights_cpp(const NumericVector& x, double grid_lower,
                int grid_length, double grid_diff)
{
  double closest_grid_i;
  NumericVector w(grid_length);  
  unsigned long long n_x = x.length();
  
  for (unsigned long long i = 0; i < n_x; i++) {
    closest_grid_i = floor((x[i] - grid_lower) / grid_diff + 0.5);
    w[closest_grid_i] += 1;
  }

  return w;
}

double
weighted_kernel_sum(double i, double max_nonzero_i, double grid_length,
                    const NumericVector& weights,
                    const double *kernel_evals)
{
  int below, above;
  double res = weights[i] * *kernel_evals;
  for (int j = 1; j <= max_nonzero_i; j++) {
    below = i - j;
    above = i + j;
    if (below >= 0)
      res += weights[below] * *(kernel_evals + j);
    if (above < grid_length)
      res += weights[above] * *(kernel_evals + j);
  }

  return res;
}
                    

// [[Rcpp::export]]
NumericVector
ep_density_binning_cpp(const NumericVector& x, double h,
                       const NumericVector& grid)
{
  int i;
  int n_x = x.length();
  int grid_length = grid.length();
  double grid_lower = grid[0];
  double grid_diff = grid[1] - grid[0];
  NumericVector y(grid_length);
  NumericVector weights = bin_weights_cpp(x, grid_lower, grid_length, grid_diff);
  int max_nonzero_i = floor(h / grid_diff);
  double kernel_evals[max_nonzero_i + 1];
  kernel_evals[0] = 0.75; // Ep kernel in 0.
  for (i = 1; i <= max_nonzero_i; i++)
    kernel_evals[i] = ep_kernel_single_cpp((i * grid_diff) / h);  
  for (i = 0; i < grid_length; i++)
    y[i] += weighted_kernel_sum(i, max_nonzero_i, grid_length,
                                weights, kernel_evals) / (h * n_x);
  return y;
}

// [[Rcpp::export]]
NumericVector
range_cpp(const NumericVector& x)
{
  unsigned long long n_x = x.length();
  NumericVector range(2);
  range[0] = x[0];
  range[1] = x[0];
  for (unsigned long long i = 1; i < n_x; i++) {
    if (x[i] < range[0])
      range[0] = x[i];
    else if (x[i] > range[1])
      range[1] = x[i];
  }
  return range;
}
