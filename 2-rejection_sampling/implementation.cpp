#include <math.h>
#include <Rcpp.h>
#include <cstdio>
using namespace Rcpp;

/* Gaussian envelope */

double
log_q_over_p_cpp(double y, const NumericVector& x, double sum_xz)
{
  unsigned long long n = x.length();
  double result = y * y / 2 + y * sum_xz;
  
  for (unsigned long long i = 0; i < n; i++)
    result -= exp(y * x[i]);
  return result;
}

// [[Rcpp::export]]
NumericVector
sim_loop_cpp(unsigned long long n, double alpha_log,
             const NumericVector& x, double sum_xz)
{
  NumericVector y(n);
  double y_norm;
  unsigned long long count = 0;
  
  while (count < n) {
    y_norm = R::rnorm(0,1);
    if (log(R::runif(0,1)) <= alpha_log + log_q_over_p_cpp(y_norm, x, sum_xz))
      y[count++] = y_norm;
  }
  return y;
}

/* Adaptive envelopes */

double
logf(double y, const NumericVector& x, double sum_xz)
{
  unsigned long long n = x.length();
  double result = y * sum_xz;
  
  for (unsigned long long i = 0; i < n; i++)
    result -= exp(y * x[i]);
  return result;
}

double
dlogf(double y, const NumericVector& x, double sum_xz)
{
  unsigned long long n = x.length();
  double result = sum_xz;
  
  for (unsigned long long i = 0; i < n; i++)
    result -= x[i] * exp(y * x[i]);
  return result;
}

double *
get_a(const NumericVector& t, unsigned long long n_t,
      const NumericVector& x, double sum_xz)
{
  double *a = (double *) malloc(n_t * sizeof(double));
  
  for (unsigned long long i = 0; i < n_t; i++)
    *(a + i) = dlogf(t[i], x, sum_xz);
  return a;
}

double *
get_b(double *a, const NumericVector& t, unsigned long long n_t,
      const NumericVector& x, double sum_xz)
{
  double *b = (double *) malloc(n_t * sizeof(double));
  for (unsigned long long i = 0; i < n_t; i++)
    *(b + i) = logf(t[i], x, sum_xz) - *(a + i) * t[i];
  return b;
}

double *
get_z(double *a, double *b, unsigned long long n_ab)
{
  double *z = (double *) malloc((n_ab + 1) * sizeof(double));
  
  *z = 0;
  for (unsigned long long i = 1; i < n_ab; i++)
    *(z + i) = (*(b + i) - *(b + i - 1)) / (*(a + i - 1) - *(a + i));
  *(z + n_ab) = 1;
  return z;
}

/* Replaces contents of array with the cumulative sum */
void
replace_with_cumsum(double *array, unsigned long long array_length)
{
  for (unsigned long long i = 1; i < array_length; i++)
    *(array + i) += *(array + i - 1);
}

/* Returns pointer to array of length n_ab where element 0 is Q_0 = 0, element 1
   is Q_1, and so on, until element n_ab which is Q_m = c */
double *
get_Q_and_c(double *a, double *b, double *z, unsigned long long n_ab)
{
  double *Q_and_c = (double *) malloc((n_ab + 1) * sizeof(double));
  
  *Q_and_c = 0;
  for (unsigned long long i = 0; i < n_ab; i++)
    *(Q_and_c + i + 1) = exp(*(b + i))
      * (exp(*(a + i) * *(z + i + 1))
         - exp(*(a + i) * *(z + i))) / *(a + i);
  replace_with_cumsum(Q_and_c, n_ab + 1);
  return Q_and_c;
}

unsigned long long
find_interval(double point, double *grid, unsigned long long n_grid)
{
  if (*(grid + n_grid - 1) < point)
    return n_grid - 1;
      
  for (unsigned long long i = 0; i < n_grid - 1; i++)
    if (*(grid + i) < point && point <= *(grid + i + 1))
      return i;

  return -1;                    /* Error code */
}

double
get_accepted_proposal(double *a, double *b, double c,
                      double *z, double *Q, unsigned long long n_ab,
                      const NumericVector& x, double sum_xz)
{
  double cu, proposal;
  unsigned long long i;
  
  do {
    cu = c * R::runif(0, 1);
    i = find_interval(cu, Q, n_ab);
    proposal = log(*(a + i) * exp(-*(b + i)) * (cu - *(Q + i))
                   + exp(*(a + i) * *(z + i))) / *(a + i);
  } while (log(R::runif(0, 1)) >
           logf(proposal, x, sum_xz) - *(a + i) * proposal - *(b + i));

  return proposal;
}

// [[Rcpp::export]]
NumericVector
ctest_cpp(const NumericVector& y)
{
  double Q[] = {0, 0.3, 0.6, 1};
  unsigned long long n = y.length();
  NumericVector result(n);  
  for (unsigned long long i = 0; i < n; i++)
    result[i] = find_interval(y[i], Q, 4);
  return result;
}

// [[Rcpp::export]]
NumericVector
sim_adapt_more_cpp(unsigned long long n, const NumericVector& t,
                   const NumericVector& x, double sum_xz)
{
  unsigned long long n_t = t.length();
  double *a = get_a(t, n_t, x, sum_xz);
  double *b = get_b(a, t, n_t, x, sum_xz);
  double *z = get_z(a, b, n_t);
  double *Q_and_c = get_Q_and_c(a, b, z, n_t);
  double c = *(Q_and_c + n_t);
  NumericVector simulations(n);  
  for (unsigned long long i = 0; i < n; i++)
    simulations[i] = get_accepted_proposal(a, b, c, z, Q_and_c, n_t, x, sum_xz);
  return simulations;
}

/* Binary search implementation for benchmarking */

/* Uses binary search to find the interval containing point. */
unsigned long long
find_interval_binary_search(double point, double *grid, unsigned long long n_grid)
{
  unsigned long long i_min = 0;
  unsigned long long i_max = n_grid - 1;
  unsigned long long i;

  if (*(grid + i_max) <= point)
    return i_max;
  
  while (i_max - i_min > 1) {
    i = floor((i_max + i_min) / 2);
    if (*(grid + i) == point)
      return i;
    else if (*(grid + i) > point)
      i_max = i;      
    else
      i_min = i;
  }

  return i_min;
}

double
get_accepted_proposal_binary_search(double *a, double *b, double c,
                                    double *z, double *Q, unsigned long long n_ab,
                                    const NumericVector& x, double sum_xz)
{
  double cu, proposal;
  unsigned long long i;
  
  do {
    cu = c * R::runif(0, 1);
    i = find_interval_binary_search(cu, Q, n_ab);
    proposal = log(*(a + i) * exp(-*(b + i)) * (cu - *(Q + i))
                   + exp(*(a + i) * *(z + i))) / *(a + i);
  } while (log(R::runif(0, 1)) >
           logf(proposal, x, sum_xz) - *(a + i) * proposal - *(b + i));

  return proposal;
}

// [[Rcpp::export]]
NumericVector
ctest_cpp_binary_search(const NumericVector& y)
{
  double Q[] = {0, 0.3, 0.6, 1};
  unsigned long long n = y.length();
  NumericVector result(n);  
  for (unsigned long long i = 0; i < n; i++)
    result[i] = find_interval_binary_search(y[i], Q, 4);
  return result;
}

// [[Rcpp::export]]
NumericVector
sim_adapt_more_cpp_binary_search(unsigned long long n, const NumericVector& t,
                                 const NumericVector& x, double sum_xz)
{
  unsigned long long n_t = t.length();
  double *a = get_a(t, n_t, x, sum_xz);
  double *b = get_b(a, t, n_t, x, sum_xz);
  double *z = get_z(a, b, n_t);
  double *Q_and_c = get_Q_and_c(a, b, z, n_t);
  double c = *(Q_and_c + n_t);
  NumericVector simulations(n);  
  for (unsigned long long i = 0; i < n; i++)
    simulations[i] = get_accepted_proposal(a, b, c, z, Q_and_c, n_t, x, sum_xz);
  return simulations;
}
