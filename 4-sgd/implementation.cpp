#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat
crossprod(const arma::mat& A, const arma::mat& B)
{
  return A.t() * B;
}

arma::colvec
inv_logit(const arma::colvec& x)
{
  return exp(x) / (1 + exp(x));
}

arma::colvec
p(const arma::colvec& beta, const arma::mat& Phi_i)
{
  return inv_logit(Phi_i * beta);
}

// [[Rcpp::export]]
arma::mat
unpen_grad_cpp(const arma::colvec& beta, arma::uvec i,
               double lambda, const arma::colvec& y,
               const arma::mat& Phi)
{
  i -= 1;                       // 0-indexing instead of 1-indexing.
  return crossprod(Phi.rows(i), p(beta, Phi.rows(i)) - y.elem(i)) / i.n_elem;
}

/* Complete C++ implementation */

/* Version with 0-indexing for use in C++ */
arma::colvec
unpen_grad(const arma::colvec& beta, arma::uvec i,
           double lambda, const arma::colvec& y,
           const arma::mat& Phi)
{
  return crossprod(Phi.rows(i), p(beta, Phi.rows(i)) - y.elem(i)) / i.n_elem;
}

arma::colvec
grad_cpp(const arma::colvec& beta, arma::uvec i,
         double lambda, const arma::colvec& y,
         const arma::mat& Phi, const arma::mat& Omega)
{
  if (lambda > 0)
    return unpen_grad(beta, i, lambda, y, Phi) + 2 * lambda * Omega;
  else
    return unpen_grad(beta, i, lambda, y, Phi);
}

void
batch(arma::colvec& par, double learning_rate, int batch_size,
      double lambda, arma::uvec& epoch_order, const arma::colvec& y,
      const arma::mat& Phi, const arma::mat& Omega)
{
  unsigned long long lower_i = 0;
  unsigned long long upper_i = batch_size - 1;
  unsigned long long n_batches = floor(y.n_elem / batch_size);
  for (unsigned long long j = 0; j < n_batches; j++) {
    par = par - learning_rate *
      grad_cpp(par, epoch_order.subvec(lower_i, upper_i), lambda, y, Phi, Omega);
    lower_i += batch_size;
    upper_i += batch_size;
  }
}

// [[Rcpp::export]]
arma::colvec
sgd_cpp(arma::colvec par, const arma::colvec& learning_rates,
        int n_iter, int batch_size, double lambda, const arma::colvec& y,
        const arma::mat& Phi, const arma::mat& Omega)
{
  arma::uvec epoch_order(y.n_elem);
  for (int k = 0; k < n_iter; k++) {
    epoch_order = arma::randperm(y.n_elem);
    batch(par, learning_rates[k], batch_size, lambda, epoch_order, y, Phi, Omega);
  }
  return par;
}

