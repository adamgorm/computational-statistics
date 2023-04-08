library(splines)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("implementation.cpp")

## Forces evaluation of all arguments
force_all_args <- function()
    as.list(parent.frame())

inv_logit <- function(x)
    exp(x) / (1 + exp(x))

decay_scheduler <- function(gamma0 = 1, a = 1, K = 1, gamma1, n1)
{
    force_all_args()
    if (!missing(gamma1) && !missing(n1))
        K <- n1^a * gamma1 / (gamma0 - gamma1)
    b <- gamma0 * K
    function(n) b / (K + n^a)
}

get_knots <- function(inner_knots)
    sort(c(rep(range(inner_knots), 3), inner_knots))

Omega <- function(inner_knots)
{
    knots <- sort(c(rep(range(inner_knots), 3), inner_knots))
    d <- diff(inner_knots)  # The vector of knot differences; b - a
    g_ab <- splineDesign(knots, inner_knots, derivs = 2)
    knots_mid <- inner_knots[-length(inner_knots)] + d / 2
    g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
    g_a <- g_ab[-nrow(g_ab), ]
    g_b <- g_ab[-1, ]
    (crossprod(d * g_a,  g_a) +
     4 * crossprod(d * g_ab_mid, g_ab_mid) +
     crossprod(d * g_b, g_b)) / 6
}

Phi <- function(x_vec, knots)
    splineDesign(knots, x_vec)

d2f_two_norm_squared <- function(beta, inner_knots)
    crossprod(beta, Omega(inner_knots) %*% beta)

p_old_slow <- function(beta, x_vec, knots)
    inv_logit(Phi(x_vec, knots) %*% beta)

p <- function(beta, phi_mat)
    inv_logit(phi_mat %*% beta)

grad_loglik <- function(beta, lambda, x_vec, y_vec, knots, inner_knots)
    - crossprod(Phi(x_vec, knots),
                y_vec - p(beta, x_vec, knots)) / length(x_vec) +
        lambda * 2 * grad_d2f_two_norm_squared(beta, inner_knots)

get_H <- function(lambda, x_vec, y_vec, knots, inner_knots)
{
    force_all_args()
    N <- length(x_vec)
    Omega_mat <- Omega(inner_knots)
    Phi_mat <- Phi(x_vec, knots)
    function(beta) {
        p_beta <- p(beta, Phi_mat)
        -1/N * (crossprod(y_vec, log(p_beta)) + 
                crossprod(1 - y_vec, log(1 - p_beta))) + 
            lambda * crossprod(beta, Omega_mat %*% beta)
    }
}

get_grad_H_slow <- function(lambda, x_vec, y_vec, knots, inner_knots)
{
    force_all_args()
    Phi_mat <- Phi(x_vec, knots)
    Omega_mat <- Omega(inner_knots)
    function(beta, i) {
        Phi_i <- matrix(Phi_mat[i,], nrow = length(i))
        -crossprod(Phi_i,
                   y_vec[i] - p(beta, Phi_i)) /
            length(x_vec[i]) +
            2 * lambda * Omega_mat %*% beta
    }
}

get_grad_H <- function(lambda, x_vec, y_vec, knots, inner_knots)
{
    force_all_args()
    Phi_mat <- Phi(x_vec, knots)
    unpen_grad <- function(beta, i) {
        n_i <- length(i)
        Phi_i <- matrix(Phi_mat[i, ], nrow = n_i)
        -crossprod(Phi_i, y_vec[i] - p(beta, Phi_i)) / n_i
    }
    if (lambda > 0) {
        two_lambda_Omega <- 2 * lambda * Omega(inner_knots)
        function(beta, i)
            unpen_grad(beta, i) + two_lambda_Omega %*% beta
    } else
        unpen_grad
}

get_grad_H_cpp_wrap <- function(lambda, x_vec, y_vec, knots, inner_knots)
{
    force_all_args()
    Phi_mat <- Phi(x_vec, knots)
    if (lambda > 0) {
        two_lambda_Omega <- 2 * lambda * Omega(inner_knots)
        function(beta, i)
            unpen_grad_cpp(beta, i, lambda, y_vec, Phi_mat) +
                two_lambda_Omega %*% beta
    } else
        function(beta, i)
            unpen_grad_cpp(beta, i, lambda, y_vec, Phi_mat)
}

sgd <- function(par, grad, n_obs, decay_schedule, epoch = batch,
                n_iter = 100, sampler = sample, cb = NULL, ...)
{
    learning_rates <- decay_schedule(1:n_iter)
    for(k in 1:n_iter) {
        if(!is.null(cb))
            cb()
        epoch_sample <- sampler(n_obs)
        par <- epoch(par, epoch_sample, learning_rates[k], grad, ...)
    }
    if (!is.null(cb))
        cb()
    par
}

## batch epoch

batch <- function(par, epoch_sample, learning_rate, grad,
                  batch_size = 50, ...)
{
    n_batches <- floor(length(epoch_sample) / batch_size)
    for(j in 0:(n_batches - 1)) {
        i <- epoch_sample[(j * batch_size + 1):
                          (j * batch_size + batch_size)]
        par <- par - learning_rate * grad(par, i, ...)
    }
    par
}

sgd_cpp_wrap <- function(par, x, y, decay_schedule, inner_knots,
                         lambda = 0, n_iter = 100, batch_size = 1)
{
    learning_rates <- decay_schedule(1:n_iter)
    Om <- Omega(inner_knots)
    Phi_mat <- Phi(x, get_knots(inner_knots))
    sgd_cpp(par, learning_rates, n_iter, batch_size, lambda, y, Phi_mat, Om)
}

momentum <- function() {
    rho <- 0
    function(par, epoch_sample, learning_rate, grad,
             batch_size = 50, mom_memory = 0.95, ...)
    {
        M <- floor(length(epoch_sample) / batch_size)
        for(j in 0:(M - 1)) {
            i <- epoch_sample[(j * batch_size + 1):
                      (j * batch_size + batch_size)]
            ## Using '<<-' assigns the value to rho in the enclosing environment
            rho <<- mom_memory * rho + (1 - mom_memory) * grad(par, i, ...)  
            par <- par - learning_rate * rho
        }
        par
    }
}



## Newton Raphson

get_newton_update <- function(n_param, x, y, lambda)
{
    force_all_args()
    n_sim <- length(x)
    inner_knots <- seq(min(x), max(x), length.out = n_param - 2)
    knots <- get_knots(inner_knots)
    Phi_mat <- Phi(x, knots)
    tPhi <- t(Phi_mat)
    Omega_mat <- Omega(inner_knots)
    
    function(beta)
    {
        p_vec <- as.vector(p(beta, Phi_mat))
        p_prod <- p_vec * (1 - p_vec)
        W <- diag(p_prod)
        W_inv <- diag(1 / p_prod)
        z <- Phi_mat %*% beta + W_inv %*% (y - p(beta, Phi_mat))
        solve(crossprod(Phi_mat, W %*% Phi_mat) + n_sim * lambda * Omega_mat) %*%
            crossprod(Phi_mat, W %*% z)
    }
}


newton <- function(init_guess, x, y, max_iter = 100, epsilon = 1e-5,
                   lambda = 0, cb = NULL) {
    n_param <- length(init_guess)
    par_new <- init_guess
    newton_update <- get_newton_update(n_param, x, y, lambda)
    for (i in 1:max_iter) {
        if (!is.null(cb))
            cb()
        par_old <- par_new
        par_new <- newton_update(par_old)
        if(sum((par_new - par_old)^2) <= epsilon * (sum(par_new^2) + epsilon)) 
            break
    }
    if (!is.null(cb))
        cb()
    par_new
}

##### Old versions

sgd_old <- function(n_iter, beta_init, lambda, knots, inner_knots, x, y,
                    decay_schedule = decay_scheduler(), batch_size = 1, cb = NULL)
{
    beta_new <- beta_init
    n <- length(x)
    for (i in 1:n_iter) {
        beta_old <- beta_new
        indexes <- sample(n, batch_size, replace = TRUE)
        x_vec <- x[indexes]
        y_vec <- y[indexes]
        beta_new <- beta_old - decay_schedule(i) *
            grad_loglik(beta_old, lambda, x_vec, y_vec,
                   knots, inner_knots)
        if (!is.null(cb)) cb()
    }
    beta_new
}

sgd_single_old <- function(n_iter, beta_init, lambda, knots, inner_knots, x, y,
                           decay_schedule = decay_scheduler(), cb = NULL)
{
    beta_new <- beta_init
    n <- length(x)
    for (i in 1:n_iter) {
        beta_old <- beta_new
        indexes <- sample(n, batch_size, replace = TRUE)
        x_vec <- x[indexes]
        y_vec <- y[indexes]
        beta_new <- beta_old - decay_schedule(i) *
            grad_loglik(beta_old, lambda, x_vec, y_vec,
                   knots, inner_knots)
        if (!is.null(cb)) cb()
    }
    beta_new
}

## SGD NRH

sgd_nrh <- function(
  par, 
  grad,              # Function of parameter and observation index
  N,                 # Sample size
  gamma,             # Decay schedule or a fixed learning rate
  maxiter = 100,     # Max epoch iterations
  sampler = sample,  # How data is resampled. Default is a random permutation
  cb = NULL, 
  ...
) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter) 
  for(k in 1:maxiter) {
    samp <- sampler(N)   
    if(!is.null(cb)) cb()
    for(j in 1:N) {
      i <- samp[j]
      par <- par - gamma[k] * grad(par, i, ...)
    }
  }
  par
}

epoch_simple <- function(par, samp, gamma, grad)
{
    N <- length(samp)
    for(j in 1:N) {
        i <- samp[j]
        par <- par - gamma * grad(par, i)
    }
    par
}

## Old implementation without epochs
## The tuning parameters must be changed for this, since they will now update in
## each step, instead of in each epoch.
sgd_batch <- function(par, n_obs, batch_size, decay_schedule, grad,
                      n_iter = 100, cb = NULL, ...)
{
    learning_rates <- decay_schedule(1:n_iter)
    for(k in 1:n_iter) {
        if(!is.null(cb))
            cb()
        i <- sample(n_obs, batch_size, replace = FALSE)
        par <- par - decay_schedule(k) * grad(par, i, ...)
    }
    par
}

sgd_old <- function(n_iter, beta_init, lambda, knots, inner_knots, x, y,
                    decay_schedule = decay_scheduler(), batch_size = 1, cb = NULL)
{
    beta_new <- beta_init
    n <- length(x)
    for (i in 1:n_iter) {
        beta_old <- beta_new
        indexes <- sample(n, batch_size, replace = TRUE)
        x_vec <- x[indexes]
        y_vec <- y[indexes]
        beta_new <- beta_old - decay_schedule(i) *
            grad_loglik(beta_old, lambda, x_vec, y_vec,
                   knots, inner_knots)
        if (!is.null(cb)) cb()
    }
    beta_new
}
