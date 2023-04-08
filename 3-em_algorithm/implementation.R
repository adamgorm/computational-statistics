library(numDeriv)

### Generalized solution

## Forces evaluation of all arguments
force_all_args <- function()
    as.list(parent.frame())

## The relative tolerance of the optimization affects the stopping criterion of
## the EM algorithm. If the relative tolerance of the optimization step is very
## small, then the algorithm will take smaller steps and thus reach the max
## allowed relative tolerance of the EM algorithm earlier as well (so it will
## not give as precise an estimate).
get_update_estimates_general <- function(Q)
{
    force_all_args()
    function(old_est)
        optim(old_est, function(par) -Q(par, old_est),
              control = list(reltol = 1e-16))$par
}

## Either Q or update_estimates must be supplied. If Q is supplied, numerical
## optimization is used to update the estimates.
EM <- function(epsilon, par_init, update_estimates = NULL, Q = NULL,
               callback = NULL)
{
    if (is.null(update_estimates))
        update_estimates <- get_update_estimates_general(Q)
    epsilon_squared <- epsilon^2
    new_est <- par_init
    not_converged <- TRUE
    while (not_converged) {
        old_est <- new_est
        new_est <- update_estimates(old_est)
        not_converged <- crossprod(new_est - old_est) >
            epsilon_squared * (crossprod(old_est) + epsilon)^2
        if (!is.null(callback))
            callback()
    }
    new_est
}

## Fisher information
## Q's second argument must be names 'par_prime'

fisher_def <- function(est, loglik, x)
    hessian(function(par) loglik(par, x = x), x = est)

fisher1 <- function(est, Q)
    -jacobian(function(par) grad(Q, par, par_prime = par), est)

D2Q <- function(est, Q)
    hessian(Q, est, par_prime = est)

fisher2 <- function(est, Q)
    -D2Q(est, Q) -
        jacobian(function(par_prime) grad(Q, est, par_prime = par_prime), est)

fisher3 <- function(est, Q, update_estimates = NULL)
{
    if (is.null(update_estimates))
        update_estimates <- get_update_estimates_general(Q)
    -(diag(length(est)) - t(jacobian(update_estimates, est))) %*%
        D2Q(est, Q)
}

emp_fisher <- function(mu, sigma2, x, nu)
{
    w_means <- (nu + 1) * 0.5 / (0.5 * (1 + ((x - mu)^2 /(sigma2 * nu))))
    dQ_dmu      <- w_means * ((x - mu) / (nu * sigma2))
    dQ_dsigma2 <- -1/(2*sigma2) + (w_means * (x - mu)^2) / (2 * nu * sigma2^2 ) 
    grad_list     <- cbind(mu = dQ_dmu, sigma2 = dQ_dsigma2)
    crossprod(grad_list)
}


### Functions specific to this model

sim_y <- function(n, mu, sigma2, nu)
{
    w <- rchisq(n, df = nu)
    x <- rnorm(n, mean = mu, sd = sqrt(nu * sigma2 / w))
    list(
        x = x,
        w = w
    )
}

## Maximum Likelihood Estimator
mle <- function(x, w, nu)
{
    mu_hat <- sum(w * x) / sum(w)
    c(mu_hat = mu_hat,
      sigma2_hat = mean(w * (x - mu_hat)^2) / nu)
}

get_update_estimates_specific <- function(x, nu)
{
    force_all_args()
    n <- length(x)
    function(old_est)
    {
        beta_inverse <- (1 + (x - old_est[1])^2 / (nu * old_est[2]))^(-1)
        mu_new <- sum(x * beta_inverse) / sum(beta_inverse)
        sigma2_new <- sum((x - mu_new)^2 * (1 + nu) * beta_inverse) / (n * nu)
        c(mu = mu_new, sigma2 = sigma2_new)
    }
}

## Only implemented up to additive constant, since we only need the derivatives
getQ <- function(x, nu)
{
    force(x)
    force(nu)
    function(par, par_prime)
        - length(x) * log(par[2]) / 2 -
            sum((x - par[1])^2 / (1 + (x - par_prime[1])^2 / (nu * par_prime[2]))) *
            (1 + nu) / (2 * nu * par[2])
}

### Gradient ascent

## Marginal log-likelihood and marginal gradient up to a constant.

get_marginal_loglik <- function(x, nu)
    function(par)
        -n/2 * log(par[2]) -
            (nu + 1)/2 * sum(log(1 + (x - par[1])^2 / (nu * par[2])))

get_marginal_grad <- function(x, nu)
    function(par)
        c(
            (nu + 1) * sum((x - par[1]) / ((x - par[1])^2 + nu * par[2])),
            nu/2 * sum( ((x - par[1])^2 - par[2]) /
                        (par[2] * (nu * par[2] + (x - par[1])^2)) )
        )

grad_ascent <- function(init_guess, objective, grad,
                        gamma = 0.01, epsilon = 1e-16,
                        callback = NULL)
{
    small_relative_ascent <- function(new_est, old_est)
    {
        objective_diff <- objective(new_est) - objective(old_est)
        objective_diff >= 0 &&
            objective_diff <= epsilon * (abs(objective(new_est)) + epsilon)
    }

    converged <- FALSE
    new_est <- init_guess
    while (!converged) {
        old_est <- new_est
        gr <- grad(old_est)
        new_est <- old_est + gamma * gr
        if (small_relative_ascent(new_est, old_est))
            converged <- TRUE
        if (!is.null(callback))
            callback()
    }

    new_est
}
