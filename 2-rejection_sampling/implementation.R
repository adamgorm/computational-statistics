library(tidyverse)

dat <- read.csv("Poisson.csv", header = TRUE)
x <- dat$x
sum_xz <- sum(x * dat$z)

f <- Vectorize(function(y)
    prod(exp(y * dat$z * dat$x - exp(y * dat$x))))

### Gaussian envelope

alpha_log_gauss <- scan("alpha_log.txt") # Precalculated

log_q_over_p_single <- function(y)
    y^2 / 2 + y * sum_xz - sum(exp(y * x))

not_rejected_log_single <- function(y,
                             alpha_log = alpha_log_gauss,
                             log_q_over_p = log_q_over_p_single)
    log(runif(1)) <= alpha_log + log_q_over_p(y)

r_abs_norm <- function(n)
    abs(rnorm(n))

sim_loop <- function(n,
                     r_proposal = r_abs_norm,
                     not_rejected = not_rejected_log_single)
{
    y <- numeric(n)
    count <- 1
    while (count <= n) {
        y_proposal <- r_proposal(1)
        if (not_rejected(y_proposal)) {
            y[count] <- y_proposal
            count <- count + 1
        }
    }
    y
}

log_q_over_p_more <- Vectorize(log_q_over_p_single)

not_rejected_more <- function(y,
                              alpha_log = alpha_log_gauss,
                              log_q_over_p = log_q_over_p_more)
    log(runif(length(y))) <= alpha_log + log_q_over_p(y)

## Tries to simulate n values by estimating the probability of rejection based
## on the total number of accepted and rejected proposals so far.
sim_approx <- function(n,
                       total_accepted,
                       total_rejected,
                       r_proposal = r_abs_norm,
                       not_rejected = not_rejected_more)
{
    if (total_accepted == 0)
        n_proposal <- n
    else
        n_proposal <- n * (total_accepted + total_rejected) / total_accepted
    y <- r_proposal(n_proposal)
    accepted <- y[not_rejected(y)]
    n_accepted <- length(accepted)
    list(
        accepted = accepted,
        n_accepted = n_accepted,
        n_rejected = n_proposal - n_accepted
        
    )
}

sim_smart <- function(n,
                      r_proposal = r_abs_norm,
                      not_rejected = not_rejected_more)
{
    y <- numeric(n)
    total_accepted <- 0
    total_rejected <- 0
    while (total_accepted < n) {
        sim <- sim_approx(n,
                          total_accepted,
                          total_rejected,
                          r_proposal,
                          not_rejected)
        if (sim$n_accepted > 0) {
            y[(1 + total_accepted):(total_accepted + sim$n_accepted)] <-
                sim$accepted
            total_accepted <- total_accepted + sim$n_accepted
        }
        total_rejected <- total_rejected + sim$n_rejected
    }
    y[1:n]
}

### Adaptive envelopes

logf <- Vectorize(function(y)
    y * sum_xz - sum(exp(y * x)))

dlogf <- Vectorize(function(y)
    sum_xz - sum(x * exp(y * x)))

## For the following, the index of Q and z are shifted by 1 relative to the
## notes. So my Q[2] corresponds to Q[1] from the notes (since I need Q[1] to
## correspond to Q[0] from the notes)
    
get_accepted_proposal <- function(a, b, c, z, Q)
{
    reject <- TRUE
      while(reject) {
          cu <- c * runif(1)
          i <- findInterval(cu, Q)
          proposal <- log(a[i] * exp(-b[i]) * (cu - Q[i]) + exp(a[i] * z[i])) / a[i]
          reject <- log(runif(1)) > logf(proposal) - a[i] * proposal - b[i]
      }
    proposal
}

sim_adapt_two <- function(n, x1, x2)
{
    a <- dlogf(c(x1, x2))
    b <- logf(c(x1,x2)) - a * c(x1,x2)
    z <- c(0, (b[2] - b[1]) / (a[1] - a[2]), 1)
    Q <- c(0, exp(b[1]) * (exp(a[1] * z[2]) - 1) / a[1])
    c <- Q[2] + exp(b[2]) * (exp(a[2] * z[3]) - exp(a[2] * z[2])) / a[2]
    replicate(n, get_accepted_proposal(a, b, c, z, Q))
}

### Adaptive envelopes with more than 2 points on the grid

## I use t instead of x, since x is already used by the data set.

get_z <- function(a, b)
{
    N <- length(a)
    c(
        0,
        (b[2:N] - b[1:(N-1)]) / (a[1:(N-1)] - a[2:N]),
        1
    )
}

get_Q_and_c <- function(a, b, z)
{
    N <- length(a)
    Qc <- cumsum(c(0,
                   exp(b[1:N]) * (exp(a[1:N] * z[2:(N+1)])
                       - exp(a[1:N] * z[1:N])) / a[1:N]))
    list(
        Q = Qc[1:N],
        c = Qc[N+1]
    )
}

sim_adapt_more <- function(n, t)
{
    a <- dlogf(t)
    b <- logf(t) - a * t
    z <- get_z(a, b)
    Q_and_c <- get_Q_and_c(a, b, z)
    replicate(n, get_accepted_proposal(a, b, Q_and_c$c, z, Q_and_c$Q))
}

### More generic solution

## The user has to supply logf.
## This only applies to get_accepted_proposal and sim_adapt_more.
get_accepted_proposal_generic <- function(a, b, c, z, Q, logf)
{
    reject <- TRUE
    while(reject) {
        cu <- c * runif(1)
        i <- findInterval(cu, Q)
        proposal <- log(a[i] * exp(-b[i]) * (cu - Q[i]) + exp(a[i] * z[i])) / a[i]
        reject <- log(runif(1)) > logf(proposal) - a[i] * proposal - b[i]
    }
    proposal
}

## We make dlogf an optional argument.
## If it is NULL, then we use numerical differentiation.
sim_adapt_more_generic <- function(n, t, logf, dlogf = NULL)
{
    if (is.null(dlogf))
        a <- grad(logf, t, method = "simple")
    else
        a <- dlogf(t)
    b <- logf(t) - a * t
    z <- get_z(a, b)
    Q_and_c <- get_Q_and_c(a, b, z)
    replicate(n, get_accepted_proposal(a, b, Q_and_c$c, z, Q_and_c$Q, logf))
}

