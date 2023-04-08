library(Rcpp)

sourceCpp("implementation.cpp")

sqrt2 <- sqrt(2)
sqrt2pi <- sqrt(2 * pi)

## Forces evaluation of all arguments
force_all_args <- function()
    as.list(parent.frame())

## Epanechnikov kernel in x
ep_kernel <- function(x)
{
    nonzero_index <- abs(x) <= 1
    result <- numeric(length(x))
    result[nonzero_index] <- (1 - x[nonzero_index]^2) * 3/4
    result
}

## Computes the kernel density estimate in a single point x based on observations
## x_obs and with bandwidth h
ep_density_single <- function(x, x_obs, h)
    sum(ep_kernel_cpp((x - x_obs) / h)) / (h * length(x_obs))

## Computes the kernel density estimate based on observations x_obs and with
## bandwidth h. The density is estimated in m points.
ep_density <- function(x_obs, h, m = 512)
{
    grid <- seq(min(x_obs) - 3 * h, max(x_obs) + 3 * h, length.out = m)
    list(
        x = grid,
        y = sapply(grid, function(x) ep_density_single(x, x_obs, h))
    )
}

ep_density_cpp_wrap <- function(x_obs, h, m = 512)
{
    grid <- seq(min(x_obs) - 3 * h, max(x_obs) + 3 * h, length.out = m)
    list(
        x = grid,
        y = ep_density_cpp(x_obs, h, grid)
    )
}

my_density <- function(x, bw = "plug-in", binning = TRUE, ...)
{
    if (bw == "plug-in")
        h <- bandwidth_plug_in_precalculated(x)
    else if (bw == "cv")
        h <- bandwidth_cv(x, ...)
    else if (bw == "silverman")
        h <- bandwidth_silverman(x)
    else if (is.function(bw))
        h <- bw(x, ...)
    else if (is.numeric(bw))
        h <- bw
    else
        stop("Invalid bw.")

    if (binning)
        xy <- ep_density_binning_cpp_wrap(x, h, ...)
    else
        xy <- ep_density_cpp_wrap(x, h, ...)

    structure(
        list(
            x = xy$x,
            y = xy$y,
            bw = h
        ),
        class = "my_density"
    )
}

plot.my_density <- function(my_density_object, ...)
    plot(my_density_object$x,
         my_density_object$y,
         type = "l",
         xlab = "x",
         ylab = "Density",
         main = paste("Density estimate with bandwidth",
                      round(my_density_object$bw, digits = 2)),
         ...)

## Silverman's estimate of the standard deviation.
## More robust to outliers than the regular sample standard deviation.
sd_silverman <- function(x_obs)
    min(sd(x_obs),
        IQR(x_obs) / 1.34)

## Silverman's rule of thumb
bandwidth_silverman <- function(x_obs)
    0.9 * sd_silverman(x_obs) * length(x_obs)^(-0.2)

#### AMISE plug-in bandwidth selection

## Fourth derivative of standard Gaussian density.
d4_std_gaussian_density <- function(x)
{
    x_sq <- x^2
    exp(-x_sq/2) * ((x_sq - 6) * x_sq + 3) / sqrt2pi
}

### diff_vec implementations to compare (used for estimate of the two-norm of
### f0'')

## Returns vector of all differences of pairs.
diff_vec_expand.grid <- function(x_obs)
{
    x_grid <- expand.grid(Xi = x_obs, Xj = x_obs)
    x_grid[, 1] - x_grid[, 2]
}

## Same as above, but uses "outer" instead of "expand.grid" and returns matrix.
diff_matrix_outer <- function(x_obs)
    outer(x_obs, x_obs, FUN = "-")

## Returns vector version of matrix result from above.
diff_vec_outer <- function(x_obs)
    as.vector(diff_matrix_outer(x_obs))

## Same as above, but uses for loop
diff_vec_loop <- function(x_obs)
{
    n <- length(x_obs)
    res <- numeric(n^2)
    count <- 1
    for (i in 1:n) {
        for (j in 1:n) {
            res[count] <- x_obs[i] - x_obs[j]
            count <- count + 1
        }
    }
    res
}

diff_half_outer <- function(x_obs)
{
    diff_mat <- diff_matrix_outer(x_obs)
    diff_mat[upper.tri(diff_mat)]
}

diff_half_combn <- function(x_obs)
    combn(x_obs, 2, FUN = diff)

## Estimate of the squared 2-norm of f0''.
d2f0_two_norm_squared_estimate <- function(x_obs,
                                           r = bandwidth_silverman(x_obs),
                                           diff_vec = diff_vec_expand.grid)
    sum(d4_std_gaussian_density(diff_vec(x_obs) / (sqrt2 * r))) /
        (length(x_obs)^2 * (sqrt2 * r)^5)

## Uses that the Gaussian density is symmetric and the upper triangle of the
## diff_matrix is equal to minus the lower triangle.
d2f0_two_norm_squared_estimate_half_matrix <- function(x_obs,
                                                       r = bandwidth_silverman(x_obs))
{
    diff_mat <- diff_matrix_outer(x_obs)
    rsqrt2 <- r * sqrt2
    (2 * sum(d4_std_gaussian_density(
            diff_mat[upper.tri(diff_mat)] / rsqrt2)) +
     length(x_obs) * d4_std_gaussian_0) / (length(x_obs)^2 * rsqrt2^5)
}

d2f0_two_norm_squared_estimate_combn <- function(x_obs,
                                                 r = bandwidth_silverman(x_obs))
{
    rsqrt2 <- r * sqrt2
    (2 * sum(d4_std_gaussian_density(combn(x_obs, 2, FUN = diff) / rsqrt2)) +
        length(x_obs) * d4_std_gaussian_0) / (length(x_obs)^2 * rsqrt2^5)
}

d2f0_two_norm_squared_estimate_cpp_wrap <- function(x_obs,
                                                    r  = bandwidth_silverman(x_obs))
    d2f0_two_norm_squared_estimate_cpp(x_obs, r)

d4_std_gaussian_0 <- d4_std_gaussian_density(0)

K_two_norm_squared <- 0.6
sigma4_K <- 0.2^2

## Calculates the plug-in estimate of the AMISE-optimal bandwidth.
bandwidth_plug_in <- function(x_obs)
{
    (K_two_norm_squared /
     (d2f0_two_norm_squared_estimate(x_obs) *
      sigma4_K *
      length(x_obs)))^0.2
}

## Same as above, but the expression is simplified.
bandwidth_plug_in_precalculated <- function(x_obs)
    (15 / (d2f0_two_norm_squared_estimate_cpp_wrap(x_obs) * length(x_obs)))^0.2

#### Cross validation bandwidth selection

## Returns function to optimize for cross validation
get_cv_objective_naive <- function(x_obs, k)
{
    force_all_args()
    n <- length(x_obs)
    ## All the indices with a 3 is I_3, and so on.
    fold_indices <- rep(1:k, times = ceiling(n / k))[sample(n)]

    ## Function to maximize
    function(h)
    {
        dens_estimates <- numeric(n)
        for (i in 1:n) {
            ## All x-values that are not in the same fold as the i'th
            ## observation.
            x_other_folds <- x_obs[fold_indices != fold_indices[i]]
            dens_estimates[i] <- ep_density_single(x_obs[i], x_other_folds, h)
        }
        ## Return log-likelihood.
        sum(log(dens_estimates))
    }
}


## Returns function to maximize for cross validation
## Implementation with caching
get_cv_objective <- function(x_obs, k)
{
    force_all_args()
    n <- length(x_obs)
    ## All the indices with a 3 is I_3, and so on.
    fold_indices <- rep(1:k, times = ceiling(n / k))[sample(n)]    
    other_folds_list <- vector("list", k)

    function(h)
    {
        dens_estimates <- numeric(n)
        for (i in 1:n) {
            ## All x-values that are not in the same fold as the i'th
            ## observation (I^{-i}).
            x_other_folds <- other_folds_list[[fold_indices[i]]]
            if (is.null(x_other_folds)) {
                x_other_folds <- x_obs[fold_indices != fold_indices[i]]
                other_folds_list[[fold_indices[i]]] <<- x_other_folds
            }
            dens_estimates[i] <- ep_density_single(x_obs[i],
                                                   x_other_folds,
                                                   h)
        }
        sum(log(dens_estimates))
    }
}
    

## Selects h by k-fold cross validation
bandwidth_cv <- function(x_obs, k = 10)
{
    optimize(
        f = get_cv_objective(x_obs, k),
        interval = c(0, max(x_obs) - min(x_obs)),
        maximum = TRUE
    )$maximum
}

#### Implementations used for benchmarking

logF12 <- log(read.table("../data/infrared.txt", header = TRUE)$F12)

## Fourth derivative of standard Gaussian density.
d4_std_gaussian_density_slow <- function(x)
    exp(-x^2/2) * (x^4 - 6 * x^2 + 3) / sqrt2pi

## Returns function to maximize for cross validation
bandwidth_cv_complete <- function(x_obs, k)
{
    n <- length(x_obs)
    ## All the indices with a 3 is I_3, and so on.
    fold_indices <- rep(1:k, times = ceiling(n / k))[sample(n)]    
    other_folds_list <- vector("list", k)
    for (i in 1:k)
        other_folds_list[[i]] <- x_obs[fold_indices != i]
        

    f <- function(h)
    {
        dens_estimates <- numeric(n)
        for (i in 1:n) {
            ## All x-values that are not in the same fold as the i'th
            ## observation (I^{-i}).            
            dens_estimates[i] <- ep_density_single(x_obs[i],
                                                   other_folds_list[[fold_indices[i]]],
                                                   h)
        }
        sum(log(dens_estimates))
    }

    optimize(
        f = f,
        interval = c(0, max(x_obs) - min(x_obs)),
        maximum = TRUE
    )$maximum
}


bandwidth_cv_complete_cache <- function(x_obs, k)
{
    n <- length(x_obs)
    ## All the indices with a 3 is I_3, and so on.
    fold_indices <- rep(1:k, times = ceiling(n / k))[sample(n)]    
    other_folds_list <- vector("list", k)

    f <- function(h)
    {
        dens_estimates <- numeric(n)
        for (i in 1:n) {
            ## All x-values that are not in the same fold as the i'th
            ## observation (I^{-i}).
            x_other_folds <- other_folds_list[[fold_indices[i]]]
            if (is.null(x_other_folds)) {
                print(i)
                x_other_folds <- x_obs[fold_indices != fold_indices[i]]
                other_folds_list[[fold_indices[i]]] <<- x_other_folds
            }
            dens_estimates[i] <- ep_density_single(x_obs[i],
                                                   x_other_folds,
                                                   h)
        }
        sum(log(dens_estimates))
    }

    optimize(
        f = f,
        interval = c(0, max(x_obs) - min(x_obs)),
        maximum = TRUE
    )$maximum
}

## Returns function to maximize for cross validation
bandwidth_cv_complete_naive <- function(x_obs, k)
{
    n <- length(x_obs)
    ## All the indices with a 3 is I_3, and so on.
    fold_indices <- rep(1:k, times = ceiling(n / k))[sample(n)]    

    f <- function(h)
    {
        dens_estimates <- numeric(n)
        for (i in 1:n) {
            ## All x-values that are not in the same fold as the i'th
            ## observation (I^{-i}).
            dens_estimates[i] <- ep_density_single(
                x_obs[i],
                x_obs[fold_indices != fold_indices[i]],
                      h)
        }
        sum(log(dens_estimates))
    }

    optimize(
        f = f,
        interval = c(0, max(x_obs) - min(x_obs)),
        maximum = TRUE
    )$maximum
}

## Returns function to maximize for cross validation
## Naive implementation
get_cv_objective_smart_but_no_cache <- function(x_obs, k)
{
    n <- length(x_obs)
    ## All the indices with a 3 is I_3, and so on.
    fold_indices <- rep(1:k, times = ceiling(n / k))[sample(n)]    
    other_folds_list <- vector("list", k)

    function(h)
    {
        dens_estimates <- numeric(n)
        for (i in 1:n) {
            ## All x-values that are not in the same fold as the i'th
            ## observation (I^{-i}).            
            x_other_folds <- other_folds_list[[fold_indices[i]]]
            if (is.null(x_other_folds)) {
                print(i)
                x_other_folds <- x_obs[fold_indices != fold_indices[i]]
                other_folds_list[[fold_indices[i]]] <- x_other_folds
            }
            dens_estimates[i] <- ep_density_single(x_obs[i],
                                                   x_other_folds,
                                                   h)
        }
        sum(log(dens_estimates))
    }
}

## Binning

bin_weights <- function(x, grid_lower, grid_length, grid_diff)
{
  w <- numeric(grid_length)
  for(i in seq_along(x)) {
    closest_grid_i <- ceiling((x[i] - grid_lower) / grid_diff + 0.5)
    w[closest_grid_i] <- w[closest_grid_i] + 1
  }
  w
}

ep_density_binning_cpp_wrap <- function(x, h, grid_length = 512)
{
    range_x <- range_cpp(x)
    grid <- seq(range_x[1] - 3 * h, range_x[2] + 3 * h, length.out = grid_length)
    list(
        x = grid,
        y = ep_density_binning_cpp(x, h, grid)
    )
}
    

ep_density_binning <- function(x, h, grid_length = 512)
{
    range_x <- range_cpp(x)
    grid_lower <- range_x[1] - 3 * h
    grid_upper <- range_x[2] + 3 * h
    grid <- seq(grid_lower, grid_upper, length.out = grid_length)
    grid_diff <- grid[2] - grid[1]
    max_nonzero_i <- floor(h / grid_diff)
    kernel_evals <- ep_kernel_cpp((grid[1:(max_nonzero_i + 1)] - grid_lower) / h)
    kernel_vec <- c(rev(kernel_evals[-1]), kernel_evals)
    weights <- c(rep(0, max_nonzero_i),
                 bin_weights_cpp(x, grid_lower, grid_length, grid_diff),
                 rep(0, max_nonzero_i))
    y <- numeric(grid_length)
    for (i in (1 + max_nonzero_i):(grid_length + max_nonzero_i))
        y[i - max_nonzero_i] <- sum(weights[(i - max_nonzero_i):
                                            (i + max_nonzero_i)] *
                                    kernel_vec)

    list(x = grid, y = y / (h * length(x)))
}

ep_density_binning_toeplitz <- function(x, h, grid_length = 512)
{
    n <- length(x)
    range_x <- range_cpp(x)
    grid_lower <- range_x[1] - 3 * h
    grid_upper <- range_x[2] + 3 * h    
    grid <- seq(grid_lower, grid_upper, length.out = grid_length)
    grid_diff <- grid[2] - grid[1]
    kernel_evals <- ep_kernel_cpp((grid - grid_lower) / h)
    weights <- bin_weights_cpp(x, grid_lower, grid_length, grid_diff)
    list(
        x = grid,
        y = colSums(weights * toeplitz(kernel_evals)) / (h * length(x))
    )
}

## Epanechnikov kernel in x
ep_kernel <- function(x)
{
    nonzero_index <- abs(x) <= 1
    result <- numeric(length(x))
    result[nonzero_index] <- (1 - x[nonzero_index]^2) * 3/4
    result
}
