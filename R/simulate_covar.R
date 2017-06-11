if(F) {
  rho <- 0
  var <- c(5, 15)
  mu  <- c(11.3, 36.1)
  N   <- 8  # number of individuals per group
  G   <- 20 # number of groups

  # given parameter correlation (rho) and variance of each variable, calculate
  # the variance-covariance matrix
  simulated_Sigma <- sqrt(var) %*% t(sqrt(var)) * matrix(c(1, rho, rho, 1),
                                                         nrow = 2)
  # sample parameters
  mus <- MASS::mvrnorm(G, mu, simulated_Sigma, empirical = TRUE)
  y <- matrix(rep(mus, each = N), ncol = 2)

  # draw random population densities between 0 and K
  den <- rep(rep(c(1, 2, 4/3, 10), each=N/4), each = G)

  # simulate data
  off_lo <- rnorm(nrow(y), mean = y[,1], sd = 0.1)[den == 1]
  off_hi <- rnorm(nrow(y), mean = y[,2], sd = 0.1)[den == 10]

  # matrix of zeros for means of multivariate normal
  zeros <- matrix(0, nrow = G, ncol = 2)

  family <- rep(1:G, each = N/4)

  data <- list(N = N,
               G = G,
               family = family,
               lo_den_off   = off_lo,
               hi_den_off   = off_hi,
               zeros = zeros)

  # set parallel options
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  inits <- function() {
    list(beta_A = 11,  beta_B = 36)
  }

  fits <- rstan::stan(file     = system.file("simulate_covar_fit.stan",
                                             package = "correlatedtraits"),
                      #init     = inits,
                      data     = data,
                      verbose  = TRUE,
                      refresh  = 500,
                      chains   = 1,
                      save_dso = TRUE)

  (quantile(rstan::extract(fits, pars = c('a'), permuted = TRUE)$a, 0.5))
  (quantile(rstan::extract(fits, pars = c('est_B'), permuted = TRUE)$est_B, 0.5))
  (quantile(rstan::extract(fits, pars = c('sda'), permuted = TRUE)$sda, 0.5))
  (quantile(rstan::extract(fits, pars = c('sdb'), permuted = TRUE)$sdb, 0.5))
  (quantile(rstan::extract(fits, pars = c('sigmaa'), permuted = TRUE)$sigmaa, 0.5)^2)
  (quantile(rstan::extract(fits, pars = c('sigmab'), permuted = TRUE)$sigmab, 0.5)^2)
  cov(y)
  cov(cbind(off_lo, off_hi))

  out <- rstan::extract(fits, pars = c('S', 'est_A', 'est_B'), permuted = TRUE)

  params <- data.frame(actual = c(mu[1], mu[2]),
                       '2.5%'  = c(quantile(out$est_A, 0.025),
                                   quantile(out$est_B, 0.025)),
                       '50%'   = c(quantile(out$est_A, 0.5),
                                   quantile(out$est_B, 0.5)),
                       '97.5%' = c(quantile(out$est_A, 0.975),
                                   quantile(out$est_B, 0.975)))

  rownames(params) <- c("A", "B")
  colnames(params) <- c("actual", "2.5%", "50%", "97.5%")
  params

  mat_names <- c("A", "B")

  colnames(simulated_Sigma) <- mat_names
  rownames(simulated_Sigma) <- mat_names


  estimated_Sigma <- matrix(c(quantile(out$S[, 1, 1], 0.025),
                              quantile(out$S[, 1, 2], 0.025),
                              quantile(out$S[, 2, 1], 0.025),
                              quantile(out$S[, 2, 2], 0.025)),
                            nrow = 2)

  simulated_Rho <- data.frame(matrix(c(1, rho, rho, 1), nrow = 2))
  colnames(simulated_Rho) <- mat_names
  rownames(simulated_Rho) <- mat_names

  estimated_Rho <- data.frame(cov2cor(estimated_Sigma))
  colnames(estimated_Rho) <- mat_names
  rownames(estimated_Rho) <- mat_names

  estimated_Sigma <- data.frame(estimated_Sigma)
  colnames(estimated_Sigma) <- mat_names
  rownames(estimated_Sigma) <- mat_names

  (ret <- list(parameters = params,
              Sigma      = list(simulated = simulated_Sigma,
                                estimated = estimated_Sigma),
              Rho        = list(simulated = simulated_Rho,
                                estimated = estimated_Rho)))
  cov(mus)
  cov(y)
  cov(cbind(off_lo, off_hi))
  (quantile(rstan::extract(fits, pars = c('sigmab'), permuted = TRUE)$sigma, 0.975))

}
