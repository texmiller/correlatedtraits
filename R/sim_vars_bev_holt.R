if(F) {
  library(correlatedtraits)
  library(MCMCglmm)
  library(rstan)
  library(mvtnorm)

  rm(list=ls())
  #set.seed(42)

  # this code (mostly) from this post:
  #http://discourse.mc-stan.org/t/multivariate-animal-model/631/13

  # number of individuals to simulate
  n_dam <- 3
  n_sir <- 28

  # number of parents
  n_par <- n_sir + n_dam * n_sir

  # Pedigrees
  # individuals where both traits were measured
  ped <- sim_halfsib(n_sir, n_dam, 8)

  # combine pedigrees
  ped$animal <- 1:nrow(ped)

  # G matrix
  G <- matrix(c(0.2, 0.05, 0.05, 0.1), 2, 2)

  # Residual matrix
  E <- matrix(c(0.4, 0.10, 0.10, 0.2), 2, 2)

  # Simulated breeding values
  a <- rbv(ped, G)

  # Fixed effects
  beta <- matrix(c(1.8, log(30)), 1, 2)

  Intercept = rep(1, nrow(a))
  X <- cbind(Intercept)

  # Correlated noise
  e <- rmvnorm(nrow(a), sigma = E)

  # Simulated data
  Y <- X %*% beta + a + e

  den <- c(rep(NA, n_par), rep(rep(1/c(1, 3, 5, 10), each = 2), times = n_dam * n_sir))
  D_exp <- exp(Y[, 1])
  r_exp <- exp(Y[, 2])
  b     <- 10

  #simple negative binomial
  D_obs <- rpois(nrow(Y), D_exp)
  H_obs <- rnbinom(nrow(Y), mu = r_exp / (1 + b * den), size = 10)

  # censor traits that weren't measured
  D_obs[1:n_par] <- NA
  H_obs[1:n_par] <- NA
  #D_obs[(nrow(ped_Dr) + nrow(ped_D) + 1):nrow(ped)] <- NA
  #r_obs[(nrow(ped_Dr) + 1):(nrow(ped_Dr) + nrow(ped_D))] <- NA

  D_exp[1:n_par] <- NA
  r_exp[1:n_par] <- NA
  #D_exp[(nrow(ped_Dr) + nrow(ped_D) + 1):nrow(ped)] <- NA
  #r_exp[(nrow(ped_Dr) + 1):(nrow(ped_Dr) + nrow(ped_D))] <- NA

  dat_i = cbind(cumsum(!is.na(D_obs)), cumsum(!is.na(H_obs)))
  dat_i[is.na(D_obs), 1] <- 0
  dat_i[is.na(H_obs), 2] <- 0

  # stan does not allow NAs; remove them
  D_obs <- na.omit(D_obs)
  H_obs <- na.omit(H_obs)
  den   <- na.omit(den)

  D_exp <- na.omit(D_exp)
  r_exp <- na.omit(r_exp)

  # Create relationship matrix
  inv.phylo <- MCMCglmm::inverseA(ped, scale = TRUE)
  A <- solve(inv.phylo$Ainv)
  A <- (A + t(A))/2 # Not always symmetric after inversion
  rownames(A) <- rownames(inv.phylo$Ainv)
  isSymmetric((A))

  stan_data = list(M   = ncol(Y),
                   N   = nrow(Y),
                   ND  = length(D_obs),
                   NH  = length(H_obs),
                   D   = c(D_obs),
                   H   = c(H_obs),
                   iD  = which(dat_i[, 1]>0),
                   iH  = which(dat_i[, 2]>0),
                   den = den,
                   A   = as.matrix(A)
  )

  # set parallel options
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  fits <- rstan::stan(file     = system.file("multivariate_animal_model_bev_holt.stan",
                                             package = "correlatedtraits"),
                      #init     = inits,
                      data     = stan_data,
                      verbose  = TRUE,
                      refresh  = 100,
                      chains   = 3,
                      control  = list(adapt_delta   = 0.99))

  model = rstan::extract(fits)
  rstan::traceplot(fits, pars = "G")
  rstan::traceplot(fits, pars = "E")
  rstan::traceplot(fits, pars = "P")
  rstan::summary(fits, pars = "G")$summary; colMeans(model$G); G
  rstan::summary(fits, pars = "E")$summary; colMeans(model$E); E
  rstan::summary(fits, pars = "P")$summary; colMeans(model$P); G + E

  hist(colMeans(model$mu_D) - log(D_exp), breaks = seq(-10, 10, by = 0.2))
  hist(colMeans(model$mu_r) - log(r_exp), breaks = seq(-10, 10, by = 0.2))

  summary(fits, pars = "beta")$summary; t(colMeans(model$beta)); beta

  summary(fits, pars = "b")$summary; 10
}
