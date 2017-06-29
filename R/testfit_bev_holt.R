if(F) {
  library(correlatedtraits)

  rm(list=ls())

  set.seed(42)

  # set data directory
  dir <- "D:/Dropbox/1 Projects/Academic/2016 genetic covariance invasions/data/"

  # set paths to data files
  female_fd_path <- paste0(dir, "fecundity_dispersal_raw.csv")
  female_f_path  <- paste0(dir, "fecundity_no_dispersal_raw.csv")
  dates_path     <- paste0(dir, "experiment_dates.csv")
  male_d_path    <- paste0(dir, "male_dispersal.csv")

  # clean data
  data <-  wrangle_beetle_data(female_fd_path,
                               female_f_path,
                               male_d_path,
                               dates_path)

  data <- data[data$sex == "f", ]
  data <- data[data$t > 0, ]

  # Pedigrees
  # individuals where both traits were measured
  ped_DH <- data[!is.na(data$dist), ]

  # individuals where only low density growth was measured
  ped_H  <- data[is.na(data$dist), ]

  # build pedigree
  ped        <- rbind(ped_DH, ped_H)
  ped        <- ped[!is.na(ped$sire), ]
  ped$dam    <- (ped$sire - 1) * 3 + ped$dam + max(ped$sire)
  ped$animal <- NA
  ped        <- ped[, c("animal", "dam", "sire", "dist", "beans", "f")]

  parent_ids <- unique(ped[, c("dam", "sire")])
  parents <- data.frame(animal = NA, dam = parent_ids$dam, sire = parent_ids$sire, dist = NA, beans = NA, f = NA)
  parents$dam <- NA
  parents$sire <- NA

  ped <- rbind(parents, ped)
  ped$animal <- 1:nrow(ped)

  # Create relationship matrix
  inv.phylo <- MCMCglmm::inverseA(ped[, 1:3], scale = TRUE)
  A <- solve(inv.phylo$Ainv)
  A <- (A + t(A))/2 # Not always symmetric after inversion
  rownames(A) <- rownames(inv.phylo$Ainv)
  isSymmetric((A))

  # number of parents
  n_par <- nrow(parent_ids)

  # observations
  D_obs <- ped$dist
  H_obs <- ped$f
  den   <- c(1/na.omit(ped$beans))

  # censor traits that weren't measured
  dat_i = cbind(cumsum(!is.na(D_obs)), cumsum(!is.na(H_obs)))
  dat_i[is.na(D_obs), 1] <- 0
  dat_i[is.na(H_obs), 2] <- 0

  # stan does not allow NAs; remove them
  D_obs <- abs(na.omit(D_obs))
  H_obs <-   c(na.omit(H_obs))

  stan_data = list(M  = 2,
                   N   = nrow(ped),
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
                      data     = stan_data,
                      iter     = 10000,
                      verbose  = TRUE,
                      refresh  = 500,
                      chains   = 4,
                      control  = list(adapt_delta   = 0.99))
  save.image(file="bev_holt.rda")

  model = rstan::extract(fits)
  rstan::traceplot(fits, pars = "G")
  rstan::traceplot(fits, pars = "E")
  rstan::traceplot(fits, pars = "P")
  rstan::summary(fits, pars = "G")$summary; colMeans(model$G)
  rstan::summary(fits, pars = "E")$summary; colMeans(model$E)
  rstan::summary(fits, pars = "P")$summary; colMeans(model$P)

  rstan::traceplot(fits, pars = "corrG")
  rstan::traceplot(fits, pars = "corrE")
  rstan::summary(fits, pars = "corrG")$summary; colMeans(model$corrG)
  rstan::summary(fits, pars = "corrE")$summary; colMeans(model$corrE)

  rstan::traceplot(fits, pars = "h2")
  rstan::summary(fits, pars = "h2")$summary; colMeans(model$h)

  rstan::traceplot(fits, pars = "beta")
  rstan::summary(fits, pars = "beta")$summary

  rstan::traceplot(fits, pars = "b")
  rstan::summary(fits, pars = "b")$summary

  plot(D_obs, colMeans(model$sim_D), ylim=c(0, 25), xlim=c(0, 25))
  lines(c(-100,100),c(-100,100),'l',col='red')

  plot(H_obs, colMeans(model$sim_H), ylim=c(0, 40), xlim=c(0, 40))
  lines(c(-100,100),c(-100,100),'l',col='red')

  shinystan::launch_shinystan(fits)
}

