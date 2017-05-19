# Maximum likelihood fit to a Poisson distribution
fit_p <- function(params, y) {
  -sum(log(dpois(y, lambda = params)))
}

# Maximum likelihood fit to a Negative Binomial distribution
fit_nb <- function(params, y) {
  -sum(log(dnbinom(y, mu = params[1], size = params[2])))
}

# Maximum likelihood fit to a Poisson Inverse Gaussian distribution
fit_pig <- function(params, y) {
  -sum(log(
    gamlss.dist::dSICHEL(y, mu = params[1], sigma = params[2], nu = -0.5)
  ))
}

# Maximum likelihood fit to a Sichel distribution
fit_sichel <- function(params, y) {
  -sum(log(
    gamlss.dist::dSICHEL(y, mu = params[1], sigma = params[2], nu = params[3])
  ))
}

# get frequencies of count data
get_frequency <- function(x) {
  x <- hist(x, breaks = seq(-0.5, 40.5, by = 1), plot = FALSE)
  x <- x$counts / sum(x$counts)
  x[which(x == 0)] <- NA
  return(x)
}

# calculate AIC
AIC_calc <- function(x) {
  2 * x$value + 2 * length(x$par)
}

# make AIC table
AIC_table <- function(models, aics) {
  out      <- data.frame(Model = models, AIC = aics)
  out      <- out[order(out$AIC), ]
  out$dAIC <- out$AIC - out$AIC[1]
  return(out)
}

# plot distribution fits
plot_distribution_fits <- function(fits, obs, title) {
  p <- ggplot2::ggplot(data = fits, aes(x = x, y = y)) +
    theme_classic() +
    ggtitle(title) +
    scale_y_continuous(name = "Frequency", expand = c(0, 0)) +
    scale_x_continuous(name = "Distance (in patches)", expand = c(0, 0)) +
    geom_line(aes(col = d, linetype = d), size = 1) +
    geom_point(data = na.omit(data.frame(x = obs$x, y = obs$frequency)),
               aes(x = x, y = y), pch = 19, size = 1.4) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "dashed"),
                          guide = guide_legend(title = "Distribution")) +
    scale_color_manual(values = c('#4daf4a', '#e41a1c', '#984ea3', '#377eb8'),
                       guide = guide_legend(title = "Distribution")) +
    theme(legend.position = c(0.8, 0.8))

  return(p)
}

#' Test multiple kernel fits to dispersal data.
#'
#' \code{test_kernel_fits} takes output from \code{\link{wrangle_beetle_data}}
#'    and runs tests to determine the best-fitting dispersal kernel for all
#'    beetles as a whole, for females independently, and for males
#'    independently. The following kernels are tested: the Poisson distribution
#'    (\code{\link[stats]{Poisson}}), the negative binomial distribution
#'    (\code{\link[stats]{NegBinomial}}), the Sichel distribution
#'    (\code{\link[gamlss.dist]{SICHEL}}), and the Poisson Inverse Gaussian
#'    distribution (a special case of the Sichel distribution where
#'    \eqn{\nu = -0.5}).
#'
#' @inheritParams fit_Bev_Holt
#' @return A list of lists. Each top-level list (\code{$all}, \code{$females},
#'     or \code{$males}) contains the following objects:
#'     \itemize{
#'         \item{\code{$dist}}{ - A vector of the distances that were measured
#'             for each individual in the data subset.}
#'         \item{\code{$x}}{ - A vector of patch distances relevant to the study
#'             (\code{x = 0:40}).}
#'         \item{\code{$AIC}}{ - An AIC table ranking each distribution by its
#'             AIC score when fit to the data in \code{$dist}.}
#'         \item{\code{$frequency}}{ - The frequency of each distance measured in
#'             $dist, measured over the patches in \code{$x}.}
#'         \item{\code{$plot}}{ - A \code{\link[ggplot2]{ggplot}} object showing
#'             the data in \code{$frequency} along with the best-fit plot of
#'             each candidate distribution.}
#'         \item{\code{$<distribution>}}{ - The results of maximum-likelihood
#'             fits of each distribution to the data in $dist, where
#'             \code{<distribution>} is either \code{poisson}, \code{nbinom},
#'             \code{pig}, or \code{sichel}. The structure of these results
#'             exaclty follows the output from \code{\link[stats]{optim}}.
#'             For example, to get parameter estimates for fits to the Poisson
#'             distribution, using data from female beetles, use
#'             \code{output$females$poisson$par}.}
#'         \item{\code{$<distribution>$fit}}{ - The mass density function for a
#'             given distribution over \code{$x} (i.e., the line drawn in
#'             \code{$plot}.}
#'     }
#' @export
test_kernel_fits <- function(clean_data) {
  # SUBSET DATA ------------------------------------------------------------------
  # all the data, males only, and females only
  all     <- list(dist = as.numeric(abs(na.omit(
    clean_data$dist
    ))))
  females <- list(dist = as.numeric(abs(na.omit(
    clean_data[clean_data$sex == 'f', ]$dist
    ))))
  males   <- list(dist = as.numeric(abs(na.omit(
    clean_data[clean_data$sex == 'm', ]$dist
    ))))

  # FREQUENCY OF EACH DISTANCE -------------------------------------------------
  all$frequency     <- get_frequency(all$dist)
  females$frequency <- get_frequency(females$dist)
  males$frequency   <- get_frequency(males$dist)

  # MAXIMUM-LIKELIHOOD FITS ----------------------------------------------------
  # Poisson distribution
  all$poisson     <- optim(par = 5, fn = fit_p, y = all$dist,
                           method = "Br", lower = 0, upper = 10)
  females$poisson <- optim(par = 5, fn = fit_p, y = females$dist,
                           method = "Br", lower = 0, upper = 10)
  males$poisson   <- optim(par = 5, fn = fit_p, y = males$dist,
                           method = "Br", lower = 0, upper = 10)

  # Negative Binomial distribution
  all$nbinom     <- optim(par = c(5, 1), fn = fit_nb, y = all$dist)
  females$nbinom <- optim(par = c(5, 1), fn = fit_nb, y = females$dist)
  males$nbinom   <- optim(par = c(5, 1), fn = fit_nb, y = males$dist)

  # Poisson Inverse Gaussian distribution
  all$pig     <- optim(par = c(5, 0.1), fn = fit_pig, y = all$dist)
  females$pig <- optim(par = c(5, 0.1), fn = fit_pig, y = females$dist)
  males$pig   <- optim(par = c(5, 0.1), fn = fit_pig, y = males$dist)

  # Sichel
  all$sichel     <- optim(par = c(5, 1, 1), fn = fit_sichel, y = all$dist)
  females$sichel <- optim(par = c(5, 1, 1), fn = fit_sichel, y = females$dist)
  males$sichel   <- optim(par = c(5, 1, 1), fn = fit_sichel, y = males$dist)

  # CALCULATE DENSITY ----------------------------------------------------------
  # patch number
  x <- 0:40

  # Poisson distribution
  all$poisson$fit     <- dpois(x = x,
                               lambda = all$poisson$par[1])
  females$poisson$fit <- dpois(x = x,
                               lambda = females$poisson$par[1])
  males$poisson$fit   <- dpois(x = x,
                               lambda = males$poisson$par[1])

  # Negative Binomial distribution
  all$nbinom$fit     <- dnbinom(x = x,
                                mu   = all$nbinom$par[1],
                                size = all$nbinom$par[2])
  females$nbinom$fit <- dnbinom(x = x,
                                mu   = females$nbinom$par[1],
                                size = females$nbinom$par[2])
  males$nbinom$fit   <- dnbinom(x = x,
                                mu   = males$nbinom$par[1],
                                size = males$nbinom$par[2])

  # Poisson Inverse Gaussian distribution
  all$pig$fit     <- gamlss.dist::dSICHEL(x = x,
                                          mu    = all$pig$par[1],
                                          sigma = all$pig$par[2],
                                          nu    = -0.5)
  females$pig$fit <- gamlss.dist::dSICHEL(x = x,
                                          mu    = females$pig$par[1],
                                          sigma = females$pig$par[2],
                                          nu    = -0.5)
  males$pig$fit   <- gamlss.dist::dSICHEL(x = x,
                                          mu    = males$pig$par[1],
                                          sigma = males$pig$par[2],
                                          nu    = -0.5)

  # Sichel distribution
  all$sichel$fit     <- gamlss.dist::dSICHEL(x = x,
                                             mu    = all$sichel$par[1],
                                             sigma = all$sichel$par[2],
                                             nu    = all$sichel$par[3])
  females$sichel$fit <- gamlss.dist::dSICHEL(x = x,
                                             mu    = females$sichel$par[1],
                                             sigma = females$sichel$par[2],
                                             nu    = females$sichel$par[3])
  males$sichel$fit   <- gamlss.dist::dSICHEL(x = x,
                                             mu    = males$sichel$par[1],
                                             sigma = males$sichel$par[2],
                                             nu    = males$sichel$par[3])


  # PLOT DATA ------------------------------------------------------------------
  # make data frames for plotting
  a_plot_df <- rbind(data.frame(x = x, y = all$poisson$fit, d = "Poisson"),
                     data.frame(x = x, y = all$nbinom$fit,  d = "NBinom"),
                     data.frame(x = x, y = all$pig$fit,     d = "PIG"),
                     data.frame(x = x, y = all$sichel$fit,  d = "Sichel"))


  f_plot_df <- rbind(data.frame(x = x, y = females$poisson$fit, d = "Poisson"),
                     data.frame(x = x, y = females$nbinom$fit,  d = "NBinom"),
                     data.frame(x = x, y = females$pig$fit,     d = "PIG"),
                     data.frame(x = x, y = females$sichel$fit,  d = "Sichel"))


  m_plot_df <- rbind(data.frame(x = x, y = males$poisson$fit, d = "Poisson"),
                     data.frame(x = x, y = males$nbinom$fit,  d = "NBinom"),
                     data.frame(x = x, y = males$pig$fit,     d = "PIG"),
                     data.frame(x = x, y = males$sichel$fit,  d = "Sichel"))

  # add patch number to output lists
  all$x     <- x
  females$x <- x
  males$x   <- x

  # plot data and fits
  all$plot     <- plot_distribution_fits(a_plot_df, all,     "All Beetle Data")
  females$plot <- plot_distribution_fits(f_plot_df, females, "Female Data Only")
  males$plot   <- plot_distribution_fits(m_plot_df, males,   "Male Data Only")


  # MODEL SELECTION ------------------------------------------------------------
  all$AIC     <- AIC_table(models = c("Poisson",
                                      "Negative Binomial",
                                      "Poisson Inverse Gaussian",
                                      "Sichel"),
                           aics = c(AIC_calc(all$poisson),
                                    AIC_calc(all$nbinom),
                                    AIC_calc(all$pig),
                                    AIC_calc(all$sichel)))

  females$AIC <- AIC_table(models = c("Poisson",
                                      "Negative Binomial",
                                      "Poisson Inverse Gaussian",
                                      "Sichel"),
                           aics = c(AIC_calc(females$poisson),
                                    AIC_calc(females$nbinom),
                                    AIC_calc(females$pig),
                                    AIC_calc(females$sichel)))

  males$AIC   <- AIC_table(models = c("Poisson",
                                      "Negative Binomial",
                                      "Poisson Inverse Gaussian",
                                      "Sichel"),
                           aics = c(AIC_calc(males$poisson),
                                    AIC_calc(males$nbinom),
                                    AIC_calc(males$pig),
                                    AIC_calc(males$sichel)))

  return(list(all = all, females = females, males = males))
}
