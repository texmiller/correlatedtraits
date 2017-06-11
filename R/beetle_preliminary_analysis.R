#' Fit raw beetle data to a specified growth function..
#'
#' This code takes data returned from \code{\link{wrangle_beetle_data}}, fits a
#'     user-specified growth function, plots the resulting curve, and returns
#'     relevant parameter values as well as a transformed dataset.
#'     \code{fit_growth} transforms the data -- which has constant female
#'     density but variable bean density -- to have constant bean density and
#'     variable female density.
#'
#' The form of the Beverton-Holt equation is taken from Otto and Day (2007):
#'     \deqn{
#'         N_{t+1} = \frac{R_{0} N_t}{1 + \alpha N_t}
#'     }{
#'         N(t+1) = R_0 * N(t)/ (1 + \alpha N(t))
#'     }
#'
#' The logistic equation uses the discrete-time logistic formula:
#'     \deqn{
#'         N_{t+1} = N_{t} + N_{t}r_{d} (1 - \frac{N_t}{K})
#'     }{
#'         N(t+1) = N(t) + N(t)r * (1 - N(t)/K)
#'     }
#'
#' Data were originally collected by varying resource density (i.e., number of
#'     beans) and keeping female density constant (one female per trial). Thus,
#'     fecundity was measured at the following female-to-bean ratios:
#'     \deqn{
#'         \frac{1}{1}\space,
#'         \frac{1}{3}\space,
#'         \frac{1}{5}\space,
#'         \frac{1}{10}\space
#'     }{
#'         1/1, 1/3, 1/5, 1/10
#'     }
#'     Assuming the relationship that drives density dependence is this ratio of
#'     females-to-beans (and not the absolute number of females or beans), we
#'     can re-scale these fractions to yield a common denominator (bean density)
#'     among all fecundity trials. For female density \eqn{F} and bean
#'     density \eqn{B}, the fraction \eqn{\frac{F}{B}}{F/B} can be re-scaled to
#'     have any denominator by multiplying by some proportion equal to 1:
#'     \deqn{
#'         \frac{\rho}{\rho} \times \frac{F}{B} = \frac{F}{B}
#'     }{
#'         (\rho F)/(\rho B) = F/B
#'     }
#'     We can choose a proportion \eqn{\rho} such that the observed bean density
#'     \eqn{B} is scaled to some desired bean density, \eqn{B_{new}}{B_new}:
#'     \deqn{
#'         \rho = \frac{\frac{B_{new}}{B}}{\frac{B_{new}}{B}} = 1
#'     }{
#'         \rho = (B_new/B) / (B_new/B) = 1
#'     }
#'     For example, when \eqn{B_{new} = 10}{B_new = 10} and \eqn{B = 5},
#'     \eqn{\rho = 2/2}.
#'
#'     In this analysis \eqn{B_{new} = 10}{B_new = 10}, which gives the
#'     following female-to-bean ratios at constant resource densities:
#'     \deqn{
#'         \frac{10}{10}\space,
#'         \frac{3. \overline{33}}{10}\space,
#'         \frac{2}{10}\space,
#'         \frac{1}{10}\space
#'     }{
#'         10/10, 3.33/10, 2/10, 1/10
#'     }
#'
#' @param data The cleaned and compiled data returned from
#'     \code{\link{wrangle_beetle_data}}.
#' @param show_plot Should the plot be drawn? Default is \code{TRUE}.
#' @param xmax The maximum x-axis value for the plot.
#' @param model Which model should be used to fit the growth data? Currenly,
#'     either \code{Beverton-Holt} or \code{logistic}.
#' @return A plot of the cleaned data and the fitted growth function, the
#'     result of the model fit, parameter estimates for the growth rate,
#'      carrying capacity, and any other model parameters, and the transformed
#'      data set.
#'
#' @references Otto and Day (2007). A Biologist's Guide to Mathematical Modeling
#'     in Ecology and Evolution. Princeton University Press. Page 185.
#' @references Gotelli (2001). A Primer of Ecology, Third Edition. Sinauer Associates,
#'     Inc. Page 35.
#' @export
fit_growth <- function(data, show_plot = TRUE, xmax = 30, model) {

  # remove NA rows from data
  tmp <- na.omit(data)

  # Modify the fraction of females per bean so that all resource densities are
  # constant.
  data <- data.frame(x = 10 / tmp$beans, y = tmp$f * 10 / tmp$beans)

  if(grepl(model, "Beverton-Holt", ignore.case = TRUE)){
    # fit a Beverton-Holt model
    # The alpha (a) in the equation below is the inverse (1/alpha) from Otto and
    # Day (2007), because the model won't converge otherwise.

    # get starter values for `nls` using the log of y; this function is a bit
    # finicky and occasionally won't run otherwise.
    start <- nls(log(y) ~ log((R0 * x) / (1 + x / alpha)), data,
                 list(R0 = 1, alpha = 1))

    # use starter values to get estimates for untransformed data.
    fit <- nls(y ~ (R0 * x) / (1 + x / alpha), data, coef(start))

    # Get coefficients and calculate K
    a <- coef(fit)[2]
    params <- list(R0     = coef(fit)[1],
                   K      = coef(fit)[2] * (coef(fit)[1] - 1),
                   alpha  = coef(fit)[2])
  } else if(grepl(model, "logistic", ignore.case = TRUE)){

    # fit a discrete-time logistic model
    fit <- nls(y ~ K * x * exp(r) / (K + x * (exp(r) - 1)),
                data, list(r = 2.42, K = 36.1))

    # Get coefficients
    params <- list(r = coef(fit)[1],
                   K = coef(fit)[2])

  } else if(grepl(model, "monomolecular", ignore.case = TRUE)){

    fit <- nls(y ~ K * (1 - exp(-r * x)),
               data, list(r = 1, K = 36.1))

    # Get coefficients
    params <- list(r = coef(fit)[1],
                   K = coef(fit)[2])

  } else if(grepl(model, "Ricker", ignore.case = TRUE)){

    fit <- nls(y ~ x * exp(r * (1 - x / K)),
               data, list(r = 1, K = 36.1))

    # Get coefficients
    params <- list(r = coef(fit)[1],
                   K = coef(fit)[2])
  } else if(grepl(model, "r-k Skellam", ignore.case = TRUE)){

    fit <- nls(y ~  K * (1 - (1 - r / K)^x),
               data, list(r = 11, K = 36))

    # Get coefficients
    params <- list(r = coef(fit)[1],
                   K = coef(fit)[2])

  } else {
    stop("The specified model is not available. See ?fit_growth")
  }

  # calculate AIC
  AIC <- 2 * length(coef(fit)) + nrow(data) * log(fit$m$deviance())

  # format parameters for plotting
  r <- as.numeric(format(params[[1]], digits = 3))
  K <- as.numeric(format(params[[2]], digits = 3))

  # make predictions on a range of values from 0:100 individuals
  data_predict <- data.frame(x = 0:100,
                             y = predict(fit,
                                         newdata = data.frame(x = 0:100)))

  # plot results
  p <-
    ggplot2::ggplot(
      data = data,
      ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_classic() +
    ggplot2::geom_jitter(
      width = 0.2,
      height = 1,
      alpha = 0.5,
      pch = 19,
      colour = "grey50") +
    ggplot2::geom_line(
      data = data_predict[1:(xmax + 1), ],
      ggplot2::aes(x, y),
      color = 'red',
      size = 1) +
    ggplot2::geom_hline(
      yintercept = K,
      linetype = "dashed") +
    ggplot2::geom_abline(
      intercept = 0,
      slope = mean(data[data$x == 1, ]$y),
      linetype = "dashed") +
    ggplot2::annotate(
      "text",
      x = 3,
      y = 90,
      size = 3,
      label = paste0("italic(r) == ", r),
      parse = TRUE) +
    ggplot2::scale_y_continuous(
      name   = expression(italic(N[t+1])),
      limits = c(0, 105),
      expand = c(0, 0),
      breaks = c(seq(0, 100, by = 25), K),
      labels = c(seq(0, 100, by = 25),
      bquote(italic("K ") ~ "=" ~ .(K)))) +
    ggplot2::scale_x_continuous(
      name   = expression(italic(N[t])),
      limits = c(0, xmax),
      expand = c(0, 0) ,
      breaks = c(0, 1, 2,  3.33,  seq(10, xmax, length = 3)),
      labels = c(0, 1, 2, "3.33", seq(10, xmax, length = 3))) +
    ggplot2::stat_summary(
      fun.y = mean,
      geom = "point",
      size = 2)

  if(show_plot) {
    print(p)
  }

  return(list(plot   = p,
              fit    = fit,
              AIC    = AIC,
              params = params,
              data   = data)
  )
}
