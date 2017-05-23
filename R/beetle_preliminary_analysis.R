#' Fit raw beetle data to a Beverton-Holt recruitment curve.
#'
#' This code takes data returned from \code{\link{wrangle_beetle_data}}, fits a
#'     Beverton-Holt recruitement curve, plots the resulting curve, and returns
#'     relevant parameter values as well as a transformed dataset.
#'     \code{fit_Bev_Holt} transforms the data -- which has constant female
#'     density but variable bean density -- to have constant bean density and
#'     variable female density.
#'
#' The form of the Beverton-Holt equation is taken from Otto and Day (2007):
#'     \deqn{
#'         N_{t+1} = \frac{\lambda N_t}{1 + \alpha N_t}
#'     }{
#'         N(t+1) = \lambda N(t)/ (1 + \alpha N(t))
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
#' @param clean_data The cleaned and compiled data returned from
#'     \code{\link{wrangle_beetle_data}}.
#' @param show_plot Should the plot be drawn? Default is \code{TRUE}.
#' @param xmax The maximum x-axis value for the plot.
#' @return A plot of the cleaned data and the fitted Beverton-Holt function, the
#'     result of the model fit, Beverton-Holt parameters \code{lambda},
#'      \code{K}, and \code{alpha}, and the transformed data set.
#'
#' @references Otto and Day (2007). A Biologist's Guide to Mathematical Modeling
#'     in Ecology and Evolution. Princeton University Press. Page 185.
#' @export
fit_Bev_Holt <- function(clean_data, show_plot = TRUE, xmax = 30) {

  # remove NA rows from clean_data
  tmp <- na.omit(clean_data)

  # Modify the fraction of females per bean so that all resource densities are
  # constant.
  data <- data.frame(x = 10 / tmp$beans, y = tmp$f * 10 / tmp$beans)

  # fit a Beverton-Holt Model.
  # The alpha (a) in the equation below is the inverse (1/alpha) from Otto and
  # Day (2007), because the model won't converge otherwise.

  # get starter values for `nls` using the log of y; this function is a bit
  # finicky and occasionally won't run otherwise.
  start <- nls(log(y) ~ log((lambda * x) / (1 + x / alpha)), data,
               list(lambda = 1, alpha = 1))

  # use starter values to get estimates for untransformed data.
  logfit <- nls(y ~ (lambda * x) / (1 + x / alpha), data, coef(start))

  # make predictions on a range of values from 0:100 individuals
  data_predict <- data.frame(x = 0:100,
                             y = predict(logfit, newdata = data.frame(x=0:100)))

  # Get coefficients and calculate K
  a <- coef(logfit)[2]
  L <- as.numeric(format(coef(logfit)[1], digits = 3))
  K <- as.numeric(format((L - 1) * a, digits = 3))

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
    ggplot2::annotate(
      "text",
      x = 3,
      y = 90,
      size = 3,
      label = paste0("italic(lambda) == ", L),
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
              fit    = logfit,
              lambda = as.numeric(coef(logfit)[1]),
              K      = as.numeric((coef(logfit)[1] - 1) * coef(logfit)[2]),
              alpha  = as.numeric(coef(logfit)[2]),
              data   = data)
  )
}
