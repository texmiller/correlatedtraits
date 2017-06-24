#' Simulate a half-sib pedigree.
#'
#' This function is slightly modified from the 'simPedHS' function in the
#'     'nadiv' package.
#'
#' @param s Number of sires.
#' @param d Number of dams per sire.
#' @param n Number of offspring per dam.
#' @export
sim_halfsib = function (s, d, n)
{
  if (n < 2) {
    stop("must have more than 1 offspring per family (n > or = 2)")
  }
  sires     <- seq(1, s, 1)
  dams      <- seq(s + 1, (s * d) + s, 1)
  offspring <- seq((s * d) + s + 1, (s * d * n) + (s * d) + s, 1)

  ped <- data.frame(animal = c(sires, dams, offspring),
                    dam    = c(rep(NA,   length(sires)),
                               rep(NA,   length(dams)),
                               rep(dams, each = n)),
                    sire   = c(rep(NA, length(sires)),
                               rep(NA, length(dams)),
                               rep(sires, each = d * n)))
  return(ped)
}
