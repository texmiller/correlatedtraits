#' a function to transform clean_data into a form that can be passed to Stan.
#' @export
stan_data <- function(data) {

  # remove NAs (males) from data
  data <- na.omit(data)

  # Update data so that dam IDs are unique
  N_sire <- max(data$sire, na.rm=T)                       # Maximum sire ID
  N_dam  <- max(data$dam,  na.rm=T)                       # Maximum dams ID
  data$dam <- N_sire + (data$sire - 1) * N_dam + data$dam # Make unique dam IDs

  # make sire and dam IDs continuous, starting at 1
  data$sire <- as.numeric(as.factor(data$sire))
  data$dam  <- as.numeric(as.factor(data$dam))

  # Add offspring IDs to data
  data$oid <- (max(data$dam) + 1):(max(data$dam) + nrow(data))

  # set some metadata variables
  N   <- nrow(data)                 # total number of observations
  NS  <- length(unique(data$sire))  # number of sires
  ND  <- length(unique(data$dam))   # number of dams

  ###... Setup Model Matrices --------------------------------------------------
  # vector to connect sires to dams
  DS <- unique(data[, c('dam', 'sire')])$sire

  # model matrix to connect offspring to dams
  ID <- model.matrix(oid ~ as.factor(dam) - 1, data = data)

  # matrix of zeros to use as mean of sire effects
  zeros <- matrix(0, nrow = NS, ncol = 3)

  # matrix of ones to change betas from a scaler to a vector of scalers
  ones <- matrix(1, nrow = ND, ncol = 1)

  first <- matrix(c(1, 0, 0), nrow = 3, ncol = 1)
  secnd <- matrix(c(0, 1, 0), nrow = 3, ncol = 1)
  third <- matrix(c(0, 0, 1), nrow = 3, ncol = 1)

  # data that will be fit to the model
  dis <- abs(data$dist)
  off <- fit_Bev_Holt(data, show_plot = FALSE)$data$y
  den <- matrix(fit_Bev_Holt(data, show_plot = FALSE)$data$x, ncol = 1)

  return(list(N  = N,
              NS = NS,
              ND = ND,
              DS = DS,
              ID = ID,
              dis = dis,
              off = off,
              den = den,
              ones  = ones,
              zeros = zeros,
              first = first,
              secnd = secnd,
              third = third))
}
