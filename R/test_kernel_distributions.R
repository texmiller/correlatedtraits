# fit to a Poisson distribution
fit_p <- function(params, y) {
  -sum(log(dpois(y, lambda = params)))
}

# fit to Conway-Maxwell-Poisson distribution
fit_cmp <- function(params, y) {
  -sum(log(compoisson::dcom(y, lambda = params[1], nu = params[2])))
}

# fit to a Negative Binomial distribution
fit_nb <- function(params, y) {
  -sum(log(dnbinom(y, mu = params[1], size = params[2])))
}

# fit to a Poisson Inverse-Gaussian distribution
fit_pig <- function(params, y) {
  -sum(log(
    gamlss.dist::dSICHEL(y, mu = params[1], sigma = params[2], nu = -0.5)
  ))
}

# fit to a Sichel distribution
fit_sichel <- function(params, y) {
  -sum(log(
    gamlss.dist::dSICHEL(y, mu = params[1], sigma = params[2], nu = params[3])
  ))
}

# get frequency of count data
get_frequency <- function(x) {
  x <- hist(x, breaks = seq(-0.5, 40.5, by = 1), plot = FALSE)
  x <- x$counts / sum(x$counts)
  x[which(x == 0)] <- NA
  return(x)
}

# calculate AIC
myAIC <- function(x) {
  2 * x$value + 2 * length(x$par)
}

# make AIC table
myAICtab <- function(models, aics) {
  out      <- data.frame(Model = models, AIC = aics)
  out      <- out[order(out$AIC), ]
  out$dAIC <- out$AIC - out$AIC[1]
  return(out)
}

# Make three data sets: all the data, males only, and females only
a <- abs(na.omit(data$dist))
f <- abs(na.omit(data[data$sex == 'f', ]$dist))
m <- abs(na.omit(data[data$sex == 'm', ]$dist))

# Number of patches
patch <- 0:40

## fit to data
# Poisson distribution
a_p <- optim(par = 5, fn = fit_p, y = a, method = "Br", lower = 0, upper = 10)
f_p <- optim(par = 5, fn = fit_p, y = f, method = "Br", lower = 0, upper = 10)
m_p <- optim(par = 5, fn = fit_p, y = m, method = "Br", lower = 0, upper = 10)

# Conway-Maxwell-Poisson distribution
a_c <- optim(par = c(5, 1), fn = fit_cmp, y = a, method = "L-B", lower = 1E-6)
f_c <- optim(par = c(5, 1), fn = fit_cmp, y = f, method = "L-B", lower = 1E-5)
m_c <- optim(par = c(5, 1), fn = fit_cmp, y = m, method = "L-B", lower = 1E-5)

# Negative Binomial distribution
a_n <- optim(par = c(5, 1), fn = fit_nb, y = a)
f_n <- optim(par = c(5, 1), fn = fit_nb, y = f)
m_n <- optim(par = c(5, 1), fn = fit_nb, y = m)

# Poisson Inverse-Gaussian distribution
a_i <- optim(par = c(5, 0.1), fn = fit_pig, y = a)
f_i <- optim(par = c(5, 0.1), fn = fit_pig, y = f)
m_i <- optim(par = c(5, 0.1), fn = fit_pig, y = m)

# Sichel
a_s <- optim(par = c(5, 1, 1), fn = fit_sichel, y = a)
f_s <- optim(par = c(5, 1, 1), fn = fit_sichel, y = f)
m_s <- optim(par = c(5, 1, 1), fn = fit_sichel, y = m)

# plot fits
# ALL DATA
plot(patch, get_frequency(a), pch = 20, ylim = c(0, 0.2), main = 'All Beetles')
lines(patch,
      dpois(x = patch, lambda = a_p$par[1]),
      col = "blue", lwd = 2)
lines(patch,
      compoisson::dcom(x = patch, lambda = a_c$par[1], nu = a_c$par[2]),
      col = "green", lwd = 2)
lines(patch,
      gamlss.dist::dSICHEL(x = patch, mu = a_i$par[1], sigma = a_i$par[2],
                           nu = -0.5),
      col = "purple", lwd = 2)
lines(patch,
      gamlss.dist::dSICHEL(x = patch, mu = a_s$par[1], sigma = a_s$par[2],
                           nu = a_s$par[3]),
      col="orange", lwd = 2)
lines(patch,
      dnbinom(x = patch, mu = a_n$par[1], size = a_n$par[2]),
      col = "red", lwd = 2)
points(patch, get_frequency(a), pch = 20)
legend("topright", c("Poisson", "CMP", "NegBinom", "PIG", "Sichel"), lwd=2,
       col = c("blue", "green", "red", "purple", "orange"))

# ALL FEMALES
plot(patch, get_frequency(f), pch = 20, ylim = c(0, 0.2), main = 'All Females')
lines(patch,
      dpois(x = patch, lambda = f_p$par[1]),
      col = "blue", lwd = 2)
lines(patch,
      compoisson::dcom(x = patch, lambda = f_c$par[1], nu = f_c$par[2]),
      col = "green", lwd = 2)
lines(patch,
      gamlss.dist::dSICHEL(x = patch, mu = f_i$par[1], sigma = f_i$par[2],
                           nu = -0.5),
      col = "purple", lwd = 2)
lines(patch,
      gamlss.dist::dSICHEL(x = patch, mu = f_s$par[1], sigma = f_s$par[2],
                           nu = f_s$par[3]),
      col="orange", lwd = 2)
lines(patch,
      dnbinom(x = patch, mu = f_n$par[1], size = f_n$par[2]),
      col = "red", lwd = 2)
points(patch, get_frequency(f), pch = 20)
legend("topright", c("Poisson", "CMP", "NegBinom", "PIG", "Sichel"), lwd=2,
       col = c("blue", "green", "red", "purple", "orange"))

# ALL MALES
plot(patch, get_frequency(m), pch = 20, ylim = c(0, 0.2), main = 'All Males')
lines(patch,
      dpois(x = patch, lambda = m_p$par[1]),
      col = "blue", lwd = 2)
lines(patch,
      compoisson::dcom(x = patch, lambda = m_c$par[1], nu = m_c$par[2]),
      col = "green", lwd = 2)
lines(patch,
      gamlss.dist::dSICHEL(x = patch, mu = m_i$par[1], sigma = m_i$par[2],
                           nu = -0.5),
      col = "purple", lwd = 2)
lines(patch,
      gamlss.dist::dSICHEL(x = patch, mu = m_s$par[1], sigma = m_s$par[2],
                           nu = m_s$par[3]),
      col="orange", lwd = 2)
lines(patch,
      dnbinom(x = patch, mu = m_n$par[1], size = m_n$par[2]),
      col = "red", lwd = 2)
points(patch, get_frequency(m), pch = 20)
legend("topright", c("Poisson", "CMP", "NegBinom", "PIG", "Sichel"), lwd=2,
       col = c("blue", "green", "red", "purple", "orange"))

# AICs
# ALL
a_AIC <- myAICtab(models = c("Poisson",
                             "Conway-Maxwell-Poisson",
                             "Negative Binomial",
                             "Poisson Inverse-Gaussian",
                             "Sichel"),
                  aics = c(myAIC(a_p),
                           myAIC(a_c),
                           myAIC(a_n),
                           myAIC(a_i),
                           myAIC(a_s)))

# FEMALES
f_AIC <- myAICtab(models = c("Poisson",
                             "Conway-Maxwell-Poisson",
                             "Negative Binomial",
                             "Poisson Inverse-Gaussian",
                             "Sichel"),
                  aics = c(myAIC(f_p),
                           myAIC(f_c),
                           myAIC(f_n),
                           myAIC(f_i),
                           myAIC(f_s)))

# MALES
m_AIC <- myAICtab(models = c("Poisson",
                             "Conway-Maxwell-Poisson",
                             "Negative Binomial",
                             "Poisson Inverse-Gaussian",
                             "Sichel"),
                  aics = c(myAIC(m_p),
                           myAIC(m_c),
                           myAIC(m_n),
                           myAIC(m_i),
                           myAIC(m_s)))
a_AIC
f_AIC
m_AIC
