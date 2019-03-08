## Data simulation for test of maternal effects estimation in the animal model
## For revision of Ochocki et al Am Nat manuscript
## 3-8-2019
library(tidyverse)

## parameters of experimental design (matching our actual experiment)
n_sires <- 50
n_dams_per_sire <- 3
n_fullsib_females <- 8

# consider only one trait (dispersal)

## biological parameters
## three scenarios are meant to vary the relative amount of sire and dam variance, holding
## total phenotypic variance (V_G_d + V_M_d + V_E_d) constant
## additive genetic variance
V_G <- c(0.23,0.115,0.0)
## maternal variance
V_M <- c(0.0,0.115,0.23)
## environmental variance
V_E <- 0.18
## grand trait mean
mu <- 1.63

sim_dat <-
  tibble(sireID = rep(1:n_sires,each=n_dams_per_sire*n_fullsib_females),
       damID = rep(1:(n_sires*n_dams_per_sire),each = n_fullsib_females),
       crossID = interaction(sireID,damID),
       beetleID = interaction(sireID,damID,rep(1:8,times = n_sires*n_dams_per_sire)),
       a_sire1 = rep(rnorm(n=n_sires,mean=0,sd=sqrt(V_G[1])),each=n_dams_per_sire*n_fullsib_females),
       a_sire2 = rep(rnorm(n=n_sires,mean=0,sd=sqrt(V_G[2])),each=n_dams_per_sire*n_fullsib_females),
       a_sire3 = rep(rnorm(n=n_sires,mean=0,sd=sqrt(V_G[3])),each=n_dams_per_sire*n_fullsib_females),
       a_dam1 = rep(rnorm(n=n_dams_per_sire*n_sires,mean=0,sd=sqrt(V_G[1])),each=n_fullsib_females),
       a_dam2 = rep(rnorm(n=n_dams_per_sire*n_sires,mean=0,sd=sqrt(V_G[2])),each=n_fullsib_females),
       a_dam3 = rep(rnorm(n=n_dams_per_sire*n_sires,mean=0,sd=sqrt(V_G[3])),each=n_fullsib_females),
       midparent1 = (a_sire1 + a_dam1)/2, ## a little unsure about this wrt Eq. 7 - should I sample this instead?
       midparent2 = (a_sire2 + a_dam2)/2,
       midparent3 = (a_sire3 + a_dam3)/2,
       m_dam1 = rep(rnorm(n=n_dams_per_sire*n_sires,mean=0,sd=sqrt(V_M[1])),each=n_fullsib_females),
       m_dam2 = rep(rnorm(n=n_dams_per_sire*n_sires,mean=0,sd=sqrt(V_M[2])),each=n_fullsib_females),
       m_dam3 = rep(rnorm(n=n_dams_per_sire*n_sires,mean=0,sd=sqrt(V_M[3])),each=n_fullsib_females),
       e_i = rnorm(n=n_dams_per_sire*n_sires*n_fullsib_females,sd=sqrt(V_E)),
       lambda1 = mu + midparent1 + m_dam1 + e_i,
       lambda2 = mu + midparent2 + m_dam2 + e_i,
       lambda3 = mu + midparent3 + m_dam3 + e_i
  )
## need to add observation-level random variables as a separate, vectorized step
sim_dat$y1 <- rpois(n=nrow(sim_dat),lambda=exp(sim_dat$lambda1))
sim_dat$y2 <- rpois(n=nrow(sim_dat),lambda=exp(sim_dat$lambda2))
sim_dat$y3 <- rpois(n=nrow(sim_dat),lambda=exp(sim_dat$lambda3))

## write to file
write_csv(sim_dat,"animal_model_sim_dat.csv")
