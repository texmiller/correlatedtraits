## 4/2/2019
## New figures to visualize data and biplots of fitted quantities from animal model
library(tidyverse)
library(mvtnorm)
library(RColorBrewer)
library(colorspace)
library(colorRamps)

## parameters of experimental design (matching our actual experiment)
n_sires <- 50
n_dams_per_sire <- 3
n_fullsib_females <- 8

# real parameter estimates from Brad's updated analysis
## Dispersal
mu_d <- 1.642
V_Gd <- 0.03
V_Ed <- 0.251
V_Md <- 0.102
## Fertility
mu_r <- 2.737
V_Gr <- 0.033
V_Er <- 0.298
V_Mr <- 0.018
K <- 3.564
## Covariances
Cov_G <- -0.007
Cov_E <- -0.05
Cov_M <- -0.011

## assemble simulated quantities
sireID <- rep(1:n_sires,each=n_dams_per_sire*n_fullsib_females)
damID <- rep(1:(n_sires*n_dams_per_sire),each = n_fullsib_females)
crossID <- interaction(sireID,damID)
beetleID <- interaction(sireID,damID,rep(1:n_fullsib_females,times = n_sires*n_dams_per_sire))
#sires
a_j <- rmvnorm(n = n_sires, mean = c(0,0), sigma = matrix(c(V_Gd, Cov_G, Cov_G, V_Gr),nrow=2))
ad_j <- rep(a_j[,1],each=n_dams_per_sire*n_fullsib_females)
ar_j <- rep(a_j[,2],each=n_dams_per_sire*n_fullsib_females)
#dams
a_k <- rmvnorm(n = n_sires*n_dams_per_sire, mean = c(0,0), sigma = matrix(c(V_Gd, Cov_G, Cov_G, V_Gr),nrow=2))
ad_k <- rep(a_k[,1],each=n_fullsib_females)
ar_k <- rep(a_k[,2],each=n_fullsib_females)
#additive effects
ad_jk <- (ad_j+ad_k)/2
ar_jk <- (ar_j+ar_k)/2
#maternal effects
m_k <- rmvnorm(n = n_sires*n_dams_per_sire, mean = c(0,0), sigma = matrix(c(V_Md, Cov_M, Cov_M, V_Mr),nrow=2))
md_k <- rep(m_k[,1],each=n_fullsib_females)
mr_k <- rep(m_k[,2],each=n_fullsib_females)
#environmental effects
e_i <- rmvnorm(n = n_sires*n_dams_per_sire*n_fullsib_females, mean = c(0,0), sigma = matrix(c(V_Ed, Cov_E, Cov_E, V_Er),nrow=2))
ed_i <- e_i[,1]
er_i <- e_i[,2]
#expected values
loglambda_d <- mu_d + ad_jk + md_k + ed_i
loglambda_r <- mu_r + ar_jk + mr_k + er_i
#data realizations
y_d <- rpois(n = n_sires*n_dams_per_sire*n_fullsib_females, lambda = exp(loglambda_d))
B <- sample(c(1,3,5,10),size=n_sires*n_dams_per_sire*n_fullsib_females,replace=T)
y_r <- rpois(n = n_sires*n_dams_per_sire*n_fullsib_females, lambda = exp(loglambda_r) / (1 + (exp(loglambda_r)-1)/(K*B)))


## Graphics
win.graph(width = 8, height = 8)
par(mfrow=c(2,2),mar=c(5,5,2,1))
plot(jitter(y_d),jitter(y_r/B),pch=1,cex=1.2,cex.lab=1.2,
     xlab = "Observed dispersal distance (# patches)",
     ylab = "Observed fertility (beetles/female/bean)")
title("A",adj=0)
plot(ad_jk,ar_jk,pch=16,cex=1.2,cex.lab=1.2,
     xlab = expression(paste("Additive genetic effects on dispersal ( ",a[jk]^d," )")),
     ylab = expression(paste("Additive genetic effects on fertility ( ",a[jk]^r," )")))
title("B",adj=0)
plot(md_k,mr_k,pch=16,cex=1.2,cex.lab=1.2,
     xlab = expression(paste("Maternal effects on dispersal ( ",m[k]^d," )")),
     ylab = expression(paste("Maternal effects on fertility ( ",m[k]^r," )")))
title("C",adj=0)
plot(ed_i,er_i,pch=16,cex=1.2,cex.lab=1.2,
     xlab = expression(paste("Environmental effects on dispersal ( ",e[i]^d," )")),
     ylab = expression(paste("Environmental effects on fertility ( ",e[i]^r," )")))
title("D",adj=0)
