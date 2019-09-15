## 4/2/2019
## New figures to visualize data and biplots of fitted quantities from animal model
library(tidyverse)
library(mvtnorm)
library(RColorBrewer)
library(colorspace)
library(colorRamps)
library(rstan)
library(bayesplot)

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


# real data and model outputs ---------------------------------------------
# load the workspace that Brad shared with Tom (17 April 2019 ) - it's big
load("D:/Dropbox/Manuscripts/Brad demography dispersal correlation/bev_holt_maternal_effects.rda")
load("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/bev_holt_maternal_effects.rda")

## check on raw data - looks right
ggplot(ped)+
  geom_histogram(aes(x=abs(dist)))
ggplot(ped)+
  geom_point(aes(x=beans,y=f))
ggplot(ped)+
  geom_histogram(aes(x=f/beans))
## Is there a correlation in the raw data?
cor.test(abs(ped$dist), ped$f/ped$beans)

breeding_vals <- unique(data.frame(rstan::summary(fits, pars = "a")$summary)) 
breeding_vals$trait <- rep(c("d","r"),times=nrow(breeding_vals)/2)

maternal_vals <- unique(data.frame(rstan::summary(fits, pars = "m")$summary))
maternal_vals$trait <- rep(c("d","r"),times=nrow(maternal_vals)/2)

enviro_vals <- unique(data.frame(rstan::summary(fits, pars = "d")$summary))
enviro_vals$trait <- rep(c("d","r"),times=nrow(enviro_vals)/2)

win.graph(width = 8, height = 8)
par(mfrow=c(2,2),mar=c(5,5,2,1))

plot(jitter(abs(ped$dist)),jitter(ped$f/ped$beans),pch=1,cex=1.2,cex.lab=1.2,lwd=2,
     xlab = "Observed dispersal distance (# patches)",
     ylab = "Observed fertility (beetles/female/bean)")
title("A",adj=0)

plot(breeding_vals[seq(1,nrow(breeding_vals),by=2),"mean"],
     breeding_vals[seq(2,nrow(breeding_vals),by=2),"mean"],pch=16,cex=1.2,cex.lab=1.2,
     xlab = expression(paste("Additive genetic effects on dispersal ( ",a[jk]^d," )")),
     ylab = expression(paste("Additive genetic effects on fertility ( ",a[jk]^r," )")))
title("B",adj=0)

plot(maternal_vals[seq(1,nrow(maternal_vals),by=2),"mean"],
     maternal_vals[seq(2,nrow(maternal_vals),by=2),"mean"],pch=16,cex=1.2,cex.lab=1.2,
     xlab = expression(paste("Maternal effects on dispersal ( ",m[k]^d," )")),
     ylab = expression(paste("Maternal effects on fertility ( ",m[k]^r," )")))
title("C",adj=0)

plot(enviro_vals[seq(1,nrow(enviro_vals),by=2),"mean"],
     enviro_vals[seq(2,nrow(enviro_vals),by=2),"mean"],pch=16,cex=1.2,cex.lab=1.2,
     xlab = expression(paste("Environmental effects on dispersal ( ",e[i]^d," )")),
     ylab = expression(paste("Environmental effects on fertility ( ",e[i]^r," )")))
title("D",adj=0)


# Sleuthing the ridges ----------------------------------------------------

# did some estimates not converge?? no, these look good
par(mfrow=c(1,3))
hist(breeding_vals$Rhat)
hist(maternal_vals$Rhat)
hist(enviro_vals$Rhat)

## are there some really wide credible intervals?
hist(maternal_vals$X97.5.-maternal_vals$X2.5.)
abline(v=1,col="red")

## do the big cred intervals fall on the ridge?
mat_d <- maternal_vals %>% 
  filter(trait=="d") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > 1, 1, 0))
mat_r <- maternal_vals %>% 
  filter(trait=="r") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > 1, 1, 0))  
mat_all <- full_join(mat_d,mat_r,by="id")

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[mat_all$bigwidth.x==1],mat_all$mean.y[mat_all$bigwidth.x==1],col="red")
points(mat_all$mean.x[mat_all$bigwidth.y==1],mat_all$mean.y[mat_all$bigwidth.y==1],col="blue")

## well that's interesting. How many estimates does this affect?
mat_d %>% filter(bigwidth==1) %>% summarise(n())
## how many dams have offspring with fertility data but no dispersal data?
ped %>% filter(is.na(dist) & !is.na(f)) %>% select(dam) %>% summarise(n())
ped %>% filter(!is.na(dist) & !is.na(f)) %>% select(dam) %>% distinct() %>% summarise(n())
ped %>% filter(is.na(dist) & !is.na(f)) %>% select(dam) %>% distinct() %>% summarise(n())
ped %>% filter(!is.na(dist) & is.na(f)) %>% select(dam) %>% distinct() %>% summarise(n())
ped %>% select(dam) %>% distinct() %>% summarise(n())
ped %>% group_by(dam) %>% summarise(mean_d = mean(dist),
                                    mean_r = mean(f)) %>% group_by(is.na(mean_d)) %>% summarise(n())
#how many dams have NO dispersal data
ped %>% filter(!is.na(dist)) %>% group_by(dam) %>% summarise(n = n()) %>% filter(n==0)

## do again with enviro
hist(enviro_vals$X97.5.-enviro_vals$X2.5.)
abline(v=1.9,col="red")

env_d <- enviro_vals %>% 
  filter(trait=="d") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > 1.9, 1, 0))
env_r <- enviro_vals %>% 
  filter(trait=="r") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > 1.9, 1, 0))  
env_all <- full_join(env_d,env_r,by="id")

plot(env_all$mean.x,env_all$mean.y)
points(env_all$mean.x[env_all$bigwidth.x==1],env_all$mean.y[env_all$bigwidth.x==1],col="red")
points(env_all$mean.x[env_all$bigwidth.y==1],env_all$mean.y[env_all$bigwidth.y==1],col="blue")

env_d %>% filter(bigwidth==1) %>% summarise(n())
env_r %>% filter(bigwidth==1) %>% summarise(n())
env_all %>% filter(bigwidth.x==1 & bigwidth.y==1) %>% summarise(n())

## were there also problems with genetic effects?

## do again with enviro
hist(breeding_vals$X97.5.-breeding_vals$X2.5.)
abline(v=.85,col="red")

bv_d <- breeding_vals %>% 
  filter(trait=="d") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > .8, 1, 0))
bv_r <- breeding_vals %>% 
  filter(trait=="r") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > .8, 1, 0))  
bv_all <- full_join(bv_d,bv_r,by="id")

plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[bv_all$bigwidth.x==1],bv_all$mean.y[bv_all$bigwidth.x==1],col="red")
points(bv_all$mean.x[bv_all$bigwidth.y==1],bv_all$mean.y[bv_all$bigwidth.y==1],col="blue")

bv_d %>% filter(bigwidth==1) %>% summarise(n())
bv_r %>% filter(bigwidth==1) %>% summarise(n())
bv_all %>% filter(bigwidth.x==1 & bigwidth.y==1) %>% summarise(n())

## can I find the ridge just based on where we have missing data?
maternal_vals_all <- (data.frame(rstan::summary(fits, pars = "m")$summary))
maternal_vals_all$trait <- rep(c("d","r"),times=nrow(maternal_vals_all)/2)
mat_d <- maternal_vals_all %>% 
  filter(trait=="d") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > 1, 1, 0))
mat_r <- maternal_vals_all %>% 
  filter(trait=="r") %>% 
  mutate(id = row_number(),
         width = X97.5.-X2.5.,
         bigwidth = ifelse(width > 1, 1, 0))  
mat_all <- full_join(mat_d,mat_r,by="id")
dim(ped);dim(mat_all)

## indices with missing data
nodata <- which(is.na(ped$dist) & is.na(ped$f))
fonly <- which(is.na(ped$dist) & !is.na(ped$f))
## most moms had at least some dispersal data, but can I tally how many?
mom_dist <- ped %>% group_by(dam) %>% select(dist) %>% filter(!is.na(dist)) %>% summarise(n_d = n())
mom_f <- ped %>% group_by(dam) %>% select(f) %>% filter(!is.na(f)) %>% summarise(n_f = n())
mom_dist_f <- full_join(mom_dist,mom_f,by="dam") %>% mutate(dist_f_diff = n_d - n_f) %>% mutate(dist_f_diff = ifelse(is.na(n_d),-n_f,dist_f_diff))
hist(mom_dist_f$dist_f_diff)
diffs <- full_join(ped,mom_dist_f,by="dam")$dist_f_diff

win.graph()
par(mfrow=c(1,3))
plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[fonly],bv_all$mean.y[fonly],col="red")
points(bv_all$mean.x[nodata],bv_all$mean.y[nodata],col="blue")
abline(0,1/rstan::summary(fits, pars = "corrG")$summary[,"mean"][2],lwd=2,lty=2)

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[fonly],mat_all$mean.y[fonly],col="red")
points(mat_all$mean.x[nodata],mat_all$mean.y[nodata],col="blue")
abline(0,1/rstan::summary(fits, pars = "corrME")$summary[,"mean"][2],lwd=2,lty=2)

plot(env_all$mean.x,env_all$mean.y)
points(env_all$mean.x[fonly],env_all$mean.y[fonly],col="red")
points(env_all$mean.x[nodata],env_all$mean.y[nodata],col="blue")
abline(0,1/rstan::summary(fits, pars = "corrE")$summary[,"mean"][2],lwd=2,lty=2)


win.graph()
par(mfrow=c(2,3))

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[diffs <= -1],mat_all$mean.y[diffs <= -1],col="red")

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[diffs <= -2],mat_all$mean.y[diffs <= -2],col="red")

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[diffs <= -3],mat_all$mean.y[diffs <= -3],col="red")

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[diffs <= -4],mat_all$mean.y[diffs <= -4],col="red")

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[diffs <= -5],mat_all$mean.y[diffs <= -5],col="red")

plot(mat_all$mean.x,mat_all$mean.y)
points(mat_all$mean.x[diffs <= -6],mat_all$mean.y[diffs <= -6],col="red")


par(mfrow=c(2,3))
plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[diffs <= -1],bv_all$mean.y[diffs <= -1],col="red")

plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[diffs <= -2],bv_all$mean.y[diffs <= -2],col="red")

plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[diffs <= -3],bv_all$mean.y[diffs <= -3],col="red")

plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[diffs <= -4],bv_all$mean.y[diffs <= -4],col="red")

plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[diffs <= -5],bv_all$mean.y[diffs <= -5],col="red")

plot(bv_all$mean.x,bv_all$mean.y)
points(bv_all$mean.x[diffs <= -6],bv_all$mean.y[diffs <= -6],col="red")


# New plot idea -----------------------------------------------------------
# not sure I want to show the biplots after all

corrs = rstan::extract(fits, pars=c("corrE","corrME","corrG"))

ggplot(ped)+
  geom_point(aes(x=abs(dist), y=f/beans))
                 
mcmc_areas(data.frame(corrs), pars=c("corrG.1.2","corrME.1.2","corrE.1.2"))+ 
  ggplot2::scale_y_discrete(labels = c(expression(rho[G]), expression(rho[M]),expression(rho[E])))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14,face="bold"))

## fuck ggplot

pdf("C:/Users/tm9/Desktop/git local/correlatedtraits_tom/Ochocki_correlated_traits/Final publication files/Figures/Fig1_dingbatsF.pdf",
    width = 8, height = 4, useDingbats = F)
par(mfrow=c(1,2),mar=c(5,5,2,1))

plot(jitter(abs(ped$dist)),jitter(ped$f/ped$beans),pch=1,cex=1,cex.lab=1,lwd=1,
     xlab = "Observed dispersal distance (# patches)",
     ylab = "Observed fertility (beetles/female/bean)")
title("A",adj=0)

plot(density(data.frame(corrs)$corrG.2.1),ylim=c(0,5),xlim=c(-1,1),type="n",main=" ",
     xlab="Correlation",ylab="Posterior probability density",cex.lab=1)
abline(v=0,col="grey")
lines(density(data.frame(corrs)$corrG.2.1),lwd=3,col=alpha("#1b9e77",.5))
arrows(median(data.frame(corrs)$corrG.2.1),0,
      median(data.frame(corrs)$corrG.2.1),
       density(data.frame(corrs)$corrG.2.1)[["y"]][which.min(abs(density(data.frame(corrs)$corrG.2.1)[["x"]]-median(data.frame(corrs)$corrG.2.1)))],
      lwd=5,col="#1b9e77",code=0)

lines(density(data.frame(corrs)$corrME.2.1),lwd=3,col=alpha("#d95f02",.5))
arrows(median(data.frame(corrs)$corrME.2.1),0,
       median(data.frame(corrs)$corrME.2.1),
       density(data.frame(corrs)$corrME.2.1)[["y"]][which.min(abs(density(data.frame(corrs)$corrME.2.1)[["x"]]-median(data.frame(corrs)$corrME.2.1)))],
       lwd=5,col="#d95f02",code=0)

lines(density(data.frame(corrs)$corrE.2.1),lwd=3,col=alpha("#7570b3",.5))
arrows(median(data.frame(corrs)$corrE.2.1),0,
       median(data.frame(corrs)$corrE.2.1),
       density(data.frame(corrs)$corrE.2.1)[["y"]][which.min(abs(density(data.frame(corrs)$corrE.2.1)[["x"]]-median(data.frame(corrs)$corrE.2.1)))],
       lwd=5,col="#7570b3",code=0)
title("B",adj=0)

legend("topright",legend=c(expression(rho[G]),expression(rho[M]),expression(rho[E])),
       lwd=4,
       col=c("#1b9e77","#d95f02","#7570b3"),
       bty="n")

dev.off()

## over/under zero?
sum(data.frame(corrs)$corrG.2.1 < 0) / sum(data.frame(corrs)$corrG.2.1 > 0)
sum(data.frame(corrs)$corrME.2.1 < 0) / sum(data.frame(corrs)$corrME.2.1 > 0)
sum(data.frame(corrs)$corrE.2.1 < 0) / sum(data.frame(corrs)$corrE.2.1 > 0)


# Update Table 1 ----------------------------------------------------------

params <- c("beta","G","E","ME","P","h2","K","corrG","corrE","corrME")
summ_tab <- rstan::summary(fits, pars = params)$summary

# Print the ped data frame to post on datadryad along for final publication
write.csv(ped,"C:/Users/tm9/Desktop/git local/correlatedtraits_tom/Ochocki_correlated_traits/Final publication files/QG_exp_dat.csv",
          row.names = F)
