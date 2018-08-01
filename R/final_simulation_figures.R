library(tidyverse)
library(data.table)
library(scales)
library(SDMTools)
library(swatches)
## need to load sim output then a data manipulation step
dat <- read.csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/ct_sim_out_tidy.csv")
dat <- read.csv("D:/Dropbox/Manuscripts/Brad demography dispersal correlation/ct_sim_out_tidy.csv")

extent_dat <- dat %>%
  filter(gen == 20,
         location == "rightmost1") %>%
  group_by(P_var,h2,rho_G,rho_E) %>%
  summarise(mean_extent = mean(patch),
            CV_extent = sd(patch) / mean(patch))

extent_dat_beetle <- extent_dat %>%
  filter(h2 == "default"| h2 == "none", P_var == "default")

## some miscellany I will need below
col_ramp <- read_gpl("R/RdYlBu_9.gpl",use_names = F)
correlations <- c(-.9,-.75,-.5,-.25,0,.25,.5,.75,.9)
max_extent <- max(extent_dat_beetle$mean_extent)
beetle_lines <- extent_dat_beetle %>% filter(rho_G==-0.37, rho_E==-0.16)
beetle_mean_fold_change<-log(beetle_lines$mean_extent[1]/beetle_lines$mean_extent[2])
beetle_CV_fold_change<-log(beetle_lines$CV_extent[1]/beetle_lines$CV_extent[2])

# Figure 2 ----------------------------------------------------------------
## This one shows raw extent and CV values for h2==default and h2==none

win.graph()
par(mfrow=c(2,1),mar=c(4,4,1,1))
with(extent_dat_beetle,{
  plot(rho_G,mean_extent,type="n",
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final invasion extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "default"],
          mean_extent[rho_E == correlations[i] & h2 == "default"],
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "none"],
          mean_extent[rho_E == correlations[i] & h2 == "none"],
          lwd=2,lty=2,col=col_ramp[i])
    }
  #lines(rho_G[rho_E == -.16 & h2 == "default"],
  #      mean_extent[rho_E == -.16 & h2 == "default"],
  #      lwd=2)
  #lines(rho_G[rho_E == -.16 & h2 == "none"],
  #      mean_extent[rho_E == -.16 & h2 == "none"],
  #      lwd=2,lty=2)
  points(c(-.37,-.37),beetle_lines$mean_extent,pch=21,bg=c("black","white"),lwd=2,cex=1.8)
  #arrows(-1.1,beetle_lines$mean_extent,-.37,beetle_lines$mean_extent,length=0)
  #arrows(-.37,0,-.37,beetle_lines$mean_extent,length=0)
  title("A",adj=0)

  plot(rho_G,CV_extent,type="n",
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of invasion extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "default"],
          CV_extent[rho_E == correlations[i] & h2 == "default"],
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "none"],
          CV_extent[rho_E == correlations[i] & h2 == "none"],
          lwd=2,lty=2,col=col_ramp[i])
  }
  #lines(rho_G[rho_E == -.16 & h2 == "default"],
  #      CV_extent[rho_E == -.16 & h2 == "default"],
  #      lwd=2)
  #lines(rho_G[rho_E == -.16 & h2 == "none"],
  #      CV_extent[rho_E == -.16 & h2 == "none"],
  #      lwd=2,lty=2)
  points(c(-.37,-.37),beetle_lines$CV_extent,pch=21,bg=c("black","white"),lwd=2,cex=1.8)
  title("B",adj=0)
  legend.gradient(pnts = cbind(x =c(0.2,0.5,0.2,0.5), y =c(0.2,0.2,0.1,0.1)),
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental\ncorrelation (",rho[E],")")),cex=0.8)
})


# Figure 3 ----------------------------------------------------------------

extent_dat_plot <- extent_dat %>%
  select(-CV_extent) %>%
  spread(key = h2, value = mean_extent) %>%
  mutate(fold_change_default = default / none,
         fold_change_swap = swap / none,
         fold_change_eqhi = eqhi / none,
         fold_change_eqlo = eqlo / none) %>%
  select(P_var,rho_G,rho_E,fold_change_default,fold_change_swap,fold_change_eqhi,fold_change_eqlo) %>%
  gather(fold_change_default,fold_change_swap,fold_change_eqhi,fold_change_eqlo,
         key = h2, value = fold_change)

CV_dat_plot <- extent_dat %>%
  select(-mean_extent) %>%
  spread(key = h2, value = CV_extent) %>%
  mutate(fold_change_default = default / none,
         fold_change_swap = swap / none,
         fold_change_eqhi = eqhi / none,
         fold_change_eqlo = eqlo / none) %>%
  select(P_var,rho_G,rho_E,fold_change_default,fold_change_swap,fold_change_eqhi,fold_change_eqlo) %>%
  gather(fold_change_default,fold_change_swap,fold_change_eqhi,fold_change_eqlo,
         key = h2, value = fold_change)

win.graph()
par(mfrow=c(2,4),mar=c(4,4,2,1))
with(extent_dat_plot %>% filter(P_var == "default"),{
  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default"]),
          lwd=2,col=col_ramp[i])
  }
  points(-.37,beetle_mean_fold_change,pch=16,cex=1.2)
  title("A",adj=0)
  title(expression(paste(h[d]^2 > h[r]^2)))
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_swap"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_swap"]),
          lwd=2,col=col_ramp[i])
  }
  title("B",adj=0)
  title(expression(paste(h[d]^2 < h[r]^2)))
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqhi"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqhi"]),
          lwd=2,col=col_ramp[i])
  }
  title("C",adj=0)
  title(expression(paste(h[d]^2 ,"=", h[r]^2, " (high)")))
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqlo"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqlo"]),
          lwd=2,col=col_ramp[i])
  }
  title("D",adj=0)
  abline(h=0,col="darkgray")
  title(expression(paste(h[d]^2 ,"=", h[r]^2, " (low)")))
  legend.gradient(pnts = cbind(x =c(-.8,-.6,-.8,-.6), y =c(.4,.4,.2,.2)),
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental correlation (",rho[E],")")),cex=0.8)
  })
with(CV_dat_plot %>% filter(P_var == "default"),{
  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default"]),
          lwd=2,col=col_ramp[i])
  }
  points(-.37,beetle_CV_fold_change,pch=16,cex=1.2)
  title("E",adj=0)
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_swap"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_swap"]),
          lwd=2,col=col_ramp[i])
  }
  title("F",adj=0)
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqhi"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqhi"]),
          lwd=2,col=col_ramp[i])
  }
  title("G",adj=0)
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqlo"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqlo"]),
          lwd=2,col=col_ramp[i])
  }
  title("H",adj=0)
  abline(h=0,col="darkgray")

})


# 4-panel option ----------------------------------------------------------

win.graph()
par(mfrow=c(2,2),mar=c(4,4,1,1))
with(extent_dat_beetle,{
  plot(rho_G,mean_extent,type="n",cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final invasion extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "default"],
          mean_extent[rho_E == correlations[i] & h2 == "default"],
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2 == "none"],
          mean_extent[rho_E == correlations[i] & h2 == "none"],
          lwd=2,lty=2,col=alpha(col_ramp[i],.8))
  }
  lines(rho_G[rho_E == -.16 & h2 == "default"],
        mean_extent[rho_E == -.16 & h2 == "default"],
        lwd=2)
  lines(rho_G[rho_E == -.16 & h2 == "none"],
        mean_extent[rho_E == -.16 & h2 == "none"],
        lwd=2,lty=2)
  points(c(-.37,-.37),beetle_lines$mean_extent,pch=21,bg=c("black","white"),lwd=2,cex=1.8)
  #arrows(-1.1,beetle_lines$mean_extent,-.37,beetle_lines$mean_extent,length=0)
  #arrows(-.37,0,-.37,beetle_lines$mean_extent,length=0)
  title("A",adj=0)
})

  with(extent_dat_plot %>% filter(P_var == "default"),{
    plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),cex.lab=1.2,
         xlab=expression(paste("Genetic correlation (",rho[G],")")),
         ylab="log Fold-change in final extent")
    for(i in 1:9){
      lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default"],
            log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default"]),
            lwd=2,col=alpha(col_ramp[i],.8))
    }
    lines(rho_G[rho_E == -.16 & h2 == "fold_change_default"],
          log(fold_change[rho_E == -.16 & h2 == "fold_change_default"]),
          lwd=2)
    points(-.37,beetle_mean_fold_change,pch=21,bg=c("black","white"),lwd=2,cex=1.8)
    title("B",adj=0)
    abline(h=0,col="darkgray")
  })

  with(extent_dat_beetle,{
    plot(rho_G,CV_extent,type="n",cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of invasion extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "default"],
          CV_extent[rho_E == correlations[i] & h2 == "default"],
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "none"],
          CV_extent[rho_E == correlations[i] & h2 == "none"],
          lwd=2,lty=2,col=col_ramp[i])
  }
  lines(rho_G[rho_E == -.16 & h2 == "default"],
        CV_extent[rho_E == -.16 & h2 == "default"],
        lwd=2)
  lines(rho_G[rho_E == -.16 & h2 == "none"],
        CV_extent[rho_E == -.16 & h2 == "none"],
        lwd=2,lty=2)
  points(c(-.37,-.37),beetle_lines$CV_extent,pch=21,bg=c("black","white"),lwd=2,cex=1.8)
  title("C",adj=0)
  legend.gradient(pnts = cbind(x =c(0.2,0.5,0.2,0.5), y =c(0.2,0.2,0.1,0.1)),
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental\ncorrelation (",rho[E],")")),cex=0.8)
})

with(CV_dat_plot %>% filter(P_var == "default"),{
  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
  }
  lines(rho_G[rho_E == -.16 & h2 == "fold_change_default"],
        log(fold_change[rho_E == -.16 & h2 == "fold_change_default"]),
        lwd=2)
  points(-.37,beetle_CV_fold_change,pch=21,bg=c("black","white"),lwd=2,cex=1.8)
  title("D",adj=0)
  abline(h=0,col="darkgray")
})


# Barplot for beetle system -----------------------------------------------

beetle_bars <- dat %>%
  filter(gen == 20,
         location == "rightmost1",
         rho_G == -0.37,
         rho_E == -0.16,
         h2 == "default" | h2 == "none",
         P_var == "default") %>%
  group_by(h2) %>%
  summarise(log_kids = mean(Mean_r),
            log_patches = mean(Mean_D),
            CV_extent = sd(patch) / mean(patch),
            mean_extent = mean(patch))

beetles_change <- log(beetle_bars[1,2:5]/beetle_bars[2,2:5]) %>% gather(key=measurement,value=value)
beetle_names<-c("Invasion extent", "CV of extent", "Dispersal", "Fertility")

par(mar=c(4,8,1,1))
barplot(beetles_change$value,names.arg = rev(beetle_names),xlim=c(-.5,2),horiz=T,col=rep(c("black","white"),each=2),
        xlab="log Fold-change",las=2);abline(v=0)


# Trait evolution ---------------------------------------------------------
trait_dat_figure <- dat %>%
  mutate(fold_change_D = Mean_D / 1.63,
         fold_change_r = Mean_r / 2.74)%>%
  group_by(location,P_var,h2,rho_G,rho_E,gen) %>%
  select(location,P_var,h2,rho_G,rho_E,gen,Mean_D,Mean_r,fold_change_D,fold_change_r)  %>%
  summarise(log_patches = mean(Mean_D),
            log_kids = mean(Mean_r),
            CV_patches = sd(Mean_D)/log_patches,
            CV_kids = sd(Mean_r)/log_kids,
            mean_change_D = mean(fold_change_D),
            mean_change_r = mean(fold_change_r))%>%
  filter(gen==20,location=="rightmost1")

beetle_traits <- trait_dat_figure %>% filter(h2=="default",P_var=="default",
                            rho_G==-0.37, rho_E==-0.16)

win.graph()
par(mar=c(4,4,1,1))
with(trait_dat_figure %>% filter(h2=="default",P_var=="default"),{
  plot(rho_G,log(mean_change_D),type="n",ylim=c(-0.1,0.4),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in trait value")
  abline(h=0,col="darkgray")
  legend.gradient(pnts = cbind(x =c(-0.5,-0.8,-0.5,-0.8), y =c(0.15,0.15,0.05,0.05)),
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental\ncorrelation (",rho[E],")")),cex=0.8)
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i]],
          log(mean_change_D[rho_E == correlations[i]]),
          lwd=2,col=alpha(col_ramp[i],.8))
  }
  lines(rho_G[rho_E == -.16],
        log(mean_change_D[rho_E == -.16]),
        lwd=2)
  points(beetle_traits$rho_G,log(beetle_traits$mean_change_D),
         pch=21,bg=c("black","white"),lwd=2,cex=1.4)

  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i]],
          log(mean_change_r[rho_E == correlations[i]]),
          lwd=2,col=alpha(col_ramp[i],.8))
  }
  lines(rho_G[rho_E == -.16],
        log(mean_change_r[rho_E == -.16]),
        lwd=2)
  points(beetle_traits$rho_G,log(beetle_traits$mean_change_r),
         pch=21,bg=c("black","white"),lwd=2,cex=1.4)
  text(-.9,.35,"Dispersal",adj=0,font=3)
  text(-.9,-.025,"Fertility",adj=0,font=3)
})

# Appendix Figure ---------------------------------------------------------
win.graph()
par(mfrow=c(2,4),mar=c(4,4,2,1))
with(extent_dat_plot,{
  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  points(-.37,beetle_mean_fold_change,pch=16,cex=1.2)
  title("A",adj=0)
  title(expression(paste(h[d]^2 > h[r]^2)))
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  title("B",adj=0)
  title(expression(paste(h[d]^2 < h[r]^2)))
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  title("C",adj=0)
  title(expression(paste(h[d]^2 ,"=", h[r]^2, " (high)")))
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-0.1,0.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  title("D",adj=0)
  abline(h=0,col="darkgray")
  title(expression(paste(h[d]^2 ,"=", h[r]^2, " (low)")))
  legend.gradient(pnts = cbind(x =c(-.8,-.6,-.8,-.6), y =c(.4,.4,.2,.2)),
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental correlation (",rho[E],")")),cex=0.8)
})
with(CV_dat_plot,{
  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_default" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  points(-.37,beetle_CV_fold_change,pch=16,cex=1.2)
  title("E",adj=0)
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_swap" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  title("F",adj=0)
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqhi" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
    }
  title("G",adj=0)
  abline(h=0,col="darkgray")

  plot(rho_G,log(fold_change),type="n",ylim=c(-.1,2.5),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV of extent")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "default"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "default"]),
          lwd=2,col=col_ramp[i])
    lines(rho_G[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "double"],
          log(fold_change[rho_E == correlations[i] & h2 == "fold_change_eqlo" & P_var == "double"]),
          lty=2,lwd=2,col=col_ramp[i])
  }
  title("H",adj=0)
  abline(h=0,col="darkgray")

})


# Appendix figure for trait change ----------------------------------------

win.graph()
par(mfrow=c(2,4),mar=c(4,4,2,0.5))

with(trait_dat_figure,{
  plot(rho_G,log(mean_change_D),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in dispersal")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="default" & P_var=="default"],
          log(mean_change_D[rho_E == correlations[i] & h2=="default" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="default" & P_var=="double"],
          log(mean_change_D[rho_E == correlations[i] & h2=="default" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
    }
  points(beetle_traits$rho_G,log(beetle_traits$mean_change_D),
         pch=21,bg=c("black","white"),lwd=2,cex=1.4)
  title("A",adj=0)
  title(expression(paste(h[d]^2 > h[r]^2)))

  plot(rho_G,log(mean_change_D),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab=" ")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="swap" & P_var=="default"],
          log(mean_change_D[rho_E == correlations[i] & h2=="swap" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="swap" & P_var=="double"],
          log(mean_change_D[rho_E == correlations[i] & h2=="swap" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  title("B",adj=0)
  title(expression(paste(h[d]^2 < h[r]^2)))

  plot(rho_G,log(mean_change_D),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab=" ")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="eqhi" & P_var=="default"],
          log(mean_change_D[rho_E == correlations[i] & h2=="eqhi" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="eqhi" & P_var=="double"],
          log(mean_change_D[rho_E == correlations[i] & h2=="eqhi" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  title("C",adj=0)
  title(expression(paste(h[d]^2 ,"=", h[r]^2, " (high)")))

  plot(rho_G,log(mean_change_D),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab=" ")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="eqlo" & P_var=="default"],
          log(mean_change_D[rho_E == correlations[i] & h2=="eqlo" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="eqlo" & P_var=="double"],
          log(mean_change_D[rho_E == correlations[i] & h2=="eqlo" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  title("D",adj=0)
  title(expression(paste(h[d]^2 ,"=", h[r]^2, " (low)")))

  plot(rho_G,log(mean_change_r),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in fertility")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="default" & P_var=="default"],
          log(mean_change_r[rho_E == correlations[i] & h2=="default" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="default" & P_var=="double"],
          log(mean_change_r[rho_E == correlations[i] & h2=="default" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  points(beetle_traits$rho_G,log(beetle_traits$mean_change_r),
         pch=21,bg=c("black","white"),lwd=2,cex=1.4)
  legend.gradient(pnts = cbind(x =c(-.4,-.8,-.4,-.8), y =c(0.35,0.35,0.1,0.1)),
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental\ncorrelation (",rho[E],")")),cex=0.8)
  title("E",adj=0)

  plot(rho_G,log(mean_change_r),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab=" ")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="swap" & P_var=="default"],
          log(mean_change_r[rho_E == correlations[i] & h2=="swap" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="swap" & P_var=="double"],
          log(mean_change_r[rho_E == correlations[i] & h2=="swap" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  title("F",adj=0)

  plot(rho_G,log(mean_change_r),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab=" ")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="eqhi" & P_var=="default"],
          log(mean_change_r[rho_E == correlations[i] & h2=="eqhi" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="eqhi" & P_var=="double"],
          log(mean_change_r[rho_E == correlations[i] & h2=="eqhi" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  title("G",adj=0)

  plot(rho_G,log(mean_change_r),type="n",ylim=c(-0.2,0.5),cex.lab=1.2,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab=" ")
  abline(h=0,col="darkgray")
  for(i in 1:9){
    lines(rho_G[rho_E == correlations[i] & h2=="eqlo" & P_var=="default"],
          log(mean_change_r[rho_E == correlations[i] & h2=="eqlo" & P_var=="default"]),
          lwd=2,col=alpha(col_ramp[i],.8))
    lines(rho_G[rho_E == correlations[i] & h2=="eqlo" & P_var=="double"],
          log(mean_change_r[rho_E == correlations[i] & h2=="eqlo" & P_var=="double"]),
          lty=2,lwd=2,col=alpha(col_ramp[i],.8))
  }
  title("H",adj=0)
})



