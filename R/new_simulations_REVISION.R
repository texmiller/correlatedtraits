## This code analyzes output of simulations run by Brad
## The simulations update our original AmNat simulations by including maternal effects and 
## updated estimates for other QG quantities. 
library(tidyverse)
library(scales)

# Beetle simulation --------------------------------------------------
shuffle10 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/posterior_shuffle_10.csv")
shuffle20 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/posterior_shuffle_20.csv")
sort10 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/posterior_sorting_10.csv")
sort20 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/posterior_sorting_20.csv")
evooff10 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/posterior_evo-off_prop__10.csv")
evooff20 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/posterior_evo-off_prop__20.csv")

dat_sort <- sort20; dat_shuff <- shuffle20; dat_evooff <- evooff20
#dat_sort <- sort10; dat_shuff <- shuffle10


beetle_names<-c("Invasion extent", "CV of extent", "Dispersal", "Fertility")
#beetle_bars <- c(log(mean(dat_sort$x) / mean(dat_shuff$x)),
#  log((sd(dat_sort$x)/mean(dat_sort$x)) / (sd(dat_shuff$x)/mean(dat_shuff$x))),
#  log((1.64+mean(dat_sort$mean_GD)) / (1.64+mean(dat_shuff$mean_GD))),
#  log((2.74+mean(dat_sort$mean_Gr)) / (2.74+mean(dat_shuff$mean_Gr))))
  
beetle_bars <- c(log(mean(dat_sort$x) / mean(dat_evooff$x)),
                 log((sd(dat_sort$x)/mean(dat_sort$x)) / (sd(dat_evooff$x)/mean(dat_evooff$x))),
                 log((1.64+mean(dat_sort$mean_GD)) / (1.64+mean(dat_evooff$mean_GD))),
                 log((2.74+mean(dat_sort$mean_Gr)) / (2.74+mean(dat_evooff$mean_Gr))))

win.graph()
par(mar=c(4,8,1,1))
barplot(rev(beetle_bars),
        names.arg = rev(c("Invasion extent", "CV of extent", "Dispersal", "Fertility")),
        xlim=c(-.1,1),horiz=T,col=rep(c("black","white"),each=2),
        xlab="log Fold-change",las=2);abline(v=0)



# Generalized simulations -------------------------------------------------

sim_sort10 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/simulated_sorting_10.csv") %>% 
  mutate(treatment = "sorting",
         gen = 10,
         case = ifelse(V_ED == 0.251,"A",ifelse(V_ED == 0.128,"B","C")))
sim_sort20 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/simulated_sorting_20.csv") %>% 
  mutate(treatment = "sorting",
         gen = 20,
         case = ifelse(V_ED == 0.251,"A",ifelse(V_ED == 0.128,"B","C")))
sim_shuff10 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/simulated_shuffle_10.csv") %>% 
  mutate(treatment = "shuffle",
         gen = 10,
         case = ifelse(V_ED == 0.251,"A",ifelse(V_ED == 0.128,"B","C")))
sim_shuff20 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/simulated_shuffle_20.csv") %>% 
  mutate(treatment = "shuffle",
         gen = 20,
         case = ifelse(V_ED == 0.251,"A",ifelse(V_ED == 0.128,"B","C")))
sim_evooff10 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/simulated_evo-off_prop__10.csv") %>% #10.csv") %>% 
  mutate(treatment = "evooff",
         gen = 10,
         case = ifelse(V_ED == 0.272,"A",ifelse(V_ED == 0.192,"B","C")))
sim_evooff20 <- read_csv("C:/Users/tm9/Dropbox/Manuscripts/Brad demography dispersal correlation/revision_simulations/simulated_evo-off_prop__20.csv") %>%  #20.csv") %>% 
  mutate(treatment = "evooff",
         gen = 20,
         case = ifelse(V_ED == 0.272,"A",ifelse(V_ED == 0.192,"B","C")))

## re-assign the no-evolution C case to equal the B case
sim_evooff20$x[sim_evooff20$case=="C"] <- sim_evooff20$x[sim_evooff20$case=="B"]
sim_evooff20$mean_GD[sim_evooff20$case=="C"] <- sim_evooff20$mean_GD[sim_evooff20$case=="B"]
sim_evooff20$mean_Gr[sim_evooff20$case=="C"] <- sim_evooff20$mean_Gr[sim_evooff20$case=="B"]

## check that, for evo off, traits have the same variance as for evo-on
(sim_sort20$V_GD+sim_sort20$V_MD+sim_sort20$V_ED)[1:3]
(sim_evooff20$V_GD+sim_evooff20$V_MD+sim_evooff20$V_ED)[1:3]

(sim_sort20$V_Gr+sim_sort20$V_Mr+sim_sort20$V_Er)[1:3]
(sim_evooff20$V_Gr+sim_evooff20$V_Mr+sim_evooff20$V_Er)[1:3]

sim_dat <- bind_rows(sim_sort10,sim_sort20,sim_shuff10,sim_shuff20,sim_evooff10,sim_evooff20) %>% 
  ## use parameter to assign the simulation case
  #mutate(case = ifelse(V_ED == 0.251 | V_ED == 0.272,"A",
                      # ifelse(V_ED == 0.128 | V_ED == 0.192,"B","C"))) %>% 
  mutate(treatment = fct_relevel(treatment,"sorting","shuffle"))

extent_results <- sim_dat %>% 
  group_by(case,treatment,gen,rho_G,rho_M,rho_E) %>%  
  summarise(mean_extent = mean(x),
            CV_extent = sd(x)/mean(x)) 


# Mean extent -------------------------------------------------------------

extent_results %>% filter(gen==20) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=mean_extent,color=as.factor(rho_M),linetype=as.factor(rho_E)))+
  facet_grid(treatment~case)+
  theme_bw()

mean_extent_results <- extent_results %>% 
  select(gen,case,rho_G,rho_M,rho_E,treatment,mean_extent) %>% 
  spread(key=treatment,value = mean_extent) %>% 
  mutate(log_fold_change_A = log(sorting / shuffle),
         log_fold_change_B = log(sorting / evooff))

mean_extent_results %>% filter(gen==20) %>% 
  gather(log_fold_change_A,log_fold_change_B,key=method,value=fold_change) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=fold_change,color=as.factor(rho_M),linetype=as.factor(rho_E)))+
  facet_grid(method~case)+
  theme_bw()

###### **** Manuscript version **** ######
cors <- c(-.9,0,.9)
cor_cols <- alpha(c('#b2182b','#f4a582','#2166ac'),.75)
bw_cor_cols <- c('#969696','#636363','#252525')
win.graph()
par(mfrow=c(2,3),mar=c(5,5,2,1))
with(mean_extent_results %>% filter(gen==20),{
  plot(rho_G,log_fold_change_B,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent");abline(h=0,col="lightgray")
  title("A",adj=0)
  title(expression(paste(V[G]," < ",V[M]," < ",V[E])))
  
  for(i in 1:3){
    for(j in 1:3){
    lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j]],
          log_fold_change_B[case=="A"  & rho_E==cors[i] & rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
    }
  }
  points(-.173,beetle_bars[1],pch=16,cex=2)
  
  plot(rho_G,log_fold_change_B,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent");abline(h=0,col="lightgray")
  title("B",adj=0)
  title(expression(paste(V[G]," = ",V[M]," = ",V[E])))
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j]],
            log_fold_change_B[case=="B"  & rho_E==cors[i] & rho_M==cors[j]],
            col=cor_cols[j],lty=i,lwd=2)
    }
  }
  
  plot(rho_G,log_fold_change_B,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in final extent");abline(h=0,col="lightgray")
  title("C",adj=0)
  title(expression(paste(V[G]," >> (",V[M],"+",V[E],")")))
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j]],
            log_fold_change_B[case=="C"  & rho_E==cors[i] & rho_M==cors[j]],
            col=cor_cols[j],lty=i,lwd=2)
    }
  }
  
})
with(CV_extent_results %>% filter(gen==20),{
  plot(rho_G,log_fold_change_B,type="n",cex.lab=1.4,ylim=c(-.1,2.6),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV");abline(h=0,col="lightgray")
  title("D",adj=0)
  title(expression(paste(V[G]," < ",V[M]," < ",V[E])))
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j]],
            log_fold_change_B[case=="A"  & rho_E==cors[i] & rho_M==cors[j]],
            col=cor_cols[j],lty=i,lwd=2)
    }
  }
  points(-.173,beetle_bars[2],pch=16,cex=2)
  
  legend("topleft",title=expression(paste(rho[E])),bty="n",
         legend=cors,lty=1:3,lwd=2,cex=1.2)
  legend("top",title=expression(paste(rho[M])),bty="n",
         legend=cors,fill=cor_cols,cex=1.2)

  plot(rho_G,log_fold_change_B,type="n",cex.lab=1.4,ylim=c(-.1,2.6),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV");abline(h=0,col="lightgray")
  title("E",adj=0)
  title(expression(paste(V[G]," = ",V[M]," = ",V[E])))
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j]],
            log_fold_change_B[case=="B"  & rho_E==cors[i] & rho_M==cors[j]],
            col=cor_cols[j],lty=i,lwd=2)
    }
  }
  
  plot(rho_G,log_fold_change_B,type="n",cex.lab=1.4,ylim=c(-.1,2.6),
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fold-change in CV");abline(h=0,col="lightgray")
  title("F",adj=0)
  title(expression(paste(V[G]," >> (",V[M],"+",V[E],")")))
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j]],
            log_fold_change_B[case=="C"  & rho_E==cors[i] & rho_M==cors[j]],
            col=cor_cols[j],lty=i,lwd=2)
    }
  }
  
})
par(fig = c(0.025,0.225, 0.675, .975), new = T)
with(extent_results %>% filter(gen==20),{
  plot(rho_G,mean_extent,type="n",xlab="",xaxt='n',yaxt='n',ylab="")
  title(xlab=expression(paste(rho[G])),line=0, ylab="Final extent",cex.lab=1.2)
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            mean_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=1)
      #lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
      #      mean_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
      #      col=bw_cor_cols[j],lty=i,lwd=1)
    }
  }
})
############################################

# Extent CV ---------------------------------------------------------------

extent_results %>% filter(gen==20) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=CV_extent,color=as.factor(rho_M),linetype=as.factor(rho_E)))+
  facet_grid(treatment~case)+
  theme_bw()

CV_extent_results <- extent_results %>% 
  select(gen,case,rho_G,rho_M,rho_E,treatment,CV_extent) %>% 
  spread(key=treatment,value = CV_extent) %>% 
  mutate(log_fold_change_A = log(sorting / shuffle),
         log_fold_change_B = log(sorting / evooff))

CV_extent_results %>% filter(gen==20) %>% 
  gather(log_fold_change_A,log_fold_change_B,key=method,value=fold_change) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=fold_change,color=as.factor(rho_M),linetype=as.factor(rho_E)))+
  facet_grid(method~case)+
  theme_bw()

# Trait change ------------------------------------------------------------
mu_D <- 1.64
mu_r <- 2.74

trait_results <- sim_dat %>% 
  group_by(case,treatment,gen,rho_G,rho_M,rho_E) %>%  
  summarise(dispersal = mean(mu_D + mean_GD),
            fertility = mean(mu_r + mean_Gr)) %>% 
  gather(dispersal,fertility,key=trait,value=trait_value)

trait_results %>% filter(gen==20, rho_M==0) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=trait_value,color=treatment,linetype=as.factor(rho_E)))+
  facet_grid(trait~case)+
  theme_bw()

dispersal_fold_change <- trait_results %>% 
  filter(trait=="dispersal") %>% 
  select(gen,case,rho_G,rho_M,rho_E,treatment,trait,trait_value) %>% 
  spread(key=treatment,value = trait_value) %>% 
  mutate(log_fold_change_A = log(sorting / shuffle),
         log_fold_change_B = log(sorting / evooff))

dispersal_fold_change %>% filter(gen==20) %>% 
  gather(log_fold_change_A,log_fold_change_B,key=method,value=fold_change) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=fold_change,color=as.factor(rho_M),linetype=as.factor(rho_E)))+
  facet_grid(method~case)+
  theme_bw()

fertility_fold_change <- trait_results %>% 
  filter(trait=="fertility") %>% 
  select(gen,case,rho_G,rho_M,rho_E,treatment,trait,trait_value) %>% 
  spread(key=treatment,value = trait_value) %>% 
  mutate(log_fold_change_A = log(sorting / shuffle),
         log_fold_change_B = log(sorting / evooff))

fertility_fold_change %>% filter(gen==20) %>% 
  gather(log_fold_change_A,log_fold_change_B,key=method,value=fold_change) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=fold_change,color=as.factor(rho_M),linetype=as.factor(rho_E)))+
  facet_grid(method~case)+
  theme_bw()

dispersal_plot <-dispersal_fold_change %>% filter(gen==20)
fertility_plot <-fertility_fold_change %>% filter(gen==20)

###### **** Manuscript version **** ######
win.graph()
par(mfrow=c(1,3))
plot(dispersal_plot$rho_G,dispersal_plot$log_fold_change_B,type="n",ylim=c(-.2,.6),
     ylab="log Fold-change in trait value",cex.lab=1.4,
     xlab=expression(paste("Genetic correlation (",rho[G],")")));abline(h=0,col="lightgray")
title("A",adj=0)
title(expression(paste(V[G]," < ",V[M]," < ",V[E])))
for(i in 1:3){
  for(j in 1:3){
    lines(dispersal_plot$rho_G[dispersal_plot$case=="A" & dispersal_plot$rho_E==cors[i] & dispersal_plot$rho_M==cors[j]],
          dispersal_plot$log_fold_change_B[dispersal_plot$case=="A"  & dispersal_plot$rho_E==cors[i] & dispersal_plot$rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
    lines(fertility_plot$rho_G[fertility_plot$case=="A" & fertility_plot$rho_E==cors[i] & fertility_plot$rho_M==cors[j]],
          fertility_plot$log_fold_change_B[fertility_plot$case=="A"  & fertility_plot$rho_E==cors[i] & fertility_plot$rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
  }
}
points(c(-.173,-.173),beetle_bars[3:4],pch=16,cex=2)
text(-.9,.1,"Dispersal",adj=0,font=3)
text(-.9,-.075,"Fertility",adj=0,font=3)
legend("topleft",title=expression(paste(rho[E])),bty="n",
       legend=cors,lty=1:3,lwd=2,cex=1.2)
legend("top",title=expression(paste(rho[M])),bty="n",
       legend=cors,fill=cor_cols,cex=1.2)

plot(dispersal_plot$rho_G,dispersal_plot$log_fold_change_B,type="n",ylim=c(-.2,.6),
     ylab="log Fold-change in trait value",cex.lab=1.4,
     xlab=expression(paste("Genetic correlation (",rho[G],")")));abline(h=0,col="lightgray")
title("B",adj=0)
title(expression(paste(V[G]," = ",V[M]," = ",V[E])))
for(i in 1:3){
  for(j in 1:3){
    lines(dispersal_plot$rho_G[dispersal_plot$case=="B" & dispersal_plot$rho_E==cors[i] & dispersal_plot$rho_M==cors[j]],
          dispersal_plot$log_fold_change_B[dispersal_plot$case=="B"  & dispersal_plot$rho_E==cors[i] & dispersal_plot$rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
    lines(fertility_plot$rho_G[fertility_plot$case=="B" & fertility_plot$rho_E==cors[i] & fertility_plot$rho_M==cors[j]],
          fertility_plot$log_fold_change_B[fertility_plot$case=="B"  & fertility_plot$rho_E==cors[i] & fertility_plot$rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
  }
}
text(-.9,.24,"Dispersal",adj=0,font=3)
text(-.9,-.13,"Fertility",adj=0,font=3)

plot(dispersal_plot$rho_G,dispersal_plot$log_fold_change_B,type="n",ylim=c(-.2,.6),
     ylab="log Fold-change in trait value",cex.lab=1.4,
     xlab=expression(paste("Genetic correlation (",rho[G],")")));abline(h=0,col="lightgray")
title("C",adj=0)
title(expression(paste(V[G]," >> (",V[M],"+",V[E],")")))
for(i in 1:3){
  for(j in 1:3){
    lines(dispersal_plot$rho_G[dispersal_plot$case=="C" & dispersal_plot$rho_E==cors[i] & dispersal_plot$rho_M==cors[j]],
          dispersal_plot$log_fold_change_B[dispersal_plot$case=="C"  & dispersal_plot$rho_E==cors[i] & dispersal_plot$rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
    lines(fertility_plot$rho_G[fertility_plot$case=="C" & fertility_plot$rho_E==cors[i] & fertility_plot$rho_M==cors[j]],
          fertility_plot$log_fold_change_B[fertility_plot$case=="C"  & fertility_plot$rho_E==cors[i] & fertility_plot$rho_M==cors[j]],
          col=cor_cols[j],lty=i,lwd=2)
  }
}
text(-.9,.37,"Dispersal",adj=0,font=3)
text(-.9,-.2,"Fertility",adj=0,font=3)
##############################################


####### *** Appendix Figure: Mean Extent *** #############
win.graph()
par(mfrow=c(3,3),mar=c(5,5,2,1))
with(extent_results %>% filter(gen==20),{
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            mean_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("A",adj=0)
  title(expression(paste(V[G]," < ",V[M]," < ",V[E])))
  legend("topleft",title=expression(paste(rho[E])),bty="n",
         legend=cors,lty=1:3,lwd=2,cex=1.2)
  legend("top",title=expression(paste(rho[M])),bty="n",
         legend=cors,fill=cor_cols,cex=1.2)
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            mean_extent[case=="B"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("B",adj=0)
  title(expression(paste(V[G]," = ",V[M]," = ",V[E])))
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            mean_extent[case=="C"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("C",adj=0)
  title(expression(paste(V[G]," >> (",V[M],"+",V[E],")")))
  
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            mean_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("D",adj=0)
  
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            mean_extent[case=="B"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("E",adj=0)
  
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            mean_extent[case=="C"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("F",adj=0)
  
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            mean_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("G",adj=0)
  
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            mean_extent[case=="B"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("H",adj=0)
  
  plot(rho_G,mean_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="Final extent (# patches invaded)")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            mean_extent[case=="C"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("I",adj=0)
  
})
####################################################

####### *** Appendix Figure: Extent CV*** #############
win.graph()
par(mfrow=c(3,3),mar=c(5,5,2,1))
with(extent_results %>% filter(gen==20),{
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            CV_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("A",adj=0)
  title(expression(paste(V[G]," < ",V[M]," < ",V[E])))
  legend("topleft",title=expression(paste(rho[E])),bty="n",
         legend=cors,lty=1:3,lwd=2,cex=1.2)
  legend("top",title=expression(paste(rho[M])),bty="n",
         legend=cors,fill=cor_cols,cex=1.2)
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            CV_extent[case=="B"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("B",adj=0)
  title(expression(paste(V[G]," = ",V[M]," = ",V[E])))
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            CV_extent[case=="C"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="sorting"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("C",adj=0)
  title(expression(paste(V[G]," >> (",V[M],"+",V[E],")")))
  
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            CV_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("D",adj=0)
  
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            CV_extent[case=="B"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("E",adj=0)
  
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            CV_extent[case=="C"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="shuffle"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("F",adj=0)
  
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="A" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            CV_extent[case=="A"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("G",adj=0)
  
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="B" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            CV_extent[case=="B"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("H",adj=0)
  
  plot(rho_G,CV_extent,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="CV of final extent")
  for(i in 1:3){
    for(j in 1:3){
      lines(rho_G[case=="C" & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            CV_extent[case=="C"  & rho_E==cors[i] & rho_M==cors[j] & treatment=="evooff"],
            col=cor_cols[j],lty=i,lwd=3)
    }
  }
  title("I",adj=0)
  
})
####################################################

############# *** Appendix Figure: Dispersal *** ########

trait_results %>% filter(gen==20, rho_M==0) %>% 
  ggplot()+
  geom_line(size=1,aes(x=rho_G,y=trait_value,color=treatment,linetype=as.factor(rho_E)))+
  facet_grid(trait~case)+
  theme_bw()

dispersal_app <- trait_results %>% filter(gen==20, rho_M==0, rho_E==0, trait=="dispersal")
fertility_app <- trait_results %>% filter(gen==20, rho_M==0, rho_E==0, trait=="fertility")

par(mfrow=c(2,3),mar=c(5,5,2,1))
with(dispersal_app,{
  plot(rho_G,trait_value,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Dispersal distance (# patches)")
  lines(rho_G[case=="A" & treatment=="sorting"],trait_value[case=="A" & treatment=="sorting"],lwd=3,col=cor_cols[1])
  lines(rho_G[case=="A" & treatment=="shuffle"],trait_value[case=="A" & treatment=="shuffle"],lwd=3,col=cor_cols[2])
  lines(rho_G[case=="A" & treatment=="evooff"],trait_value[case=="A" & treatment=="evooff"],lwd=3,col=cor_cols[3])
  legend("topleft",title="Simulation treatment:",bty="n",
         legend=c("Sorting","Shuffle","No Evolution"),lwd=3,col=cor_cols,cex=1.4)
  title("A",adj=0)
  title(expression(paste(V[G]," < ",V[M]," < ",V[E])))
  
  plot(rho_G,trait_value,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Dispersal distance (# patches)")
  lines(rho_G[case=="B" & treatment=="sorting"],trait_value[case=="B" & treatment=="sorting"],lwd=3,col=cor_cols[1])
  lines(rho_G[case=="B" & treatment=="shuffle"],trait_value[case=="B" & treatment=="shuffle"],lwd=3,col=cor_cols[2])
  lines(rho_G[case=="B" & treatment=="evooff"],trait_value[case=="B" & treatment=="evooff"],lwd=3,col=cor_cols[3])
  title("B",adj=0)
  title(expression(paste(V[G]," = ",V[M]," = ",V[E])))
  
  plot(rho_G,trait_value,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Dispersal distance (# patches)")
  lines(rho_G[case=="C" & treatment=="sorting"],trait_value[case=="C" & treatment=="sorting"],lwd=3,col=cor_cols[1])
  lines(rho_G[case=="C" & treatment=="shuffle"],trait_value[case=="C" & treatment=="shuffle"],lwd=3,col=cor_cols[2])
  lines(rho_G[case=="C" & treatment=="evooff"],trait_value[case=="C" & treatment=="evooff"],lwd=3,col=cor_cols[3])
  title("C",adj=0)
  title(expression(paste(V[G]," >> (",V[M],"+",V[E],")")))
})
with(fertility_app,{
  plot(rho_G,trait_value,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fertility (# offspring)")
  lines(rho_G[case=="A" & treatment=="sorting"],trait_value[case=="A" & treatment=="sorting"],lwd=3,col=cor_cols[1])
  lines(rho_G[case=="A" & treatment=="shuffle"],trait_value[case=="A" & treatment=="shuffle"],lwd=3,col=cor_cols[2])
  lines(rho_G[case=="A" & treatment=="evooff"],trait_value[case=="A" & treatment=="evooff"],lwd=3,col=cor_cols[3])
  title("D",adj=0)
  
  plot(rho_G,trait_value,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fertility (# offspring)")
  lines(rho_G[case=="B" & treatment=="sorting"],trait_value[case=="B" & treatment=="sorting"],lwd=3,col=cor_cols[1])
  lines(rho_G[case=="B" & treatment=="shuffle"],trait_value[case=="B" & treatment=="shuffle"],lwd=3,col=cor_cols[2])
  lines(rho_G[case=="B" & treatment=="evooff"],trait_value[case=="B" & treatment=="evooff"],lwd=3,col=cor_cols[3])
  title("E",adj=0)
  
  plot(rho_G,trait_value,type="n",cex.lab=1.4,
       xlab=expression(paste("Genetic correlation (",rho[G],")")),
       ylab="log Fertility (# offspring)")
  lines(rho_G[case=="C" & treatment=="sorting"],trait_value[case=="C" & treatment=="sorting"],lwd=3,col=cor_cols[1])
  lines(rho_G[case=="C" & treatment=="shuffle"],trait_value[case=="C" & treatment=="shuffle"],lwd=3,col=cor_cols[2])
  lines(rho_G[case=="C" & treatment=="evooff"],trait_value[case=="C" & treatment=="evooff"],lwd=3,col=cor_cols[3])
  title("F",adj=0)
  
})
