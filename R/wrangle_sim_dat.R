## Purpose: This code wrangles the raw output from additional simulations that Brad ran and transported to me via hard disk
## The goal is to pull pertinent quantities from the many thousands of output files and package them into a handy dataframe.
## Author: Tom Miller
## Date: 5/17/2018

##set wd to my local drive where I downloaded output files
setwd("D:/ct_simulations/output/complete")
# Random question for brad: Why are there different numnbers of replicates for each treatment (1043, 1008, 1040, etc.).
# I will only use the first thousand, but this is odd. -- ANSWER: I think it's because some simulations did not work
# which is why there are missing numbers in the output list (hence tryCatch in the loop)
library(stringr)
library(tidyverse)

##----------------- Define parameters of simulation output

# number of simulation replicates
reps <- 1000
# heritability treatments
h2_trts <- c("default","swap","none","eqlo","eqhi")
# variance treatments
var_trts <- c("default","double")
#initialize tibble to store output
sim_dat <- tibble()

##----------------- Loop over data files and grab data from the right-most occupied patch in the last generation (20)
for(i in 1:length(h2_trts)){
  for(j in 1:length(var_trts)){
    for(k in 1:reps){

      #sim_out <- read.csv(paste0("ct_",h2_trts[i],"_",var_trts[j],"_",str_pad(k,4,pad = "0"),".csv"),header=F)
      sim_out <- tryCatch(
        read.csv(paste0("ct_",h2_trts[i],"_",var_trts[j],"_",str_pad(k,4,pad = "0"),".csv"),header=F),
        error=function(e) e
      )
      #some of the numbers that we are looping over do not exist as data files, not sure why.
      #We are just skipping over this when it happens. It means final sample sizes will actually be <1000
      if(inherits(sim_out, "error")) next

      names(sim_out)[c(8,9,10,14,1)] <- c("rho_G","rho_E","gen","rightmost_dat","dat_names")

      hold <- sim_out %>%
        filter(gen==1 | gen==20) %>%
        select(rho_G,rho_E,gen,dat_names,rightmost_dat) %>%
        spread(key = dat_names, value = rightmost_dat) %>%
        mutate(h2 = h2_trts[i],
               P_var = var_trts[j])

      sim_dat = bind_rows(sim_dat,hold)

    }
  }
}

##----------------- Write output to file
write_csv(sim_dat,"ct_sim_out_tidy.csv")
