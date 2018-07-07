library(SDMTools)
library(swatches)
## need to load sim output (takes a while)

col_ramp <- read_gpl("RdYlBu_9.gpl",use_names = F)
correlations <- c(-.9,-.75,-.5,-.25,0,.25,.5,.75,.9)
max_extent <- max(extent_dat_beetle$mean_extent)
beetle_lines <- extent_dat_beetle %>% filter(rho_G==-0.37, rho_E==-0.16)

win.graph()
par(mfrow=c(1,2),mar=c(4,4,1,1))
with(extent_dat_beetle,{
  plot(rho_G,mean_extent,type="n",
       xlab=expression(paste("Genetic correlation (",rho,G,")")),
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
  points(c(-.37,-.37),beetle_lines$mean_extent,pch=c(16,1),cex=2)
  #arrows(-1.1,beetle_lines$mean_extent,-.37,beetle_lines$mean_extent,length=0)
  #arrows(-.37,0,-.37,beetle_lines$mean_extent,length=0)
  title("A",adj=0)
  
  plot(rho_G,CV_extent,type="n",
       xlab=expression(paste("Genetic correlation (",rho,G,")")),
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
  points(c(-.37,-.37),beetle_lines$CV_extent,pch=c(16,1),cex=2)
  title("B",adj=0)
  legend.gradient(pnts = cbind(x =c(-.8,-.6,-.8,-.6), y =c(60,60,40,40)), 
                  cols = col_ramp, limits = c(-0.9, 0.9),
                  title = expression(paste("Environmental correlation (",rho,E,")")),cex=0.8)
})
