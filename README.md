# DaubechiesSim
Generate the $\rm H\ddot{o}lder$ functions using the "DaubechiesSim" package.

## Setup
Use the following command in R to install the "DaubechiesSim" package:
```
library(devtools)
install_github("siqixu/DaubechiesSim",ref="main") 
```
## Examples
```
library(wavethresh)
library(DaubechiesSim)
library(doMC); registerDoMC(50)

T = 200     # simulation replicates
n = 10000   # sample size
s = 5       # number of covariates 
beta = 1.5
Ry = 1
filter.number = 5
resolution = 15
gammas = c(0,3,6,9,10,16)

wv.mother <- draw(filter.number = filter.number, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)

simDat <- function(s, n, t, wv.mother, filter.number, Ry, gammas, beta){
  set.seed(t)  
  x=matrix(runif(n*s,-1,1),ncol=s) 
  h_x=matrix(wv_trans_fast(0.5*x, wv.mother, filter.number, Ry, gammas, beta),ncol=s)
  simX=cbind(x,h_x)
  return(simX)
}
out <- foreach(i=1:T,.packages=c("DaubechiesSim")) %dopar% {
   simDat(s, n, t=i, wv.mother, filter.number, Ry, gammas, beta)
}
```
