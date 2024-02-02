graphics.off() 
rm(list=ls(all=TRUE)) 


library(rjags)
library(doSNOW)  #To run in
library(foreach) #parallel
library(truncnorm) #For truncated normal density

###
### Data
###


mydataframe = read.csv2( file="Legio_Swiss_120.csv" )
yresp = as.numeric(mydataframe[, "qPCR_count"])
xresp = as.numeric(mydataframe[, "qPCR_volume"])
N = length(yresp)
nobse = length(yresp) #Number of observations in each cluster
nclus <- 1 #Number of clusters
ntot <- nobse #Total sample size
clust <- sort(rep(1:nobse,1)) #cluster ID

clustl <- split(seq(ntot), clust) #cluster ID to compute fast
### Number of replications
nL <- 5000
nM <- 5000

###
### Information for jags
###

nburn <- 1000   #Burn-in period
niter <- 1000  #Number of iterations after burn-in
nthin <- 1     #Thinning factor

data_jags <- list(N=nobse, Y=yresp, X=xresp)
param_jags <- c("muOfLogY","sigmaOfLogY")


###
### Gamma model
###
# Package the data for shipping to JAGS:
dataList = list(
  y = yresp ,
  x = xresp ,
  N = length(yresp)
)

modelstring ="
model {
for( i in 1 : N ) {
#Likelihood
y[i] ~ dpois(mu[i]*x[i])
mu[i] ~ dlnorm(muOfLogY,1/sigmaOfLogY^2) 
}
# Prior:
muOfLogY ~ dunif( -30 , 30)
sigmaOfLogY ~ dexp(1)
}

" # close quote for modelstring
writeLines(modelstring,con="model.txt")
# RUN THE CHAINS
require(rjags)
parameters = c("muOfLogY","sigmaOfLogY","mu" ) 
adaptSteps = 1000         # Number of steps to "tune" the samplers.
burnInSteps = 1000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=1000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.

mcmc2 = jags.model( "model.txt" , data=dataList ,
                    n.chains=nChains , n.adapt=adaptSteps )
update(mcmc2, n.iter=nburn)
mcmc2_s <- coda.samples(mcmc2, c(param_jags,"mu"), n.iter=niter, thin=nthin)



##########
########## DIC 
##########

mu_mcmc <- unlist(mcmc2_s[,"muOfLogY"])  # alpha and beta sampled
sig_mcmc <- unlist(mcmc2_s[,"sigmaOfLogY"])    # for each iteration in MCMC
lam2_mcmc <- as.matrix(mcmc2_s[,paste0("mu[",1:ntot,"]")])

nK <- length(mu_mcmc) # Final number of iterations in MCMC

mu_p <- mean(mu_mcmc)  # Posterior estimates
sig_p <- mean(sig_mcmc)    # for the parameters
lam2_p <- colMeans(lam2_mcmc)


#####
##### mDIC via replication
#####

time_r <- proc.time()

cl <-makeCluster(4) #Number of CPU cores
registerDoSNOW(cl)

### Plug-in deviance
lp_yi_p <- foreach(m=1:nM, .combine='rbind') %dopar% {
  
  z_rep_p <- rgamma(ntot, mu_p, sig_p) #Replicated random effects
  
  lp_yij_p <- dpois(yresp, z_rep_p[clust]*xresp, log=TRUE) #log-like. for each data point
  lp_yij_p <- as.matrix(lp_yij_p)
  
  lp_yi_p <- sapply(clustl, function(a) colSums(lp_yij_p[a, ,drop=FALSE])) #log-like. for each obs. unit
  lp_yi_p
}

py_im2 <- colMeans(exp(lp_yi_p))
logL_pm2 <- log(py_im2)
Dev_pm2 <- -2*sum(logL_pm2)

var_plug2 <- 4*(colMeans(exp(2*lp_yi_p-2*matrix(rep(logL_pm2,nM),nrow=nM,byrow=TRUE)))-1)/nM
var_plug2 <- sum(var_plug2)
(sd_plug2 <- sqrt(var_plug2)) #Standard error

### Deviance over replicates
logLm2 <- foreach(k=1:nK, .combine='rbind') %dopar% {
  
  lp_yij <- matrix(0, nrow=ntot, ncol=nL) #To compute variance
  for(l in 1:nL)
  {
    z_rep <- rgamma(ntot, mu_mcmc[k], sig_mcmc[k]) #Replicated latent vars
    
    lp_yij[,l] <- dpois(yresp, z_rep[clust]*xresp, log=TRUE) #log-like. for each data point
  }
  
  lp_yi <- sapply(clustl, function(a) colSums(lp_yij[a,,drop=FALSE])) #log-like. for each obs. unit
  p_yi <- colMeans(exp(lp_yi))
  
  var_ave <- 1/nL*(colMeans(exp(2*lp_yi-2*matrix(rep(log(p_yi),nL),nrow=nL,byrow=TRUE)))-1)
  
  c(log(p_yi), sum(var_ave))
}

stopCluster(cl) #Closing multiple cores

time_r <- proc.time()-time_r; time_r/60 #8 minutes

var_ave2 <- sum(4*logLm2[,ncol(logLm2)])/nK^2
(sqrt(var_ave2)) #Standard error

logLm2 <- colMeans(logLm2[,1:(ncol(logLm2))])
Devm2 <- -2*sum(logLm2)

### DIC
(pDm2 <- Devm2 - Dev_pm2)
(DICm2 <- Devm2 + pDm2) 
sqrt(4*var_ave2+var_plug2) #Standard error

