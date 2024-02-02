graphics.off() 
rm(list=ls(all=TRUE)) 
source("DBDA2E-utilities.R")
fileNameRoot="A" 
library(rjags)
library(runjags)
#------------------------------------------------------------------------------
#THE DATA.
mydataframe = read.csv2( file="file.csv" )
y_GM = as.numeric(mydataframe[, "r_bar_GM"])
v = as.numeric(mydataframe[, "GM_v"])
N = length(y_GM)

#------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  y_GM = y_GM,
  v = v,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

w[i] = 1/v[i]
y_GM[i] ~ dnorm(theta[i], w[i])

# random effect
theta[i] ~ dexp(lambda)
}

# prior distributions
lambda ~ dunif(0.001, 10000)



}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("lambda") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( jagsModel , variable.names=parameters ,c("deviance"),
                         n.iter=nPerChain , thin=thinSteps )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

#THE DATA.
mydataframe = read.csv2( file="legio_meta_analysis_site_specific.csv" )
y_AM = as.numeric(mydataframe[, "r_bar_AM"])
v = as.numeric(mydataframe[, "AM_v"])
N = length(y_AM)

#------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  y_AM = y_AM,
  v = v,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

w[i] = 1/v[i]
y_AM[i] ~ dnorm(theta[i], w[i])

# random effect
theta[i] ~ dlnorm(mu, 1/sigma^2)
}

# prior distributions
mu ~ dunif(-10, 10)
sigma ~ dexp(1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("mu","sigma") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel2 = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel2 , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda2 = coda.samples( jagsModel2 , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )
#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain2 = as.matrix( mcmcCoda2 )
chainLength = NROW(mcmcChain2)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda2)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda2 , parName=parName)
}

# Plot CCDF ------------------------------------------------------------------------------

tiff(file = paste("Legio_meta_site_specific_combined.tiff", sep = ""), width = 2500, height = 2500, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(-1, 0, 1, 2, 3, 4)
par(mar=c(4.1, 4.1, 1.1, 1.1))

plot(point_x, point_y, col="white", pch = 19, 
     xlab = expression(bold(Log[10]  ~ mean ~ qPCR:culture ~ ratio)),
     ylab = "Cumulative probability", xaxt="n", xlim = range(myTicks),ylim=c(0, 1),
     font.lab = 2)
axis(side = 1, at = myTicks)

abline(v = 0, col="black", lwd=1, lty=2)

#Plot exp distribution ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:400, nrow = 200, ncol = 2)
pb = txtProgressBar(min = 0, max = 200, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 200, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = pexp((j / 40), rate = mcmcChain[i,"lambda"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c((seq(0.025, 5, 0.025))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 1, 0, 0.2),border = NA) 

#Plot log-normal distribution ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:400, nrow = 200, ncol = 2)
pb = txtProgressBar(min = 0, max = 200, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 200, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm((j / 40), meanlog = mcmcChain2[i,"mu"], sdlog = mcmcChain2[i,"sigma"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c((seq(0.025, 5, 0.025))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

# Lines Exp

exp_data=rexp(n=500000, rate=median(mcmcChain[,"lambda"]))
exp_x=sort(exp_data)
exp_y = ecdf(exp_data)(sort(exp_data) )
lines(exp_x,exp_y, col="forestgreen", lwd=2,lty=1)

# Lines LN

lnorm_data=rlnorm(n=500000, meanlog = median(mcmcChain2[i,"mu"]), sdlog = median(mcmcChain2[i,"sigma"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=2,lty=1)

# Data point ------------------------------------------------------------------------------

# GM
point_x_GM = sort(y_GM)
point_y_GM = ecdf(y_GM)(sort(y_GM))
points(point_x_GM, point_y_GM ,col = "darkgreen", bg = "darkgreen", pch = 19, lwd = 1, cex = 1)

# AM
point_x_AM = sort(y_AM)
point_y_AM = ecdf(y_AM)(sort(y_AM))
points(point_x_AM, point_y_AM ,col = "darkblue", bg = "darkblue", pch = 15, lwd = 1, cex = 1)


dev.off()

