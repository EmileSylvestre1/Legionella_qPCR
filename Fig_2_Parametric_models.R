graphics.off() 
rm(list=ls(all=TRUE)) 
source("DBDA2E-utilities.R")
library(rjags)
library(runjags)
#------------------------------------------------------------------------------
#THE DATA.
mydataframe = read.csv2( file="file.csv" )
Q = as.numeric(mydataframe[, "qPCR_count"])
C = as.numeric(mydataframe[, "culture_count"])
V_Q = as.numeric(mydataframe[, "qPCR_volume"])
V_C = as.numeric(mydataframe[, "culture_volume"])
cov_ln_R_emp = cov(log((Q+1)/V_Q),log((C+1)/V_C))
N = length(Q)

#------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  Q = Q ,
  C = C ,
  V_Q = V_Q ,
  V_C = V_C ,
  cov_ln_R_emp = cov_ln_R_emp,
  N = length(Q)
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood PLN qPCR
Q[i] ~ dpois(lambda[i]*V_Q[i])
lambda[i] ~ dlnorm(muOfLogQ,1/sigmaOfLogQ^2) 

#Likelihood PLN culture
C[i] ~ dpois(delta[i]*V_C[i])
delta[i] ~ dlnorm(muOfLogC,1/sigmaOfLogC^2) 
}
# Prior PLN qPCR
muOfLogQ ~ dunif( -100 , 100)
sigmaOfLogQ ~ dexp(1)

# Prior PLN culture
muOfLogC ~ dunif( -100, 100)
sigmaOfLogC ~ dexp(1)

# Ratio qPCR:culture
muOfLogR = muOfLogQ - muOfLogC
sigmaOfLogR = sqrt(abs((((sigmaOfLogQ) ^ 2) + ((sigmaOfLogC) ^ 2))-(2*cov_ln_R_emp)))
corr = cov_ln_R_emp/(sigmaOfLogQ*sigmaOfLogC)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("muOfLogQ","sigmaOfLogQ","muOfLogC","sigmaOfLogC","muOfLogR", "sigmaOfLogR", "corr") 
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
mcmcCoda = coda.samples( jagsModel , variable.names=parameters ,
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


# Plot CCDF ------------------------------------------------------------------------------

tiff(file = paste("Legio_individual.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1, 1.0E+0, 1.0E+1, 1.0E+2, 1.0E+3, 1.0E+4, 1.0E+5, 1.0E+6, 1.0E+7)
par(mar=c(4.1, 4.1, 1.1, 1.1))

plot(point_x, point_y, col="white", pch = 19, xlab = expression(paste(bold("qPCR:culture ratio, MPN/L, GU/L"))), ylab = "Cumulative probability", xaxt="n", xlim = range(myTicks),ylim=c(0, 1),
     font.lab = 2, log='x')
axis(side = 1, at = myTicks)

#Simulating uncertainty interval part 1:
HDI_output = matrix(1:322, nrow = 161, ncol = 2)
pb = txtProgressBar(min = 0, max = 161, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

#Plot log-normal distribution ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:80, nrow = 40, ncol = 2)
pb = txtProgressBar(min = 0, max = 40, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 40, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLogQ"], sdlog = mcmcChain[i,"sigmaOfLogQ"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 10, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 1, 0, 0.2),border = NA) 


#Plot log-normal distribution ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:80, nrow = 40, ncol = 2)
pb = txtProgressBar(min = 0, max = 40, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 40, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLogC"], sdlog = mcmcChain[i,"sigmaOfLogC"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 10, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:640, nrow = 320, ncol = 2)
pb = txtProgressBar(min = 0, max = 320, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 320, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (j / 64), meanlog = mcmcChain[i,"muOfLogR"], sdlog = mcmcChain[i,"sigmaOfLogR"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.125/8, 5, 0.125/8))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(1, 0, 0, 0.2),border = NA) 


#Simulating uncertainty interval below 1:
HDI_output = matrix(1:640, nrow = 320, ncol = 2)
pb = txtProgressBar(min = 0, max = 320, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 320, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (-j / 64), meanlog = mcmcChain[i,"muOfLogR"], sdlog = mcmcChain[i,"sigmaOfLogR"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ -(seq(0.125/8, 5, 0.125/8))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,2], rev(MatrixB[,3])), col = rgb(1, 0, 0, 0.2),border = NA) 



# Data point ------------------------------------------------------------------------------

lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLogR"]), sdlog=median(mcmcChain[,"sigmaOfLogR"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="red", lwd=2,lty=1)

lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLogQ"]), sdlog=median(mcmcChain[,"sigmaOfLogQ"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="forestgreen", lwd=2,lty=1)

lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLogC"]), sdlog=median(mcmcChain[,"sigmaOfLogC"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=2,lty=1)

point_x = sort(ratio1)
point_y = ecdf(ratio1)(sort(ratio1))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

point_x = sort(Q/V_Q)
point_y = ecdf(Q/V_Q)(sort(Q/V_Q))
points(point_x, point_y ,col = "black", bg = "forestgreen", pch = 23, lwd = 0.5, cex = 0.75)

point_x = sort(C/V_C)
point_y = ecdf(C/V_C)(sort(C/V_C))
points(point_x, point_y ,col = "black", bg = "blue", pch = 24, lwd = 0.5, cex = 0.75)


abline(v = 1, col="black", lwd=1, lty=2)

dev.off()
