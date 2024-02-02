graphics.off() 
rm(list=ls(all=TRUE)) 
source("DBDA2E-utilities.R")
fileNameRoot="A" 
library(rjags)
library(runjags)
#------------------------------------------------------------------------------
#THE DATA.
mydataframe = read.csv2( file="Legio_Swiss_100.csv" )
Q = as.numeric(mydataframe[, "qPCR_count"])
C = as.numeric(mydataframe[, "culture_count"])
V_Q = as.numeric(mydataframe[, "qPCR_volume"])
V_C = as.numeric(mydataframe[, "culture_volume"])
cov_ln_R_emp = cov(log(Q/V_Q),log(C/V_C))
#cov_ln_R_emp = log((cov((Q/V),(C/V))/ (mean(Q/V)*mean(C/V)))+1)
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
muOfLogQ ~ dunif( -30 , 30)
sigmaOfLogQ ~ dexp(0.1)

# Prior PLN culture
muOfLogC ~ dunif( -30, 30)
sigmaOfLogC ~ dexp(0.1)

# Ratio qPCR/culture
muOfLogR = muOfLogQ - muOfLogC
sigmaOfLogR = sqrt(abs((((sigmaOfLogQ) ^ 2) + ((sigmaOfLogC) ^ 2))-(2*cov_ln_R_emp)))


}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("muOfLogQ","sigmaOfLogQ","muOfLogC","sigmaOfLogC","muOfLogR", "sigmaOfLogR") 
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

ratio1 = (Q/V_Q)/(C/V_C)

#### Swiss 110 #########################################################################

#THE DATA.
mydataframe = read.csv2( file="Legio_Swiss_110.csv" )
Q = as.numeric(mydataframe[, "qPCR_count"])
C = as.numeric(mydataframe[, "culture_count"])
V_Q = as.numeric(mydataframe[, "qPCR_volume"])
V_C = as.numeric(mydataframe[, "culture_volume"])
#cov_ln_R_emp = cov(log(Q/V_Q),log(C/V_C))
#cov_ln_R_emp = log((cov((Q/V),(C/V))/ (mean(Q/V)*mean(C/V)))+1)
N = length(Q)

#------------------------------------------------------------------------------

# Package the data for shipping to JAGS:
dataList = list(
  Q = Q ,
  C = C ,
  V_Q = V_Q ,
  V_C = V_C ,
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
muOfLogQ ~ dunif( -30 , 30)
sigmaOfLogQ ~ dexp(0.1)

# Prior PLN culture
muOfLogC ~ dunif( -30, 30)
sigmaOfLogC ~ dexp(0.1)

# Ratio qPCR/culture
muOfLogR = muOfLogQ - muOfLogC
sigmaOfLogR = sqrt(abs((((sigmaOfLogQ) ^ 2) + ((sigmaOfLogC) ^ 2))-(2*0)))


}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("muOfLogQ","sigmaOfLogQ","muOfLogC","sigmaOfLogC","muOfLogR", "sigmaOfLogR") 
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
chainLength2 = NROW(mcmcChain2)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda2)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

ratio2 = (Q/V_Q)/(C/V_C)
########################Swiss 120 ##################################
#THE DATA.
mydataframe = read.csv2( file="Legio_Swiss_120.csv" )
Q = as.numeric(mydataframe[, "qPCR_count"])
C = as.numeric(mydataframe[, "culture_count"])
V_Q = as.numeric(mydataframe[, "qPCR_volume"])
V_C = as.numeric(mydataframe[, "culture_volume"])
cov_ln_R_emp = cov(log(Q/V_Q),log(C/V_C))
#cov_ln_R_emp = log((cov((Q/V),(C/V))/ (mean(Q/V)*mean(C/V)))+1)
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
muOfLogQ ~ dunif( -30 , 30)
sigmaOfLogQ ~ dexp(0.1)

# Prior PLN culture
muOfLogC ~ dunif( -30, 30)
sigmaOfLogC ~ dexp(0.1)

# Ratio qPCR/culture
muOfLogR = muOfLogQ - muOfLogC
sigmaOfLogR = sqrt(abs((((sigmaOfLogQ) ^ 2) + ((sigmaOfLogC) ^ 2))-(2*cov_ln_R_emp)))


}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("muOfLogQ","sigmaOfLogQ","muOfLogC","sigmaOfLogC","muOfLogR", "sigmaOfLogR") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel3 = jags.model( "model.txt" , data=dataList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel3 , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda3 = coda.samples( jagsModel3 , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain3 = as.matrix( mcmcCoda3 )
chainLength = NROW(mcmcChain3)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda3)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda3 , parName=parName)
}

ratio3 = (Q/V_Q)/(C/V_C)

# Plot CCDF ------------------------------------------------------------------------------

tiff(file = paste("Legio_Swiss_Temp.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1, 1.0E+0, 1.0E+1, 1.0E+2, 1.0E+3, 1.0E+4, 1.0E+5, 1.0E+6)
par(mar=c(4.1, 4.1, 1.1, 1.1))

plot(point_x, point_y, col="white", pch = 19, xlab = expression(paste(bold("qPCR:culture ratio"))), ylab = "Cumulative probability", xaxt="n", xlim = range(myTicks), ylim=c(0, 1),
     font.lab = 2, log='x')
axis(side = 1, at = myTicks)


#Plot Swiss 100 ------------------------------------------------------------------------------

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
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(1, 0.0784, 0.5765, 0.2),border = NA) 


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
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,2], rev(MatrixB[,3])), col = rgb(1, 0.0784, 0.5765, 0.2),border = NA) 


#Plot Swiss 110 ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:640, nrow = 320, ncol = 2)
pb = txtProgressBar(min = 0, max = 320, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 320, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (j / 64), meanlog = mcmcChain2[i,"muOfLogR"], sdlog = mcmcChain2[i,"sigmaOfLogR"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.125/8, 5, 0.125/8))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 1, 0, 0.2),border = NA) 


#Simulating uncertainty interval below 1:
HDI_output = matrix(1:640, nrow = 320, ncol = 2)
pb = txtProgressBar(min = 0, max = 320, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 320, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (-j / 64), meanlog = mcmcChain2[i,"muOfLogR"], sdlog = mcmcChain2[i,"sigmaOfLogR"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ -(seq(0.125/8, 5, 0.125/8))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,2], rev(MatrixB[,3])), col = rgb(0, 1, 0, 0.2),border = NA) 


#Plot Swiss 120 ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:640, nrow = 320, ncol = 2)
pb = txtProgressBar(min = 0, max = 320, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 320, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (j / 64), meanlog = mcmcChain3[i,"muOfLogR"], sdlog = mcmcChain3[i,"sigmaOfLogR"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.125/8, 5, 0.125/8))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0.75, 1, 0.2),border = NA) 


#Simulating uncertainty interval below 1:
HDI_output = matrix(1:640, nrow = 320, ncol = 2)
pb = txtProgressBar(min = 0, max = 320, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 320, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = plnorm(10 ^ (-j / 64), meanlog = mcmcChain3[i,"muOfLogR"], sdlog = mcmcChain3[i,"sigmaOfLogR"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ -(seq(0.125/8, 5, 0.125/8))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,2], rev(MatrixB[,3])), col = rgb(0, 0.75, 1, 0.2),border = NA) 



# Data point ------------------------------------------------------------------------------

lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLogR"]), sdlog=median(mcmcChain[,"sigmaOfLogR"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="deeppink", lwd=2,lty=1)

lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain2[,"muOfLogR"]), sdlog=median(mcmcChain2[,"sigmaOfLogR"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="forestgreen", lwd=2,lty=1)

lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain3[,"muOfLogR"]), sdlog=median(mcmcChain3[,"sigmaOfLogR"]))
lnorm_x=sort(lnorm_data)
lnorm_y = ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="deepskyblue", lwd=2,lty=1)

point_x = sort(ratio1)
point_y = ecdf(ratio1)(sort(ratio1))
points(point_x, point_y ,col = "black", bg = "deeppink", pch = 21, lwd = 0.5, cex = 0.75)

point_x = sort(ratio2)
point_y = ecdf(ratio2)(sort(ratio2))
points(point_x, point_y ,col = "black", bg = "forestgreen", pch = 23, lwd = 0.5, cex = 0.75)

point_x = sort(ratio3)
point_y = ecdf(ratio3)(sort(ratio3))
points(point_x, point_y ,col = "black", bg = "deepskyblue", pch = 24, lwd = 0.5, cex = 0.75)

abline(v = 1, col="black", lwd=1, lty=2)

dev.off()




