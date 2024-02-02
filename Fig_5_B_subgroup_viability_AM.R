graphics.off() 
rm(list=ls(all=TRUE))
if(!is.null(dev.list())) dev.off()
source("DBDA2E-utilities.R")

#### Library #####

library("metafor")
library("jmuOutlier")
library("ggplot2")
library("DescTools")


### Import data ##########
data = read.csv2( file="Legio_meta_analysis_viability.csv")
ni = as.numeric(data[, "ni"]) 
yi = as.numeric(data[, "r_bar"]) 
vi = as.numeric(data[, "v"]) 
data$ci.lb_AM = as.numeric(data[, "ci.lb"]) 
data$ci.ub_AM = as.numeric(data[, "ci.ub"])

data_viability_qPCR = subset(data, method == 'viability_qPCR') 
data_qPCR = subset(data, method == 'qPCR') 


######### Sample statistics ########


#--Viability qPCR----------------------------------------------------------------------------
y = as.numeric(data_viability_qPCR$r_bar)
v = as.numeric(data_viability_qPCR$v)
N = length(y) 


# Package the data for shipping to JAGS:
dataList = list(
  y = y,
  v = v,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

w[i] = 1/v[i]
y[i] ~ dnorm(theta[i], w[i])

# random effect
theta[i] ~ dnorm(mu, 1/tau^2)
}

# prior distributions
mu ~ dunif(-10, 10)
tau ~ dexp(1)


}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("mu","tau") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=20000       # Total number of steps in chains to save.
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

mean_viability = median(mcmcChain[,"mu"])
CI_viability = HDIofMCMC(mcmcChain[,"mu"], 0.95)

#------qPCR------------------------------------------------------------------------

y2 = as.numeric(data_qPCR$r_bar)
v2 = as.numeric(data_qPCR$v)
N2 = length(y) 


# Package the data for shipping to JAGS:
dataList = list(
  y = y2,
  v = v2,
  N = N2
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

w[i] = 1/v[i]
y[i] ~ dnorm(theta[i], w[i])

# random effect
theta[i] ~ dnorm(mu, 1/tau^2)
}

# prior distributions
mu ~ dunif(-10, 10)
tau ~ dexp(1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it
#------------------------------------------------------------------------------
# RUN THE CHAINS
require(rjags)
parameters = c("mu","tau") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=20000       # Total number of steps in chains to save.
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

mean_qPCR = median(mcmcChain2[,"mu"])
CI_qPCR = HDIofMCMC(mcmcChain2[,"mu"], 0.95)

###### Visualizations ############

#3. Forest plot 
tiff(file = paste("Legio_subgroup_analysis_viability_AM.tiff", sep="" ), width = 7000, height = 4000, units = "px", res = 1000)
## The parameter "height" from the tiff function needs to be adapted 
#depending on the number of studies available.

height = length(data$ni)+3

forest(yi, ci.lb = as.numeric(data$ci.lb_AM), ci.ub = as.numeric(data$ci.ub_AM), slab=paste(data$author, data$year, sep=", "),  xlim = c(-7, 6), ylim = c(-6, height),
       alim = c(-1,4), steps = 6,
       ilab = cbind(data$Assay, data$ni, formatC(data$mu, format="f", digits=1), formatC(data$sigma, format="f", digits=1)),
       ilab.xpos = c(-3.5, -2.3, -1.8, -1.3), cex=.75,
       mlab="", showweights=FALSE, xlab = ("Log ratio"), psize=0.8, 
       refline = 0, digits=1L, rows=c(7:4,-1:-4))

points(yi[1], 7, pch=22, col="black", bg="deeppink", lwd=1)
points(yi[2], 6, pch=22, col="black", bg="deeppink", lwd=1)
points(yi[3], 5, pch=22, col="black", bg="deeppink", lwd=1)
points(yi[4], 4, pch=22, col="black", bg="deeppink", lwd=1)

points(yi[5], -1, pch=22, col="black", bg="deeppink", lwd=1)
points(yi[6], -2, pch=22, col="black", bg="deeppink", lwd=1)
points(yi[7], -3, pch=22, col="black", bg="deeppink", lwd=1)
points(yi[8], -4, pch=22, col="black", bg="deeppink", lwd=1)

op <- par(cex=.75, font=2)
text(-7, height , "Author(s) and Year", pos = 4)
text(-3.6, height , "Assay")
text(-2.2, height , "n")
text(-2.0, height , expression(bold(mu)))
text(-1.4, height , expression(bold(sigma)))
text(4.5, height , "Log10 AM ratio [95% CI]")

### add text for the subgroups
text(-7, c(8,0), pos=4, c("qPCR","Viability qPCR"))

op <- par(cex=.75, font=1)

### add text for the subgroups
text(-7, c(2,-6), pos=4, c("RE Model for Subgroup","RE Model for Subgroup"))

addpoly(x = mean_qPCR,  ci.lb = CI_qPCR[1], ci.ub = CI_qPCR[2], row= 2, efac=0.75,
        mlab="", cex=1, col="deeppink")
addpoly(x = mean_viability,  ci.lb = CI_viability[1], ci.ub = CI_viability[2], row= -6, efac=0.75,
        mlab="", cex=1, col="deeppink")

dev.off()
