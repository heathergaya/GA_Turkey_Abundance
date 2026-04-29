library(nimble)
library(coda)
library(nimbleSCR)
library(parallel)
source('NimbleDists.R') #warning message expected

#extract k value from job array:
args <- commandArgs(trailingOnly = TRUE)
k <- as.numeric(sub("k=", "", args[grep("^k=", args)]))

stage1 <- nimbleCode({
  # SCR model
  # Uses NimbleSCR built in functions for speed
  for(k in 1:2){
    sigma[k] ~ dexp(1)
    logsigma[k] <- log(sigma[k])
    lam0[k] ~ dbeta(1,1)
    loglam0[k] <- log(lam0[k])
  }

  ## Habitat intensity for each habitat window is the same - homogeneous landscape
  habIntensity[1:numHabWindows] <- 1
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)

  for(i in 1:nind) { # loop over all transmitter animals
    for(m in 1:672){ #transmitter data, day doesn't matter, not a biased walk
      u[i,m,1] ~ dnorm(s[i,1], sd= sigma[sex.u[i]+1]) # x coord, using Normal(mean, sd)
      u[i,m,2] ~ dnorm(s[i,2], sd =sigma[sex.u[i]+1]) # y coord, using Normal(mean, sd)
    } #end k
    s[i,1:2] ~ dbernppAC( #located proportional to intensity surface, which is currently homogeneous
      lowerCoords = habLoCoords[1:numHabWindows,1:2],
      upperCoords = habUpCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)

    ## Using a modified form of dbinomLocal_normal from the nimbleSCR package to get detection probability using a binary detection formulation. Code for this function is inside the "NimbleDists.R" script
    for(t in 1:nFortnight){
      y[i,1:y_dim,t] ~ dbinomLocal_normal_p(detNums = detNums[i,t],
                                            size = trials[1:n.traps,t],
                                            p0 = lam0[sex.y[i]+1],
                                            detIndices = detIndices[i,,t],
                                            s = s[i,1:2],
                                            sigma = sigma[sex.y[i]+1],
                                            trapCoords = trapCoords[1:n.traps,1:2],
                                            localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                            localTrapsNum = nTraps[1:n.cells],
                                            resizeFactor = resizeFactor,
                                            habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                            indicator = 1,
                                            lengthYCombined = lengthYCombined)
    } #end t
  } #end i
})


allsims <- readRDS('simulation_forSMR.rds')
mysim <- allsims[[k]]  
consts <- mysim$nimconst
dat <- mysim$nimdat
inits <- mysim$niminits

## Nimble SCR will get cranky if y's column dimension is only 1. 
#If so, add in a dummy column to y and det_indices:
if(consts$y_dim == 1){
  tempy <- array(-1, c(dim(dat$y)[1], 2, dim(dat$y)[3]))
  tempy[,1,] <- dat$y
  dat$y <- tempy
  #same with detindices
  tempy2 <- array(-1, c(dim(dat$detIndices)[1],2,dim(dat$detIndices)[3]))
  tempy2[,1,] <- dat$detIndices
  dat$detIndices <- tempy2
  consts$y_dim <- 2
}

params <- c("sigma", "lam0", 'logsigma', 'loglam0')


cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c("dat", "consts", "inits", "params", "stage1", 'dbinomLocal_normal_p', 'rbinomLocal_normal_p'))
system.time(
  nim.out <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    library(nimbleSCR)
    prepnim <- nimbleModel(code = stage1, constants = consts,
                           data = dat, inits = inits, calculate = T)
    prepnim$calculate()
    mcmcnim <- configureMCMC(prepnim, monitors = params, print = T)
    nimMCMC <- buildMCMC(mcmcnim) #actually build the code for those samplers; gives some warnings because of nimbleSCR
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;
    Compnim <- compileNimble(nimMCMC, project = prepnim) # compile the samplers next
    Compnim$run(niter = 30000, nburnin = 20000, thin = 1)
    return(as.mcmc(as.matrix(Compnim$mvSamples)))
  }) #this will take awhile and not produce any noticeable output.
)
nim.out <- as.mcmc.list(nim.out)
stopCluster(cl)

# Stage 2 --------------------------------------------------------
## Grab logsigma, loglam0 parameters to make priors for the smaller models
# males = 2, females = 1
#
prior_means <- array(NA, c(2,2)) #2 variables, 2 sexes
prior_vcov <- array(NA, c(2,2,2))
for(t in 1:2){
  samps_m1 <- as.matrix(nim.out[,c(paste0("logsigma[", t, ']'), paste0("loglam0[", t, ']'))])
  xbar_m1 <- colMeans(samps_m1)
  xvar_m1 <- var(samps_m1)
  prior_means[,t] <- xbar_m1
  prior_vcov[,,t] <- xvar_m1
}

### This will unfortunately go fairly slowly:
stage2 <- nimbleCode({

  #homogenous landscape
  habIntensity[1:numHabWindows] <- 1
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)

  # Use the posterior estimates from the marked population as informative priors for the unmarked population
  # Multivariate normal distribution
  for(k in 1:2){ #2 sexes, independent
    meanlog_sig_lam0[1:2,k] ~ dmnorm(prior_means[1:2,k], cov=prior_vcov[1:2,1:2,k])
    sigma[k] <- exp(meanlog_sig_lam0[1,k])
    lam0[k] <- exp(meanlog_sig_lam0[2,k])

    psi[k] ~ dunif(0, 1) #proportion real
  }

  for(i in 1:M) { # loop over all individuals
    z[i] ~ dbern(psi[sex.y[i]+1]) #real/not real status depends on sex

    #get an activity center regardless if real/fake:
    s[i,1:2] ~ dbernppAC( #proportional to intensity surface
      lowerCoords = habLoCoords[1:numHabWindows,1:2],
      upperCoords = habUpCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  } #end i

  for(k in 1:2){ #sexes
    p[1:n.traps,k] <- get_p2sex( #get p(at least one det) summed across all individuals
      s = s[1:M, 1:2],
      M = M,
      k = k,
      sex = sex.y[1:M],
      lam0 = lam0[1:2],
      sigma = sigma[1:2],
      trapCoords = trapCoords[1:n.traps,1:2],
      ntraps = n.traps,
      z = z[1:M])
    for(t in 1:nFortnight){
      for(j in 1:n.traps) {
        nsum[j,t,k] ~ dbin(p[j,k], 14) # 14 days per fortnight of effort
      } #end m
    } #end t
  } #end k

    N <- sum(z[1:M]) #total estimated abundance
    NM <- sum(z[1:M]*sex.y[1:M]) #females are sex = 0
    NF <- N-NM
}) #end model

consts$prior_vcov <- prior_vcov
consts$prior_means <- prior_means
params2 <- c('N', 'psi', 'sigma', 'lam0', 'NM', 'NF')
inits$meanlog_sig_lam0 <- prior_means
inits$psi <- inits$psi[1:2] #can't remember why I made this 4 originally

cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c("dat", "consts", "inits", "params2", "stage2", 'get_p2sex'))
system.time(
  nim.out2 <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    library(nimbleSCR)
    prepnim <- nimbleModel(code = stage2, constants = consts,
                           data = dat, inits = inits, calculate = T)
    prepnim$calculate()
    mcmcnim <- configureMCMC(prepnim, monitors = params2, print = T)
    nimMCMC <- buildMCMC(mcmcnim) #actually build the code for those samplers; gives some warnings because of nimbleSCR
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;
    Compnim <- compileNimble(nimMCMC, project = prepnim) # compile the samplers next
    Compnim$run(niter = 20000, nburnin = 1000, thin = 5)
    return(as.mcmc(as.matrix(Compnim$mvSamples)))

  }) #this will take awhile and not produce any noticeable output.
)

stage2.res <- mcmc.list(nim.out2)
saveRDS(stage2.res, paste0('./SimAbund/SMR_',k, '.rds'))

stopCluster(cl)


###### Bands Only Simulation ####
#uses the same data except the GPS information:

BandSMR <- nimbleCode({
  for(k in 1:2){
    sigma[k] ~ dexp(1) 
    logsigma[k] <- log(sigma[k])
    lam0[k] ~ dbeta(1,1) 
    loglam0[k] <- log(lam0[k])
    psi[k] ~ dunif(0, 1) #proportion real 
  }
  
  ## Habitat intensity for each habitat window is the same - homogeneous landscape
  habIntensity[1:numHabWindows] <- 1
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(ii in 1:nind) { # loop over all tagged animals
    z1[ii] ~ dbern(psi[sex.y[ii]+1]) #real/not real status depends on sex
    s[ii,1:2] ~ dbernppAC( #located proportional to intensity surface, which is currently homogeneous
      lowerCoords = habLoCoords[1:numHabWindows,1:2],
      upperCoords = habUpCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    for(t in 1:nFortnight){
      y[ii,1:y_dim,t] ~ dbinomLocal_normal_p(detNums = detNums[ii,t],
                                            size = trials[1:n.traps,t],
                                            p0 = lam0[sex.y[ii]+1],
                                            detIndices = detIndices[ii,,t],
                                            s = s[ii,1:2],
                                            sigma = sigma[sex.y[ii]+1],
                                            trapCoords = trapCoords[1:n.traps,1:2],
                                            localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                            localTrapsNum = nTraps[1:n.cells],
                                            resizeFactor = resizeFactor,
                                            habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                            indicator = z1[ii],
                                            lengthYCombined = lengthYCombined)
    } #end t
  } #end ii
  
  ### start unmarked bit:
  for(i in 1:M) { # loop over all individuals
    z[i] ~ dbern(psi[sex.y[i]+1]) #real/not real status depends on sex
    
    #get an activity center regardless if real/fake:
    s2[i,1:2] ~ dbernppAC( #proportional to intensity surface
      lowerCoords = habLoCoords[1:numHabWindows,1:2],
      upperCoords = habUpCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  } #end i
  
  for(k in 1:2){ #sexes
    p[1:n.traps,k] <- get_p2sex( #get p(at least one det) summed across all individuals
      s = s2[1:M, 1:2],
      M = M,
      k = k,
      sex = sex.y[1:M],
      lam0 = lam0[1:2],
      sigma = sigma[1:2],
      trapCoords = trapCoords[1:n.traps,1:2],
      ntraps = n.traps,
      z = z[1:M])
    for(t in 1:nFortnight){
      for(j in 1:n.traps) {
        nsum[j,t,k] ~ dbin(p[j,k], 14) # 14 days per fortnight of effort
      } #end m
    } #end t
  } #end k
  
  N <- sum(z[1:M]) #total estimated abundance
  NM <- sum(z[1:M]*sex.y[1:M]) #females are sex = 0
  NF <- N-NM
}) #end model

consts <- mysim$nimconst
dat <- mysim$nimdat
inits <- mysim$niminits
params2 <- c('N', 'psi', 'sigma', 'lam0', 'NM', 'NF')
inits$psi <- inits$psi[1:2] #can't remember why I made this 4 originally
inits$s2 <- inits$s
inits$z1 <- inits$z

cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c("dat", "consts", "inits", "params2", "BandSMR", 'get_p2sex','dbinomLocal_normal_p', 'rbinomLocal_normal_p'))
system.time(
  nim.out3 <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    library(nimbleSCR)
    prepnim <- nimbleModel(code = BandSMR, constants = consts,
                           data = dat, inits = inits, calculate = T)
    prepnim$calculate()
    mcmcnim <- configureMCMC(prepnim, monitors = params2, print = T)
    nimMCMC <- buildMCMC(mcmcnim) #actually build the code for those samplers; gives some warnings because of nimbleSCR
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;
    Compnim <- compileNimble(nimMCMC, project = prepnim) # compile the samplers next
    Compnim$run(niter = 20000, nburnin = 1000, thin = 5)
    return(as.mcmc(as.matrix(Compnim$mvSamples)))
    
  }) #this will take awhile and not produce any noticeable output.
)

stage2.res_Band <- mcmc.list(nim.out3)
saveRDS(stage2.res_Band, paste0('./SimAbund/BandsSMR_',k, '.rds'))

stopCluster(cl)