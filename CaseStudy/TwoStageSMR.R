### Model code for Turkey abundance in Georgia

### Load packages
library(nimble)
library(coda)
library(nimbleSCR)
library(parallel)


#Stage 1 --------------------------------------------------------
#Only considering two stages of turkeys - males or females
stage1_2sex <- nimbleCode({

  for(jj in 1:nYears){
    beta0.loglam0[jj] ~ T(dnorm(log(0.1), sd = 0.5), -Inf, 0) # mean lam0
    beta0.logsigma[jj] ~ dnorm(log(20), sd = 2) # mean sigma
  }

  beta1.loglam0[1] <- 0 #difference because of sex (male)
  beta1.loglam0[2] ~ dnorm(0, sd = .5) #difference because of sex (male)

  beta1.logsigma[1] <- 0 #females are base level
  beta1.logsigma[2] ~ dnorm(0, sd = .5)

  for(ss in 1:2){ #two sexes
    # Detection parameter- lam0
    beta2.loglam0[ss] ~ dnorm(0, sd = .5)

    # spatial scale parameter - sigma
    beta2.logsigma[ss] ~ dnorm(0, sd = .5)

    # Starting values
    loglam0[ss,1] <- beta0.loglam0[1] + beta1.loglam0[ss] + beta2.loglam0[ss]
    lam0[ss,1] <- exp(loglam0[ss,1])
    logsigma[ss,1] <- beta0.logsigma[1] + beta1.logsigma[ss] + beta2.logsigma[ss]
    sigma[ss,1] <- exp(logsigma[ss,1])

    for(t in 2:(nWeeks*nYears)){
      # autoregression on noise for lam0
      loglam0[ss,t] <- beta0.loglam0[year[t]] + beta1.loglam0[ss] + beta2.loglam0[ss]*week[t] #+ epsilon.lam0[ss,week[t]]
      lam0[ss,t] <- exp(loglam0[ss,t])

      # autoregression on noise for sigma
      logsigma[ss,t] <- beta0.logsigma[year[t]] + beta1.logsigma[ss] + beta2.logsigma[ss]*week[t] # + epsilon.sigma[ss,week[t]]
      sigma[ss,t] <- exp(logsigma[ss,t])
    } #end week-years
  } #end ss

  # SCR model
  # Uses NimbleSCR built in functions for speed

  ## Habitat intensity for each habitat window is the same - homogeneous landscape for now
  habIntensity[1:numHabWindows] <- 1
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)

  for(i in 1:nind) { # loop over all transmitter birds
    for(t in ind.t[i,1]:ind.t[i,2]){ #however many weeks this bird was present on the landscape w/ working transmitter
      for(k in 1:nTelemLocs[i,t]){ #transmitter data, day doesn't matter, not a biased walk
        u[i,k,1,t] ~ dnorm(s[i,1,t], sd= sigma[sex[i],t]) # x coord, using Normal(mean, sd)
        u[i,k,2,t] ~ dnorm(s[i,2,t], sd =sigma[sex[i],t]) # y coord, using Normal(mean, sd)
      } #end k

      s[i,1:2,t] ~ dbernppAC( #located proportional to intensity surface, which is currently homogeneous
        lowerCoords = habLoCoords[1:numHabWindows,1:2],
        upperCoords = habUpCoords[1:numHabWindows,1:2],
        logIntensities = logHabIntensity[1:numHabWindows],
        logSumIntensity = logSumHabIntensity,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)


      ## Using a modified form of dbinomLocal_normal from the nimbleSCR package to get detection probability using a binary detection formulation. Code for this function is inside the "NimbleDists.R" script
      y[i,1:y_dim,t] ~ dbinomLocal_normal_p(detNums = detNums[i,t],
                                            size = trials[1:n.traps,t], #max number detections per trap per day (1)
                                            p0 = lam0[sex[i],t],
                                            detIndices = detIndices[i,,t],
                                            s = s[i,1:2,t],
                                            sigma = sigma[sex[i],t],
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

source('NimbleDists.R') #warning message expected

stage1_fornim <- readRDS('Stage1_inputs.rds')
dat <- stage1_fornim$dat
consts <- stage1_fornim$consts
consts$week <- rep(1:4, 6)
consts$year <- rep(1:6, each = 4)
inits <- stage1_fornim$inits
#Should have made better inits when creating these objects, but alas
inits$beta0.loglam0 <- c(-3.80787721, -4.29729166, -3.75941738, -4.66385864, -3.63632653, -4.81011155)
inits$beta0.logsigma <- rep(-1, 6)
inits$beta2.loglam0 <- c(.4, -.3)
inits$beta2.logsigma <- c(0, 0)

inits$s[237,,2] <- c(15, 3)
inits$s[240,,2] <- c(15, 3)


params <- c("sigma", "lam0",  'loglam0', 'logsigma',
            "beta0.loglam0", "beta1.loglam0","beta2.loglam0",
            "beta0.logsigma", "beta1.logsigma","beta2.logsigma")

source('runMCMCbites2.R')
cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c("dat", "consts", "inits", "params", "stage1_2sex", 'dbinomLocal_normal_p', 'rbinomLocal_normal_p','runMCMCbites', 'getModelState', 'restartMCMCbites', 'collectMCMCbites', 'getStateVariableNames', 'setModelState','getMCMCstate', 'setMCMCstate'))
system.time(
  nim.out <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    library(nimbleSCR)
    prepnim <- nimbleModel(code = stage1_2sex, constants = consts,
                           data = dat, inits = inits, calculate = T)
    conf <- configureMCMC(prepnim, monitors = params, print = T)#, enableWAIC = TRUE)
    Rmcmc <- buildMCMC(conf) #actually build the code for those samplers
    Compnim <- compileNimble(list(model = prepnim,
                                  mcmc = Rmcmc),
                             showCompilerOutput = F)
    Cmcmc <- Compnim$mcmc
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;

    c <- round(runif(1),4)
    runMCMCbites( mcmc = Cmcmc,                           ## Compiled MCMC
                  model = Cmodel,                         ## Compiled model code
                  conf = conf,                            ## MCMC configuration
                  bite.size = 750,                        ## Number of iterations per bite
                  bite.number = 8000,                        ## Number of MCMC bites
                  path = paste0("Turkeys_Stage1_Ts", "/chain",c),          ## Directory where MCMC outputs will be saved
                  save.rds = TRUE)                        ## Option to save the state of the model
  }) #this will take awhile and not produce any noticeable output.
)


# Stage 2 --------------------------------------------------------
# Assuming stage 1 gives you an object called stage1.res:
# Grab logsigma, loglam0 parameters to make priors for the smaller models
# males = 2, females = 1

nim.out <- readRDS('Feb20_Stage1_res.rds')
#
nWeeks <- 4 #technically fortnights, not weeks
nyears <- 6 #2020-2025
stage1.res <- readRDS('stage1_res.rds')
prior_means <- array(NA, c(2,2,nWeeks*nyears)) #2 variables, 2 sexes, nweeks*nyears
prior_vcov <- array(NA, c(2,2,2,nWeeks*nyears))
for(t in 1:(nWeeks*nyears)){
  samps_m1 <- as.matrix(stage1.res[,c(paste0("logsigma[2, ", t, ']'), paste0("loglam0[2, ", t, ']'))])
  xbar_m1 <- colMeans(samps_m1)
  xvar_m1 <- var(samps_m1)
  prior_means[,2,t] <- xbar_m1
  prior_vcov[,,2,t] <- xvar_m1

  samps_f1 <- as.matrix(stage1.res[,c(paste0("logsigma[1, ", t, ']'), paste0("loglam0[1, ", t, ']'))])
  xbar_f1 <- colMeans(samps_f1)
  xvar_f1 <- var(samps_f1)
  prior_means[,1,t] <- xbar_f1
  prior_vcov[,,1,t] <- xvar_f1
}

### This will unfortunately go slowly:
stage2_bysex <- nimbleCode({
  for(k in 1:2){
    for(t in 1:nYears){
      a0[t,k] ~ dnorm(0, sd = 1)
    }
    a3[k] ~ dnorm(0, sd = 1) #cropscape
  } #end k

  for(k in 1:2){ #males and females; males = 2, females = 1
    for(t in 1:(nWeeks*nYears)){
      for(mmm in 1:numHabWindows){
        habIntensity[mmm,k,t] <- exp(a3[k]*crop[mmm,t]+ a0[year[t],k]) #intercept + effect of crop
        logHabIntensity[mmm,k,t] <- log(habIntensity[mmm,k,t])
      }

      sumHabIntensity[k,t] <- sum(habIntensity[1:numHabWindows,k,t]) #expected number of animals in this time period
      logSumHabIntensity[k,t] <- log(sumHabIntensity[k,t])

      # Use the posterior estimates from the marked population as informative priors for the unmarked population
      # Multivariate normal distribution
      meanlog_sig_lam0[1:2,k,t] ~ dmnorm(prior_means[1:2,k,t], cov=prior_vcov[1:2,1:2,k,t])
      sigma[k,t] <- exp(meanlog_sig_lam0[1,k,t])
      lam0[k,t] <- exp(meanlog_sig_lam0[2,k,t])

      psi[k,t] <- sumHabIntensity[k,t]/M
      #psi[k, t] ~ dbeta(1, 1)

      for(i in 1:M) { # loop over all individuals
        z[i,k,t] ~ dbern(psi[k,t]) #real/not real status

        #get an activity center regardless if real/fake:
        s[i,1:2,k,t] ~ dbernppAC( #proportional to intensity surface
          lowerCoords = habLoCoords[1:numHabWindows,1:2],
          upperCoords = habUpCoords[1:numHabWindows,1:2],
          logIntensities = logHabIntensity[1:numHabWindows,k,t],
          logSumIntensity = logSumHabIntensity[k,t],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max)
      } #end i

      p[1:n.traps,k,t] <- get_p( #get p(at least one det) summed across all individuals
        s = s[1:M, 1:2,k,t],
        M = M,
        lam0 = lam0[k,t],
        sigma = sigma[k,t],
        trapCoords = trapCoords[1:n.traps,1:2],
        ntraps = n.traps,
        z = z[1:M,k,t])

      for(m in 1:n.traps) {
        nsum[m,k,t] ~ dbin(p[m,k,t], trials[m,t]) #how many days did each camera see at least one turkey of sex k
      } #end m

      N[k,t] <- sum(z[1:M,k,t]) #total turkeys of this sex in the state space

      #count up the number of activity centers in each pixel to estimate abundance at the per-pixel level:
      abund[1:numHabWindows,k,t] <- calculateDensity(s = s[1:M,1:2,k,t],
                                                     habitatGrid =  habitatGrid[1:y.max,1:x.max],
                                                     indicator = z[1:M,k,t],
                                                     numWindows = numHabWindows,
                                                     nIndividuals = M)
    } #end t (weeks)
  } #end k (sex)
}) #end model

stage2stuff <- readRDS('Stage2_Winter.rds')
dat2 <- stage2stuff$dat
consts2 <- stage2stuff$consts
inits2 <- stage2stuff$inits
inits2$a0 <- matrix(c(-1.2402034,-0.3245975, -1.5890182, 
                       -0.5807873, -2.1978113,-1.6937577, 
                       -0.6925015, -0.4840151,-1.0056680, 
                       -1.0893358, -1.1109388, -1.0021497), nrow = 6, ncol = 2, byrow = F)
inits2$a3 <- c(0.67, -.04)

consts2$prior_vcov <- prior_vcov
consts2$prior_means <- prior_means
consts2$year <- rep(1:6, each = 4)
params2 <- c('N', 'psi', 'sigma', 'lam0', 'a0' , 'a3', 'abund')
inits2$meanlog_sig_lam0 <- prior_means

source('NimbleDists.R') #warning message expected
source('runMCMCbites2.R')
cl <- makeCluster(5)
clusterExport(cl = cl, varlist = c("dat2", "consts2", "inits2", "params2", "stage2_bysex", 'get_p','runMCMCbites', 'getModelState', 'restartMCMCbites', 'collectMCMCbites', 'getStateVariableNames', 'setModelState','getMCMCstate', 'setMCMCstate'))

library(nimble)
library(coda)
library(nimbleSCR)

system.time(
  nim.out <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    library(nimbleSCR)
    prepnim <- nimbleModel(code = stage2_bysex, constants = consts2,
                           data = dat2, inits = inits2, calculate = T)
    conf <- configureMCMC(prepnim, monitors = params2, print = T)#, enableWAIC = TRUE)
    Rmcmc <- buildMCMC(conf) #actually build the code for those samplers
    Compnim <- compileNimble(list(model = prepnim,
                                  mcmc = Rmcmc),
                             showCompilerOutput = F)
    Cmcmc <- Compnim$mcmc
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;

    c <- round(runif(1),4)
    runMCMCbites( mcmc = Cmcmc,                           ## Compiled MCMC
                  model = Cmodel,                         ## Compiled model code
                  conf = conf,                            ## MCMC configuration
                  bite.size = 550,                        ## Number of iterations per bite
                  bite.number = 8000,                        ## Number of MCMC bites
                  path = paste0("Turkeys_Feb23_Ts", "/chain",c),          ## Directory where MCMC outputs will be saved
                  save.rds = TRUE)                        ## Option to save the state of the model
  }) #this will take awhile and not produce any noticeable output.
)
stopCluster(cl)
