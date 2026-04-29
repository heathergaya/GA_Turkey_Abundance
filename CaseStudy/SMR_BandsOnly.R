library(nimble)
library(coda)
library(nimbleSCR)
library(parallel)

SMR_Alt <- nimbleCode({
  ### Process model - series of closed population models ####
  for(k in 1:2){
    a3[k] ~ dnorm(0, 1)
    B2_s[k] ~ dnorm(0, 1)
    B2_l[k] ~ dnorm(0, 1)
    for(j in 1:nYears){ #yearly variation in expected abundance 
      a0[k,j] ~ dnorm(0, 1)
    }
  }
  for(j in 1:nYears){ #yearly variation in sigma/lambda
  B01[j] ~ dnorm(0, 1)
  B02[j] ~ dnorm(0, 1)
  }
  
  B1_s ~ dnorm(0, 1)
  B1_l ~ dnorm(0, 1)
  
  
  for(j in 1:nYears){ #loop over all years
      for(k in 1:2){ #males and females; males = 2, females = 1
        for(mmm in 1:numHabWindows){
          habIntensity[mmm,k,j] <- exp(a0[k,j] + a3[k]*crop[mmm,j])
          logHabIntensity[mmm,k,j] <- log(habIntensity[mmm,k,j])
        } #end mmm
        sumHabIntensity[k,j] <- sum(habIntensity[1:numHabWindows,k,j]) #expected number of animals
        logSumHabIntensity[k,j] <- log(sumHabIntensity[k,j])
        
        for(t in 1:nWeeks){
        lam0[k,t,j] <- exp(B01[j] +  B1_l*(k-1) + B2_l[k]*t) #effect of week and sex
        sigma[k,t,j] <- exp(B02[j] +  B1_s*(k-1) + B2_s[k]*t)
        
        psi[k,t,j] <- sumHabIntensity[k,j]/M #marked guys
        psi2[k,t,j] <- sumHabIntensity[k,j]/M #unmarked guys
        }} #end t and k
    for(t in 1:nWeeks){ #t again
      for(i in 1:nind) { # loop over all marked individuals
        z[i,t,j] ~ dbern(psi[sex[i],t,j]) #alive/real in year t, week j
        
      } #end i (individual) 
        } #end t
    
    for(i in 1:nind) { # loop over all marked individuals
      #get an activity center, one per year
      s[i,1:2,j] ~ dbernppAC(
        lowerCoords = habLoCoords[1:numHabWindows,1:2],
        upperCoords = habUpCoords[1:numHabWindows,1:2],
        logIntensities = logHabIntensity[1:numHabWindows,sex[i],j],
        logSumIntensity = logSumHabIntensity[sex[i],j],
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)
    } #end i (marked)
    
    for(ii in 1:M) { #all individuals, marked or not
      for(k in 1:2){
        for(t in 1:nWeeks){
          z2[ii,k,t,j] ~ dbern(psi2[k,t,j])
        } #end weeks
        
        s2[ii,k,1:2,j] ~ dbernppAC( #all activity centers, stage 2
          lowerCoords = habLoCoords[1:numHabWindows,1:2],
          upperCoords = habUpCoords[1:numHabWindows,1:2],
          logIntensities = logHabIntensity[1:numHabWindows,k,j],
          logSumIntensity = logSumHabIntensity[k,j],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max) 
      } #end k (sex)
    } #end ii (unmarked)
    
    ##### Detection Process #####
    ## starting with unmarked
    for(t in 1:nWeeks){
      for(k in 1:2){ #males and females; males = 2, females = 1
        p.all[1:n.traps,k,t,j] <- get_p( #get p(at least one det) summed across all individuals
          s = s2[1:M, k,1:2,j],
          M = M,
          lam0 = lam0[k,t,j],
          sigma = sigma[k,t,j],
          trapCoords = trapCoords[1:n.traps,1:2],
          ntraps = n.traps,
          z = z2[1:M,k,t,j])
        
        #detection of any individual, treating all as unmarked
        for(q in 1:n.traps) {
          nsum[q,k,t,j] ~ dbin(p.all[q,k,t,j], trials[q,t,j])
        } #end q
        
        N[k,t,j] <- sum(z2[1:M,k,t,j]) #total abundance estimated from unmarked 
        
        abund[1:numHabWindows,k,t,j] <- calculateDensity(s = s2[1:M,k,1:2, j],
                                                         habitatGrid =  habitatGrid[1:y.max,1:x.max],
                                                         indicator = z2[1:M,k,t,j],
                                                         numWindows = numHabWindows,
                                                         nIndividuals = M)
      } #end k (sex)
        #still inside t, j loop:
      
        #now the marked guys:
        for(i in 1:nind){ 
          #detection of tagged individuals
          y[i,1:y_dim,t,j] ~ dbinomLocal_normal_p(detNums = detNums[i,t,j], 
                                                  size = trials[1:n.traps,t,j], #max number detections
                                                  p0 = lam0[sex[i],t,j],
                                                  detIndices = detIndices[i,,t,j],
                                                  s = s[i,1:2,j],
                                                  sigma = sigma[sex[i],t,j],
                                                  trapCoords = trapCoords[1:n.traps,1:2],
                                                  localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                                  localTrapsNum = nTraps[1:n.cells],
                                                  resizeFactor = resizeFactor,
                                                  habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                  indicator = z[i,t,j]*tagged[i,t,j], #can only see if alive and already tagged
                                                  lengthYCombined = lengthYCombined)
        } #end i
      Ntot[t,j] <- sum(N[1:2, t, j])
    }# end t (weeks)
  } #end j (years)
  
})

params <- c('Ntot', 'N', 'psi', 'psi2', 'a0', 'a3', 
            'B01', 'B02', 'B1_s', 'B2_s', 'B1_l', 'B2_l',
            'abund','lam0', 'sigma')

myinputs <- readRDS('SMR2_inputs.rds')
dat <- myinputs$dat
inits <- myinputs$inits
inits$B01 <- rep(-3, 6)
inits$B02 <- rep(.5, 6)
inits$s[36,,] <- inits$s[1,,]
inits$B2_s <- c(-.05, 0)
inits$B2_l <- c(.05, 0)

consts <- myinputs$consts

source('NimbleDists.R') #warning message expected
source('runMCMCbites2.R')
cl <- makeCluster(5)
clusterExport(cl = cl, varlist = c("dat", "consts", "inits", "params", "SMR_Alt", 'get_p','runMCMCbites', 'getModelState', 'restartMCMCbites', 'collectMCMCbites', 'getStateVariableNames', 'setModelState','getMCMCstate', 'setMCMCstate', 'dbinomLocal_normal_p'))

system.time(
  nim.out <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    library(nimbleSCR)
    source('NimbleDists.R') #warning message expected
    source('runMCMCbites2.R')
    prepnim <- nimbleModel(code = SMR_Alt, constants = consts,
                           data = dat, inits = inits, calculate = T)
    prepnim$calculate()
    
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
                  bite.size = 500,                        ## Number of iterations per bite
                  bite.number = 8000,                        ## Number of MCMC bites
                  path = paste0("Turkeys_Bands_Feb20", "/chain",c),          ## Directory where MCMC outputs will be saved
                  save.rds = TRUE)                        ## Option to save the state of the model
  }) #this will take awhile and not produce any noticeable output.
)
stopCluster(cl)
