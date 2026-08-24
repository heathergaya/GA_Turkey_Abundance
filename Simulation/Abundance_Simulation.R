### Simulation for expected accuracy and precision #### 
### Uses camera array from our study ###
## Assumes one year of data, though study uses 6 ####

library(dplyr)
library(lubridate)
library(tidyr)
library(terra)
library(nimbleSCR)
library(sf)
`%notin%` <- Negate(`%in%`)


### Grab files #####
WinterCams <- read.csv('./RawData/WinterCamera_Locations.csv')
BFG <- vect('./RawData/bf grant boundary.shp')
CC <- vect('./RawData/cedarcreek_boundary2014.shp')

WinterCams <- subset(WinterCams, WinterCams$CamerID %notin% c('EB01','EB02', 'CC71', 'SA02',  'SG01',  'SG02')) #these cameras were not included in the study
traplocs <- (vect(WinterCams, geom = c('Longitude', 'Latitude')))
crs(traplocs) <- "EPSG:4326"
traps_utm <- project(traplocs, "EPSG:26917")
WinterCams <- cbind(WinterCams, geom(traps_utm)[,3:4])

ntraps <- nrow(WinterCams)
nFortnight <- 4 #4 two-week periods of sampling

### Create surface for camera trapping ####
r <- rast(crs = "EPSG:26917", extent = ext(rbind(BFG, CC)))
res(r) <- c(1000,1000)
r[] <- 1
r <- (mask(r, buffer(rbind(BFG, CC), 5000)))

coordsHabitatGridCenter <- crds(r)
colnames(coordsHabitatGridCenter) <- c("x","y")
#trapping grid 
coordsObsCenter <- WinterCams[,c("x", "y")]
colnames(coordsObsCenter) <- c("x","y")
plot(r,legend=F)
points(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"],
       pch = 16, cex = .2) 
points(coordsObsCenter[,"y"] ~ coordsObsCenter[,"x"], col="red", pch=16 ) 
habitatMask <- as.matrix(r, wide=TRUE)

## Rescale coordinates
scaledObjects <- scaleCoordsToHabitatGrid(
  coordsData = coordsObsCenter,
  coordsHabitatGridCenter = coordsHabitatGridCenter)

## Get lower and upper cell coordinates (flips raster but that's normal)
lowerAndUpperCoords <- getWindowCoords(
  scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
  scaledObsGridCenter = scaledObjects$coordsDataScaled,
  plot.check = T)

dmax = 25 #NIMBLESCR package hack that allows for sparse matrix processing 
trapLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = scaledObjects$coordsDataScaled,
                             dmax = dmax,
                             resizeFactor = 1,
                             plot.check = T)

## For plot later:
OG_array <- list(xbounds = range(scaledObjects$coordsHabitatGridCenterScaled[,1]),
                 ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2]),
                 traps = scaledObjects$coordsDataScaled)

#Simulation function #####
#simpler than real function, assumes population is the same for all 4 fortnights 
simpop <- function(sig = c(.6, .9), #female and male sigmas (scaled for nimble)
                   lam0 = c(.07, .03), #detection at activity center 
                   tagged_F = 30, #how many F are tagged 
                   tagged_M = 15 #how many M are tagged 
){ 
  
  totZ_M <- 200 # true M pop
  totZ_F <- 300 # true F pop
  
  xbounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,1])
  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2])
  
  #initial locations
  s_F <- cbind(runif(totZ_F, xbounds[1], xbounds[2]), 
               runif(totZ_F, ybounds[1], ybounds[2]))
  s_M <- cbind(runif(totZ_M, xbounds[1], xbounds[2]), 
               runif(totZ_M, ybounds[1], ybounds[2]))
  
  ## we do need our tagged friends to exist on the habitat grid though, so restrict the x values a bit:
  s_F[1:tagged_F,1] <- runif(tagged_F, 8, 20)
  s_M[1:tagged_M,1] <- runif(tagged_M, 8, 20)
  
  
  ## actual locations ###
  u_F <- array(NA, c(totZ_F, nFortnight*14*12, 2)) #14 days in a week, 12 locations per day
  u_M <- array(NA, c(totZ_M, nFortnight*14*12, 2)) 
  
  u_F[,1,] <- s_F #start at activity center
  u_M[,1,] <- s_M #start at activity center
  
  for(k in 2:(14*nFortnight*12)){
    for(i in 1:totZ_F){
      u_F[i,k,1] <- rnorm(1, s_F[i,1], sd = sig[1])
      u_F[i,k,2] <- rnorm(1, s_F[i,2], sd = sig[1])
    }
    for(i in 1:totZ_M){
      u_M[i,k,1] <- rnorm(1, s_M[i,1], sd = sig[2])
      u_M[i,k,2] <- rnorm(1, s_M[i,2], sd = sig[2])
    }
    
  } #end k
  
  plot(scaledObjects$coordsDataScaled, cex = .5, pch =4, ylim = ybounds, xlim = xbounds)
  points(u_F[1,,1],u_F[1,,2], type = 'b', cex = .2, col ='brown')
  points(u_M[1,,1],u_M[1,,2], type = 'b', cex = .2, col ='blue')
  
  ## camera capture:
  nsum_F <- array(NA, c(ntraps,nFortnight))
  lamunknown_F <- array(NA, c(totZ_F, ntraps))
  
  nday_F <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_F <- array(0, c(tagged_F, ntraps,nFortnight)) #det history
  
  
  nsum_M <- array(NA, c(ntraps,nFortnight))
  lamunknown_M <- array(NA, c(totZ_M, ntraps))
  nday_M <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_M <- array(0, c(tagged_M, ntraps,nFortnight))
  
  for(j in 1:ntraps){
    for(i in 1:totZ_F){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_F[i,1] 
      dy <- scaledObjects$coordsDataScaled[j,2] - s_F[i,2]
      d2 <- dx^2 + dy^2 #distance trap to activity center
      lamunknown_F[i,j] <- lam0[1]*exp(-d2/(2*sig[1]*sig[1])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_F[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_F[i,j])) #was I seen each day?
      }
    } #end i
    for(i in 1:totZ_M){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_M[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_M[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_M[i,j] <- lam0[2]*exp(-d2/(2*sig[2]*sig[2])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_M[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_M[i,j])) #was I seen each day?
      }
    } #end i
    for(t in 1:nFortnight){
      nsum_F[j,t] <- sum(colSums(nday_F[,j,,t] >0)) #total days in a 2-week time period when at least 1 turkey detected
      nsum_M[j,t] <- sum(colSums(nday_M[,j,,t] >0)) 
      for(ii in 1:tagged_F){
        y_F[ii,j,t] <- sum(nday_F[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      }
      for(ii in 1:tagged_M){
        y_M[ii,j,t] <- sum(nday_M[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      } #end ii
    } #end t
  } #j
  
  M <- 700
  u <- array(NA, c(tagged_F+tagged_M,14*nFortnight*12,2))
  u[1:tagged_F,,] <- u_F[1:tagged_F,,]
  u[(tagged_F+1):(tagged_M+tagged_F),,] <-   u_M[1:tagged_M,,]
  sex.u <- c(rep(0, tagged_F), rep(1, tagged_M))
  y <- array(0, c(M, ntraps, nFortnight))
  sex.y <- c(rep(0, tagged_F), rep(1, tagged_M),  #real ones
             rep(0, totZ_F-tagged_F+100),rep(1, totZ_M-tagged_M+100)) #augment up to 100 extra of each sex
  y[1:tagged_F,,] <- y_F[1:tagged_F,,]
  y[c(tagged_F+1):(tagged_F+tagged_M),,] <- y_M[1:tagged_M,,]
  
  z.init <- rbinom(M, 1, .5)
  z.init[1:(tagged_F+tagged_M)] <- 1
  
  chooseme <- sample(1:nrow(scaledObjects$coordsHabitatGridCenterScaled), M, replace = T)
  s.init <- scaledObjects$coordsHabitatGridCenterScaled[chooseme,]
  s.init[1:tagged_F,] <-  s_F[1:tagged_F,]
  s.init[c((tagged_F+1):(tagged_F+tagged_M)),] <- s_M[1:tagged_M,]
  
  nsum <- array(0, c(ntraps, nFortnight, 2)) #combine info from both sexes
  nsum[,,1] <- nsum_F
  nsum[,,2] <- nsum_M
  
  ### Prep nimble objects 
  y.full <- getSparseY(y)
  detindices <- y.full$detIndices
  detNums <- y.full$detNums
  trials <- array(14, c(ntraps, nFortnight))
  trapcoords <- scaledObjects$coordsDataScaled
  nimdat <- list(trapCoords = trapcoords,
                 trials = trials,
                 y = y.full$y,
                 detNums=detNums,
                 detIndices=detindices,
                 habLoCoords = lowerAndUpperCoords$lowerHabCoords,
                 habUpCoords = lowerAndUpperCoords$upperHabCoords,
                 habitatGrid = lowerAndUpperCoords$habitatGrid,
                 nsum = nsum,
                 u = u
  ) 
  
  nimconst <- list(
    n.traps = ntraps,
    y.max = dim(habitatMask)[1],
    x.max = dim(habitatMask)[2],
    y.maxDet = dim(trapLocal$habitatGrid)[1], 
    x.maxDet = dim(trapLocal$habitatGrid)[2], 
    resizeFactor = trapLocal$resizeFactor, 
    n.cells = dim(trapLocal$localIndices)[1], 
    maxNBDets = trapLocal$numLocalIndicesMax, 
    trapIndex = trapLocal$localIndices, 
    nTraps = trapLocal$numLocalIndices, 
    habitatIDDet = trapLocal$habitatGrid,
    lengthYCombined = y.full$lengthYCombined,
    numHabWindows = dim(lowerAndUpperCoords$lowerHabCoords)[1],
    nFortnight = nFortnight,
    nind = tagged_F+tagged_M,
    res_r = res(r)[1],
    y_dim = dim(y.full$y)[2],
    M = M,
    sex.u = sex.u,
    sex.y = sex.y
  )
  
  niminits <- list(s = s.init,
                   z = z.init,
                   lam0 = lam0,
                   sigma =sig,
                   psi = rep(.5, nFortnight)
  )
  
  ## also need to know the true N in the habitat area:
  sIDF <- array(NA, totZ_F)
  for(i in 1:totZ_F){
    sIDF[i] <- trapLocal$habitatGrid[trunc(s_F[i,2])+1, trunc(s_F[i,1])+1]
  }
  sIDM <- array(NA, totZ_M)
  for(i in 1:totZ_M){
    sIDM[i] <- trapLocal$habitatGrid[trunc(s_M[i,2])+1, trunc(s_M[i,1])+1]
  }
  Nf_true <- sum(sIDF >0)
  Nm_true <- sum(sIDM >0)
  
  print(paste0('Done with simulation!'))
  return(list(nimconst = nimconst,
              niminits = niminits,
              nimdat = nimdat,
              Nf = totZ_F,
              Nf = totZ_M,
              Nf_true = Nf_true,
              Nm_true = Nm_true,
              sig = sig,
              lam0 = lam0))
}


set.seed(6)
sim_list <- replicate(
  100,
  simpop(),
  simplify = FALSE
)

#saveRDS(sim_list, 'simulation_forSMR.rds')

### simulation #2 - uniform trap grid, 1 km ####
### Uniform trap array at 1 km spacing, but keeping 100 traps. 

#Possible traps:
grid <- rast(ext = ext(r), resolution = 1000, crs = crs(r))
site_bounds <-rbind(BFG, CC)
pts <- as.points(grid) #keep the locations of the grid we made
pts_in <- pts[site_bounds, ]
plot(site_bounds)
points(pts_in[1:58, ], col = 'blue')
points(pts_in[59:nrow(pts_in), ], col = 'darkgreen')

### This is essentially the same function, just allows the trap array to change slightly (picking 100 random traps of the uniform array) each time:
simpop2 <- function(sig = c(.6, .9), #female and male sigmas (scaled for nimble)
                   lam0 = c(.07, .03),
                   tagged_F = 30,
                   tagged_M = 15 #female and male lam0
){ 
  #we need to make sure both sites are sampled:
  BFG_traps <- geom(pts_in[sample(1:58, 25),])[,3:4] #trap coords 
  CC_traps <- geom(pts_in[sample(59:nrow(pts_in), 75),])[,3:4] 
  UniformTraps <- rbind(BFG_traps, CC_traps)
  colnames(UniformTraps) <- c('x', 'y')
  ntraps <- nrow(UniformTraps) #should always be 100
  # plot(r,legend=F)
  # plot(site_bounds, add = T)
  # points(UniformTraps,
  #        pch = 16, col="red") 
  
  habitatMask <- as.matrix(r, wide=TRUE)
  
  ## Rescale coordinates
  scaledObjects <- scaleCoordsToHabitatGrid(
    coordsData = UniformTraps,
    coordsHabitatGridCenter = coordsHabitatGridCenter)
  
  ## Get lower and upper cell coordinates (flips raster but that's normal)
  lowerAndUpperCoords <- getWindowCoords(
    scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
    scaledObsGridCenter = scaledObjects$coordsDataScaled,
    plot.check = F)
  
  dmax = 25
  trapLocal <- getLocalObjects(habitatMask = habitatMask,
                               coords = scaledObjects$coordsDataScaled,
                               dmax = dmax,
                               resizeFactor = 1,
                               plot.check = F)
  
  
  totZ_M <- 200 # true M pop
  totZ_F <- 300 # true F pop
  
  xbounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,1])
  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2])
  
  ## For plot later:
  #s1_array <- list(xbounds = range(scaledObjects$coordsHabitatGridCenterScaled[,1]),
  #                 ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2]),
  #                 traps = scaledObjects$coordsDataScaled)
  
  
  s_F <- cbind(runif(totZ_F, xbounds[1], xbounds[2]), 
               runif(totZ_F, ybounds[1], ybounds[2]))
  s_M <- cbind(runif(totZ_M, xbounds[1], xbounds[2]), 
               runif(totZ_M, ybounds[1], ybounds[2]))
  
  ## we do need our tagged friends to exist on the habitat grid though, so restrict the x values a bit:
  s_F[1:tagged_F,1] <- runif(tagged_F, 8, 20)
  s_M[1:tagged_M,1] <- runif(tagged_M, 8, 20)
  
  
  ## actual locations ###
  u_F <- array(NA, c(totZ_F, nFortnight*14*12, 2)) #14 days in a week, 12 locations per day
  u_M <- array(NA, c(totZ_M, nFortnight*14*12, 2)) 
  
  u_F[,1,] <- s_F #start at activity center
  u_M[,1,] <- s_M #start at activity center
  
  for(k in 2:(14*nFortnight*12)){
    for(i in 1:totZ_F){
      u_F[i,k,1] <- rnorm(1, s_F[i,1], sd = sig[1])
      u_F[i,k,2] <- rnorm(1, s_F[i,2], sd = sig[1])
    }
    for(i in 1:totZ_M){
      u_M[i,k,1] <- rnorm(1, s_M[i,1], sd = sig[2])
      u_M[i,k,2] <- rnorm(1, s_M[i,2], sd = sig[2])
    }
    
  } #end k
  
  plot(scaledObjects$coordsDataScaled, cex = .5, pch =4, ylim = ybounds, xlim = xbounds)
  points(u_F[1,,1],u_F[1,,2], type = 'b', cex = .2, col ='brown')
  points(u_M[1,,1],u_M[1,,2], type = 'b', cex = .2, col ='blue')
  
  ## camera capture:
  nsum_F <- array(NA, c(ntraps,nFortnight))
  lamunknown_F <- array(NA, c(totZ_F, ntraps))
  
  nday_F <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_F <- array(0, c(tagged_F, ntraps,nFortnight)) #det history
  
  
  nsum_M <- array(NA, c(ntraps,nFortnight))
  lamunknown_M <- array(NA, c(totZ_M, ntraps))
  nday_M <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_M <- array(0, c(tagged_M, ntraps,nFortnight))
  
  for(j in 1:ntraps){
    for(i in 1:totZ_F){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_F[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_F[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_F[i,j] <- lam0[1]*exp(-d2/(2*sig[1]*sig[1])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_F[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_F[i,j])) #was I seen each day?
      }
    } #end i
    for(i in 1:totZ_M){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_M[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_M[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_M[i,j] <- lam0[2]*exp(-d2/(2*sig[2]*sig[2])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_M[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_M[i,j])) #was I seen each day?
      }
    } #end i
    for(t in 1:nFortnight){
      nsum_F[j,t] <- sum(colSums(nday_F[,j,,t] >0)) #total days in a 2-week time period when at least 1 turkey detected
      nsum_M[j,t] <- sum(colSums(nday_M[,j,,t] >0)) 
      for(ii in 1:tagged_F){
        y_F[ii,j,t] <- sum(nday_F[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      }
      for(ii in 1:tagged_M){
        y_M[ii,j,t] <- sum(nday_M[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      } #end ii
    } #end t
  } #j
  
  M <- 1000
  u <- array(NA, c(tagged_F+tagged_M,14*nFortnight*12,2))
  u[1:tagged_F,,] <- u_F[1:tagged_F,,]
  u[(tagged_F+1):(tagged_M+tagged_F),,] <-   u_M[1:tagged_M,,]
  sex.u <- c(rep(0, tagged_F), rep(1, tagged_M))
  y <- array(0, c(M, ntraps, nFortnight))
  sex.y <- c(rep(0, tagged_F), rep(1, tagged_M),  #real ones
             rep(0, totZ_F-tagged_F+250),rep(1, totZ_M-tagged_M+250)) #augment up to 250 extra of each sex
  y[1:tagged_F,,] <- y_F[1:tagged_F,,]
  y[c(tagged_F+1):(tagged_F+tagged_M),,] <- y_M[1:tagged_M,,]
  
  z.init <- rbinom(M, 1, .5)
  z.init[1:(tagged_F+tagged_M)] <- 1
  
  chooseme <- sample(1:nrow(scaledObjects$coordsDataScaled), M, replace = T)
  s.init <- scaledObjects$coordsDataScaled[chooseme,]
  s.init[1:tagged_F,] <-  s_F[1:tagged_F,]
  s.init[c((tagged_F+1):(tagged_F+tagged_M)),] <- s_M[1:tagged_M,]
  
  nsum <- array(0, c(ntraps, nFortnight, 2)) #combine info from both sexes
  nsum[,,1] <- nsum_F
  nsum[,,2] <- nsum_M
  
  ### Prep nimble objects 
  y.full <- getSparseY(y)
  detindices <- y.full$detIndices
  detNums <- y.full$detNums
  trials <- array(14, c(ntraps, nFortnight))
  trapcoords <- scaledObjects$coordsDataScaled
  nimdat <- list(trapCoords = trapcoords,
                 trials = trials,
                 y = y.full$y,
                 detNums=detNums,
                 detIndices=detindices,
                 habLoCoords = lowerAndUpperCoords$lowerHabCoords,
                 habUpCoords = lowerAndUpperCoords$upperHabCoords,
                 habitatGrid = lowerAndUpperCoords$habitatGrid,
                 nsum = nsum,
                 u = u
  ) 
  
  nimconst <- list(
    n.traps = ntraps,
    y.max = dim(habitatMask)[1],
    x.max = dim(habitatMask)[2],
    y.maxDet = dim(trapLocal$habitatGrid)[1], 
    x.maxDet = dim(trapLocal$habitatGrid)[2], 
    resizeFactor = trapLocal$resizeFactor, 
    n.cells = dim(trapLocal$localIndices)[1], 
    maxNBDets = trapLocal$numLocalIndicesMax, 
    trapIndex = trapLocal$localIndices, 
    nTraps = trapLocal$numLocalIndices, 
    habitatIDDet = trapLocal$habitatGrid,
    lengthYCombined = y.full$lengthYCombined,
    numHabWindows = dim(lowerAndUpperCoords$lowerHabCoords)[1],
    nFortnight = nFortnight,
    nind = tagged_F+tagged_M,
    res_r = res(r)[1],
    y_dim = dim(y.full$y)[2],
    M = M,
    sex.u = sex.u,
    sex.y = sex.y
  )
  
  niminits <- list(s = s.init,
                   z = z.init,
                   lam0 = lam0,
                   sigma =sig,
                   psi = rep(.5, nFortnight)
  )
  
  ## also need to know the true N in the habitat area:
  sIDF <- array(NA, totZ_F)
  for(i in 1:totZ_F){
    sIDF[i] <- trapLocal$habitatGrid[trunc(s_F[i,2])+1, trunc(s_F[i,1])+1]
  }
  sIDM <- array(NA, totZ_M)
  for(i in 1:totZ_M){
    sIDM[i] <- trapLocal$habitatGrid[trunc(s_M[i,2])+1, trunc(s_M[i,1])+1]
  }
  Nf_true <- sum(sIDF >0)
  Nm_true <- sum(sIDM >0)
  
  print(paste0('Done with simulation!'))
  return(list(nimconst = nimconst,
              niminits = niminits,
              nimdat = nimdat,
              Nf = totZ_F,
              Nf = totZ_M,
              Nf_true = Nf_true,
              Nm_true = Nm_true,
              sig = sig,
              lam0 = lam0))
}


set.seed(6)
sim_list2 <- replicate(
  100,
  simpop2(),
  simplify = FALSE
)
#saveRDS(sim_list2, 'simulation_forSMR_uniformTrap.rds')

### Simuation #3 - 0.5 km spacing ####
#uniform trap array at 0.5 km spacing, but always keeping 100 traps. 

#Possible traps:
grid <- rast(ext = ext(r), resolution = 500, crs = crs(r))
site_bounds <-rbind(BFG, CC)
pts <- as.points(grid) #keep the locations of the grid we made
pts_in <- pts[site_bounds, ]
plot(site_bounds)
points(pts_in[1:224, ], col = 'blue')
points(pts_in[225:nrow(pts_in), ], col = 'darkgreen')


### This is essentially the same function as before
simpop3 <- function(sig = c(.6, .9), #female and male sigmas (scaled for nimble)
                    lam0 = c(.07, .03),
                    tagged_F = 30,
                    tagged_M = 15 #female and male lam0
){ 
  #we need to make sure both sites are sampled:
  BFG_traps <- geom(pts_in[sample(1:224, 25),])[,3:4] #trap coords 
  CC_traps <- geom(pts_in[sample(225:nrow(pts_in), 75),])[,3:4] 
  UniformTraps <- rbind(BFG_traps, CC_traps)
  colnames(UniformTraps) <- c('x', 'y')
  ntraps <- nrow(UniformTraps) #should always be 100
  # plot(r,legend=F)
  # plot(site_bounds, add = T)
  # points(UniformTraps,
  #        pch = 16, col="red") 
  
  habitatMask <- as.matrix(r, wide=TRUE)
  
  ## Rescale coordinates
  scaledObjects <- scaleCoordsToHabitatGrid(
    coordsData = UniformTraps,
    coordsHabitatGridCenter = coordsHabitatGridCenter)
  
  ## Get lower and upper cell coordinates (flips raster but that's normal)
  lowerAndUpperCoords <- getWindowCoords(
    scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
    scaledObsGridCenter = scaledObjects$coordsDataScaled,
    plot.check = F)
  
  dmax = 25
  trapLocal <- getLocalObjects(habitatMask = habitatMask,
                               coords = scaledObjects$coordsDataScaled,
                               dmax = dmax,
                               resizeFactor = 1,
                               plot.check = F)
  
  
  totZ_M <- 200 # true M pop
  totZ_F <- 300 # true F pop
  
  xbounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,1])
  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2])
  
  ## For plot later:
  # s2_array <- list(xbounds = range(scaledObjects$coordsHabitatGridCenterScaled[,1]),
  #                  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2]),
  #                  traps = scaledObjects$coordsDataScaled)
  
  s_F <- cbind(runif(totZ_F, xbounds[1], xbounds[2]), 
               runif(totZ_F, ybounds[1], ybounds[2]))
  s_M <- cbind(runif(totZ_M, xbounds[1], xbounds[2]), 
               runif(totZ_M, ybounds[1], ybounds[2]))
  
  ## we do need our tagged friends to exist on the habitat grid though, so restrict the x values a bit:
  s_F[1:tagged_F,1] <- runif(tagged_F, 8, 20)
  s_M[1:tagged_M,1] <- runif(tagged_M, 8, 20)
  
  
  ## actual locations ###
  u_F <- array(NA, c(totZ_F, nFortnight*14*12, 2)) #14 days in a week, 12 locations per day
  u_M <- array(NA, c(totZ_M, nFortnight*14*12, 2)) 
  
  u_F[,1,] <- s_F #start at activity center
  u_M[,1,] <- s_M #start at activity center
  
  for(k in 2:(14*nFortnight*12)){
    for(i in 1:totZ_F){
      u_F[i,k,1] <- rnorm(1, s_F[i,1], sd = sig[1])
      u_F[i,k,2] <- rnorm(1, s_F[i,2], sd = sig[1])
    }
    for(i in 1:totZ_M){
      u_M[i,k,1] <- rnorm(1, s_M[i,1], sd = sig[2])
      u_M[i,k,2] <- rnorm(1, s_M[i,2], sd = sig[2])
    }
    
  } #end k
  
  plot(scaledObjects$coordsDataScaled, cex = .5, pch =4, ylim = ybounds, xlim = xbounds)
  points(u_F[1,,1],u_F[1,,2], type = 'b', cex = .2, col ='brown')
  points(u_M[1,,1],u_M[1,,2], type = 'b', cex = .2, col ='blue')
  
  ## camera capture:
  nsum_F <- array(NA, c(ntraps,nFortnight))
  lamunknown_F <- array(NA, c(totZ_F, ntraps))
  
  nday_F <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_F <- array(0, c(tagged_F, ntraps,nFortnight)) #det history
  
  
  nsum_M <- array(NA, c(ntraps,nFortnight))
  lamunknown_M <- array(NA, c(totZ_M, ntraps))
  nday_M <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_M <- array(0, c(tagged_M, ntraps,nFortnight))
  
  for(j in 1:ntraps){
    for(i in 1:totZ_F){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_F[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_F[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_F[i,j] <- lam0[1]*exp(-d2/(2*sig[1]*sig[1])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_F[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_F[i,j])) #was I seen each day?
      }
    } #end i
    for(i in 1:totZ_M){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_M[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_M[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_M[i,j] <- lam0[2]*exp(-d2/(2*sig[2]*sig[2])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_M[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_M[i,j])) #was I seen each day?
      }
    } #end i
    for(t in 1:nFortnight){
      nsum_F[j,t] <- sum(colSums(nday_F[,j,,t] >0)) #total days in a 2-week time period when at least 1 turkey detected
      nsum_M[j,t] <- sum(colSums(nday_M[,j,,t] >0)) 
      for(ii in 1:tagged_F){
        y_F[ii,j,t] <- sum(nday_F[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      }
      for(ii in 1:tagged_M){
        y_M[ii,j,t] <- sum(nday_M[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      } #end ii
    } #end t
  } #j
  
  M <- 1000
  u <- array(NA, c(tagged_F+tagged_M,14*nFortnight*12,2))
  u[1:tagged_F,,] <- u_F[1:tagged_F,,]
  u[(tagged_F+1):(tagged_M+tagged_F),,] <-   u_M[1:tagged_M,,]
  sex.u <- c(rep(0, tagged_F), rep(1, tagged_M))
  y <- array(0, c(M, ntraps, nFortnight))
  sex.y <- c(rep(0, tagged_F), rep(1, tagged_M),  #real ones
             rep(0, totZ_F-tagged_F+250),rep(1, totZ_M-tagged_M+250)) #augment up to 250 extra of each sex
  y[1:tagged_F,,] <- y_F[1:tagged_F,,]
  y[c(tagged_F+1):(tagged_F+tagged_M),,] <- y_M[1:tagged_M,,]
  
  z.init <- rbinom(M, 1, .5)
  z.init[1:(tagged_F+tagged_M)] <- 1
  
  chooseme <- sample(1:nrow(scaledObjects$coordsDataScaled), M, replace = T)
  s.init <- scaledObjects$coordsDataScaled[chooseme,]
  s.init[1:tagged_F,] <-  s_F[1:tagged_F,]
  s.init[c((tagged_F+1):(tagged_F+tagged_M)),] <- s_M[1:tagged_M,]
  
  nsum <- array(0, c(ntraps, nFortnight, 2)) #combine info from both sexes
  nsum[,,1] <- nsum_F
  nsum[,,2] <- nsum_M
  
  ### Prep nimble objects 
  y.full <- getSparseY(y)
  detindices <- y.full$detIndices
  detNums <- y.full$detNums
  trials <- array(14, c(ntraps, nFortnight))
  trapcoords <- scaledObjects$coordsDataScaled
  nimdat <- list(trapCoords = trapcoords,
                 trials = trials,
                 y = y.full$y,
                 detNums=detNums,
                 detIndices=detindices,
                 habLoCoords = lowerAndUpperCoords$lowerHabCoords,
                 habUpCoords = lowerAndUpperCoords$upperHabCoords,
                 habitatGrid = lowerAndUpperCoords$habitatGrid,
                 nsum = nsum,
                 u = u
  ) 
  
  nimconst <- list(
    n.traps = ntraps,
    y.max = dim(habitatMask)[1],
    x.max = dim(habitatMask)[2],
    y.maxDet = dim(trapLocal$habitatGrid)[1], 
    x.maxDet = dim(trapLocal$habitatGrid)[2], 
    resizeFactor = trapLocal$resizeFactor, 
    n.cells = dim(trapLocal$localIndices)[1], 
    maxNBDets = trapLocal$numLocalIndicesMax, 
    trapIndex = trapLocal$localIndices, 
    nTraps = trapLocal$numLocalIndices, 
    habitatIDDet = trapLocal$habitatGrid,
    lengthYCombined = y.full$lengthYCombined,
    numHabWindows = dim(lowerAndUpperCoords$lowerHabCoords)[1],
    nFortnight = nFortnight,
    nind = tagged_F+tagged_M,
    res_r = res(r)[1],
    y_dim = dim(y.full$y)[2],
    M = M,
    sex.u = sex.u,
    sex.y = sex.y
  )
  
  niminits <- list(s = s.init,
                   z = z.init,
                   lam0 = lam0,
                   sigma =sig,
                   psi = rep(.5, nFortnight)
  )
  
  ## also need to know the true N in the habitat area:
  sIDF <- array(NA, totZ_F)
  for(i in 1:totZ_F){
    sIDF[i] <- trapLocal$habitatGrid[trunc(s_F[i,2])+1, trunc(s_F[i,1])+1]
  }
  sIDM <- array(NA, totZ_M)
  for(i in 1:totZ_M){
    sIDM[i] <- trapLocal$habitatGrid[trunc(s_M[i,2])+1, trunc(s_M[i,1])+1]
  }
  Nf_true <- sum(sIDF >0)
  Nm_true <- sum(sIDM >0)
  
  print(paste0('Done with simulation!'))
  return(list(nimconst = nimconst,
              niminits = niminits,
              nimdat = nimdat,
              Nf = totZ_F,
              Nf = totZ_M,
              Nf_true = Nf_true,
              Nm_true = Nm_true,
              sig = sig,
              lam0 = lam0))
}


set.seed(6)
sim_list3 <- replicate(
  100,
  simpop3(),
  simplify = FALSE
)
#saveRDS(sim_list3, 'simulation_forSMR_uniformTrap2.rds')

### Simuation #4 - 250km spacing, clusters ####

sample_centers <- function(xy, ncenters, min_dist){
  chosen <- integer(0)
  candidates <- sample(seq_len(nrow(xy)))
  while(length(chosen) < ncenters){
    if(length(candidates) == 0)
      stop("Couldn't place all centers.")
    i <- candidates[1]
    candidates <- candidates[-1]
    if(length(chosen) == 0){
      chosen <- i
    } else {
      d <- sqrt((xy[chosen,1]-xy[i,1])^2 +
                  (xy[chosen,2]-xy[i,2])^2)
      if(all(d >= min_dist))
        chosen <- c(chosen, i)
    }
  }
  chosen
}

sample_clustered_traps <- function(xy,
                                   ncenters,
                                   ntraps,
                                   min_center_dist = 1500,
                                   sigma = 500){
  centers <- sample_centers(xy,
                            ncenters = ncenters,
                            min_dist = min_center_dist)
  # Number of traps per center
  nt <- rep(floor(ntraps/ncenters), ncenters)
  nt[1:(ntraps %% ncenters)] <- nt[1:(ntraps %% ncenters)] + 1
  used <- integer(0)
  for(j in seq_along(centers)){
    d2 <- (xy[,1]-xy[centers[j],1])^2 +
      (xy[,2]-xy[centers[j],2])^2
    w <- exp(-d2/(2*sigma^2))
    # Don't reuse locations
    w[used] <- 0
    w <- w/sum(w)
    used <- c(used,
              sample(seq_len(nrow(xy)),
                     nt[j],
                     prob = w))
  }
  xy[used,]
}

#Possible traps:
# making small resolution for close clustering
grid <- rast(ext = ext(r), resolution = 250, crs = crs(r))
site_bounds <-rbind(BFG, CC)
pts <- as.points(grid) #keep the locations of the grid we made
pts_in <- pts[site_bounds, ]
plot(site_bounds)
points(pts_in[1:914, ], col = 'blue')
points(pts_in[915:nrow(pts_in), ], col = 'darkgreen')


### This is essentially the same function as before
simpop4 <- function(sig = c(.6, .9), #female and male sigmas (scaled for nimble)
                    lam0 = c(.07, .03),
                    tagged_F = 30,
                    tagged_M = 15 #female and male lam0
){ 
  #we need to make sure both sites are sampled:
  BFG_traps <- geom(pts_in[1:914,])[,3:4] #trap coords 
  CC_traps <- geom(pts_in[915:nrow(pts_in),])[,3:4] 

  BFG_t <- sample_clustered_traps(
    BFG_traps,
    ncenters = 3,
    ntraps = 25,
    min_center_dist = 1500,
    sigma = 800
  )
  
  CC_t <- sample_clustered_traps(
    CC_traps,
    ncenters = 6,
    ntraps = 75,
    min_center_dist = 1500,
    sigma = 800
  )
  
  mytraps <- rbind(BFG_t, CC_t) #grab both sites
  colnames(mytraps) <- c('x', 'y')
  ntraps <- nrow(mytraps) #should always be 100
  
  plot(site_bounds)
  points(BFG_t, col = 'blue', pch = 19, cex = .5)
  points(CC_t, col = 'darkgreen', pch = 19, cex = .5)
  
  habitatMask <- as.matrix(r, wide=TRUE)
  
  ## Rescale coordinates
  scaledObjects <- scaleCoordsToHabitatGrid(
    coordsData = mytraps,
    coordsHabitatGridCenter = coordsHabitatGridCenter)
  
  ## Get lower and upper cell coordinates (flips raster but that's normal)
  lowerAndUpperCoords <- getWindowCoords(
    scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
    scaledObsGridCenter = scaledObjects$coordsDataScaled,
    plot.check = F)
  
  dmax = 25
  trapLocal <- getLocalObjects(habitatMask = habitatMask,
                               coords = scaledObjects$coordsDataScaled,
                               dmax = dmax,
                               resizeFactor = 1,
                               plot.check = F)
  
  
  totZ_M <- 200 # true M pop
  totZ_F <- 300 # true F pop
  
  xbounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,1])
  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2])
  
  
  ## For plot later:
  # s3_array <- list(xbounds = range(scaledObjects$coordsHabitatGridCenterScaled[,1]),
  #                  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2]),
  #                  traps = scaledObjects$coordsDataScaled)
  
  s_F <- cbind(runif(totZ_F, xbounds[1], xbounds[2]), 
               runif(totZ_F, ybounds[1], ybounds[2]))
  s_M <- cbind(runif(totZ_M, xbounds[1], xbounds[2]), 
               runif(totZ_M, ybounds[1], ybounds[2]))
  
  ## we do need our tagged friends to exist on the habitat grid though, so restrict the x values a bit:
  s_F[1:tagged_F,1] <- runif(tagged_F, 8, 20)
  s_M[1:tagged_M,1] <- runif(tagged_M, 8, 20)
  
  
  ## actual locations ###
  u_F <- array(NA, c(totZ_F, nFortnight*14*12, 2)) #14 days in a week, 12 locations per day
  u_M <- array(NA, c(totZ_M, nFortnight*14*12, 2)) 
  
  u_F[,1,] <- s_F #start at activity center
  u_M[,1,] <- s_M #start at activity center
  
  for(k in 2:(14*nFortnight*12)){
    for(i in 1:totZ_F){
      u_F[i,k,1] <- rnorm(1, s_F[i,1], sd = sig[1])
      u_F[i,k,2] <- rnorm(1, s_F[i,2], sd = sig[1])
    }
    for(i in 1:totZ_M){
      u_M[i,k,1] <- rnorm(1, s_M[i,1], sd = sig[2])
      u_M[i,k,2] <- rnorm(1, s_M[i,2], sd = sig[2])
    }
    
  } #end k
  
  plot(scaledObjects$coordsDataScaled, cex = .5, pch =4, ylim = ybounds, xlim = xbounds)
  points(u_F[1,,1],u_F[1,,2], type = 'b', cex = .2, col ='brown')
  points(u_M[1,,1],u_M[1,,2], type = 'b', cex = .2, col ='blue')
  
  ## camera capture:
  nsum_F <- array(NA, c(ntraps,nFortnight))
  lamunknown_F <- array(NA, c(totZ_F, ntraps))
  
  nday_F <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_F <- array(0, c(tagged_F, ntraps,nFortnight)) #det history
  
  
  nsum_M <- array(NA, c(ntraps,nFortnight))
  lamunknown_M <- array(NA, c(totZ_M, ntraps))
  nday_M <- array(0, c(totZ_F, ntraps, 14, nFortnight)) #14 days in a week
  y_M <- array(0, c(tagged_M, ntraps,nFortnight))
  
  for(j in 1:ntraps){
    for(i in 1:totZ_F){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_F[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_F[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_F[i,j] <- lam0[1]*exp(-d2/(2*sig[1]*sig[1])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_F[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_F[i,j])) #was I seen each day?
      }
    } #end i
    for(i in 1:totZ_M){
      dx <- scaledObjects$coordsDataScaled[j,1] - s_M[i,1]
      dy <- scaledObjects$coordsDataScaled[j,2] - s_M[i,2]
      d2 <- dx^2 + dy^2
      lamunknown_M[i,j] <- lam0[2]*exp(-d2/(2*sig[2]*sig[2])) #each individual for this trap
      for(t in 1:nFortnight){
        nday_M[i,j,,t] <- rbinom(14, 1, p = 1-exp(-lamunknown_M[i,j])) #was I seen each day?
      }
    } #end i
    for(t in 1:nFortnight){
      nsum_F[j,t] <- sum(colSums(nday_F[,j,,t] >0)) #total days in a 2-week time period when at least 1 turkey detected
      nsum_M[j,t] <- sum(colSums(nday_M[,j,,t] >0)) 
      for(ii in 1:tagged_F){
        y_F[ii,j,t] <- sum(nday_F[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      }
      for(ii in 1:tagged_M){
        y_M[ii,j,t] <- sum(nday_M[ii,j,,t]) #total times seen this occasion (max = ndays) by this trap
        #can only see the tagged ones
      } #end ii
    } #end t
  } #j
  
  M <- 1000
  u <- array(NA, c(tagged_F+tagged_M,14*nFortnight*12,2))
  u[1:tagged_F,,] <- u_F[1:tagged_F,,]
  u[(tagged_F+1):(tagged_M+tagged_F),,] <-   u_M[1:tagged_M,,]
  sex.u <- c(rep(0, tagged_F), rep(1, tagged_M))
  y <- array(0, c(M, ntraps, nFortnight))
  sex.y <- c(rep(0, tagged_F), rep(1, tagged_M),  #real ones
             rep(0, totZ_F-tagged_F+250),rep(1, totZ_M-tagged_M+250)) #augment up to 250 extra of each sex
  y[1:tagged_F,,] <- y_F[1:tagged_F,,]
  y[c(tagged_F+1):(tagged_F+tagged_M),,] <- y_M[1:tagged_M,,]
  
  z.init <- rbinom(M, 1, .5)
  z.init[1:(tagged_F+tagged_M)] <- 1
  
  chooseme <- sample(1:nrow(scaledObjects$coordsDataScaled), M, replace = T)
  s.init <- scaledObjects$coordsDataScaled[chooseme,]
  s.init[1:tagged_F,] <-  s_F[1:tagged_F,]
  s.init[c((tagged_F+1):(tagged_F+tagged_M)),] <- s_M[1:tagged_M,]
  
  nsum <- array(0, c(ntraps, nFortnight, 2)) #combine info from both sexes
  nsum[,,1] <- nsum_F
  nsum[,,2] <- nsum_M
  
  ### Prep nimble objects 
  y.full <- getSparseY(y)
  detindices <- y.full$detIndices
  detNums <- y.full$detNums
  trials <- array(14, c(ntraps, nFortnight))
  trapcoords <- scaledObjects$coordsDataScaled
  nimdat <- list(trapCoords = trapcoords,
                 trials = trials,
                 y = y.full$y,
                 detNums=detNums,
                 detIndices=detindices,
                 habLoCoords = lowerAndUpperCoords$lowerHabCoords,
                 habUpCoords = lowerAndUpperCoords$upperHabCoords,
                 habitatGrid = lowerAndUpperCoords$habitatGrid,
                 nsum = nsum,
                 u = u
  ) 
  
  nimconst <- list(
    n.traps = ntraps,
    y.max = dim(habitatMask)[1],
    x.max = dim(habitatMask)[2],
    y.maxDet = dim(trapLocal$habitatGrid)[1], 
    x.maxDet = dim(trapLocal$habitatGrid)[2], 
    resizeFactor = trapLocal$resizeFactor, 
    n.cells = dim(trapLocal$localIndices)[1], 
    maxNBDets = trapLocal$numLocalIndicesMax, 
    trapIndex = trapLocal$localIndices, 
    nTraps = trapLocal$numLocalIndices, 
    habitatIDDet = trapLocal$habitatGrid,
    lengthYCombined = y.full$lengthYCombined,
    numHabWindows = dim(lowerAndUpperCoords$lowerHabCoords)[1],
    nFortnight = nFortnight,
    nind = tagged_F+tagged_M,
    res_r = res(r)[1],
    y_dim = dim(y.full$y)[2],
    M = M,
    sex.u = sex.u,
    sex.y = sex.y
  )
  
  niminits <- list(s = s.init,
                   z = z.init,
                   lam0 = lam0,
                   sigma =sig,
                   psi = rep(.5, nFortnight)
  )
  
  ## also need to know the true N in the habitat area:
  sIDF <- array(NA, totZ_F)
  for(i in 1:totZ_F){
    sIDF[i] <- trapLocal$habitatGrid[trunc(s_F[i,2])+1, trunc(s_F[i,1])+1]
  }
  sIDM <- array(NA, totZ_M)
  for(i in 1:totZ_M){
    sIDM[i] <- trapLocal$habitatGrid[trunc(s_M[i,2])+1, trunc(s_M[i,1])+1]
  }
  Nf_true <- sum(sIDF >0)
  Nm_true <- sum(sIDM >0)
  
  print(paste0('Done with simulation!'))
  return(list(nimconst = nimconst,
              niminits = niminits,
              nimdat = nimdat,
              Nf = totZ_F,
              Nf = totZ_M,
              Nf_true = Nf_true,
              Nm_true = Nm_true,
              sig = sig,
              lam0 = lam0))
}


set.seed(6)
sim_list4 <- replicate(
  100,
  simpop4(),
  simplify = FALSE
)
#saveRDS(sim_list4, 'simulation_forSMR_clusters.rds')

### Evaluate results ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)
library(scales)
## Grab all simulations:
sim_list <- readRDS('./Simulation/simulation_forSMR.rds')
sim_list_uniformSpacing <- readRDS('./Simulation/simulation_forSMR_uniformTrap.rds')
sim_list_uniformSmall <- readRDS('./Simulation/simulation_forSMR_uniformTrap2.rds')
sim_list_clustered <- readRDS('./Simulation/simulation_forSMR_clusters.rds')

## Skip to below commented out section to view results. The commented out section just shows how the results were turned from RDS files into a table. 

# SimOuts <- list.files('./SimAbund')
# k_vals <- readr::parse_number(SimOuts)
# n <- length(k_vals) #eventually will be 200
# Model <- ifelse(1:n %in% grep('Bands', SimOuts), 'BandSMR', 'TelemetrySMR')
# which(table(k_vals) < 2)
# 
# sim.res <- data.frame(
#   k     = integer(n),
#   N    = numeric(n),
#   NF    = numeric(n),
#   NM    = numeric(n),
#   SDF    = numeric(n),
#   SDM    = numeric(n),
#   lam01 = numeric(n),
#   lam02 = numeric(n),
#   sig1  = numeric(n),
#   sig2  = numeric(n),
#   SDFS  = numeric(n),
#   SDFL  = numeric(n),
#   SDMS  = numeric(n),
#   SDML  = numeric(n),
#   truesig1 = 0.6,
#   truesig2 = 0.9,
#   truelam01 = .07,
#   truelam02= .03,
#   N1CI  = logical(n),
#   N2CI  = logical(n),
#   Sig1CI = logical(n),
#   Sig2CI = logical(n),
#   Lam01CI = logical(n),
#   Lam02CI = logical(n),
#   TrueNF  = numeric(n),
#   TrueNM = numeric(n),
#   Model = Model
# )
# 
# for (i in 1:n){
#   me <- readRDS(file.path(paste0('./SimAbund/', SimOuts[i])))
# 
#   summ    <- summary(me)
#   sum.me  <- summ$quantiles
#   stats.me <- summ$statistics
# 
#   sim.res$k[i]     <- k_vals[i]
#   sim.res$N[i]    <- sum.me[1, 3] #total abund
#   sim.res$NF[i]    <- sum.me[2, 3] #female abund
#   sim.res$SDF[i]    <- stats.me[2, 2] #SD of abund
#   sim.res$SDM[i]    <-  stats.me[3, 2]  #male abund
#   sim.res$NM[i]    <- sum.me[3, 3]
#   sim.res$lam01[i] <- sum.me[4, 3] #lam0
#   sim.res$lam02[i] <- sum.me[5, 3]
#   sim.res$sig1[i]  <- sum.me[8, 3] #sig
#   sim.res$sig2[i]  <- sum.me[9, 3]
#   sim.res$SDFS[i]    <- stats.me[8, 2]
#   sim.res$SDFL[i]    <-  stats.me[4, 2]
#   sim.res$SDMS[i]    <- stats.me[9, 2]
#   sim.res$SDML[i]    <-  stats.me[5, 2]
# 
#   myconsts <- c(sim_list[k_vals[i]][[1]]$Nf_true, sim_list[k_vals[i]][[1]]$Nm_true)
# 
#   sim.res$N1CI[i] <- sum.me[2, 1] <= myconsts[1]  & sum.me[2, 5] >=  myconsts[1] #CI performance
#   sim.res$N2CI[i] <- sum.me[3, 1] <=  myconsts[2] & sum.me[3, 5] >=  myconsts[2]
#   sim.res$Lam01CI[i] <- sum.me[4, 1] <= .07  & sum.me[4, 5] >=  .07 #CI performance
#   sim.res$Lam02CI[i] <- sum.me[5, 1] <=  .03 & sum.me[5, 5] >=  .03
#   sim.res$Sig1CI[i] <- sum.me[8, 1] <= 0.6  & sum.me[8, 5] >=  0.6
#   sim.res$Sig2CI[i] <- sum.me[9, 1] <=  .9 & sum.me[9, 5] >=  .9
# 
#   sim.res$TrueNF[i] <- myconsts[1]
#   sim.res$TrueNM[i] <- myconsts[2]
# }

SimOuts_Uniform <- list.files('./SimAbund/SMR_2')
#k_vals <- readr::parse_number(SimOuts_Uniform)
k_vals <- as.integer(sub(".*_(\\d+)\\.rds$", "\\1", SimOuts_Uniform))
n <- length(k_vals) #eventually will be 200
Model <- ifelse(1:n %in% grep('Bands', SimOuts_Uniform), 'BandSMR', 'TelemetrySMR')
#which(table(k_vals) < 2)

sim.resUni <- data.frame(
  k     = integer(n),
  N    = numeric(n),
  NF    = numeric(n),
  NM    = numeric(n),
  SDF    = numeric(n),
  SDM    = numeric(n),
  lam01 = numeric(n),
  lam02 = numeric(n),
  sig1  = numeric(n),
  sig2  = numeric(n),
  SDFS  = numeric(n),
  SDFL  = numeric(n),
  SDMS  = numeric(n),
  SDML  = numeric(n),
  truesig1 = 0.6,
  truesig2 = 0.9,
  truelam01 = .07,
  truelam02= .03,
  N1CI  = logical(n),
  N2CI  = logical(n),
  Sig1CI = logical(n),
  Sig2CI = logical(n),
  Lam01CI = logical(n),
  Lam02CI = logical(n),
  TrueNF  = numeric(n),
  TrueNM = numeric(n),
  Model = Model,
  Type = NA
)

for (i in 1:n){
  sim.resUni$Type[i] <- ifelse(length(grep('Unipt5', SimOuts_Uniform[i])) == 1, '0.5-km Uniform', 
                               ifelse(length(grep('cluster', SimOuts_Uniform[i])) == 1,
                                      '0.25-km Clustered', '1-km Uniform'))
  
  me <- readRDS(file.path(paste0('./SimAbund/SMR_2/', SimOuts_Uniform[i])))
  
  summ    <- summary(me)
  sum.me  <- summ$quantiles
  stats.me <- summ$statistics
  
  
  sim.resUni$k[i]     <- k_vals[i]
  sim.resUni$N[i]    <- sum.me[1, 3] #total abund
  sim.resUni$NF[i]    <- sum.me[2, 3] #female abund
  sim.resUni$SDF[i]    <- stats.me[2, 2] #SD of abund
  sim.resUni$SDM[i]    <-  stats.me[3, 2]  #female abund
  sim.resUni$NM[i]    <- sum.me[3, 3]
  sim.resUni$lam01[i] <- sum.me[4, 3] #lam0
  sim.resUni$lam02[i] <- sum.me[5, 3]
  sim.resUni$sig1[i]  <- sum.me[8, 3] #sig
  sim.resUni$sig2[i]  <- sum.me[9, 3]
  sim.resUni$SDFS[i]    <- stats.me[8, 2]
  sim.resUni$SDFL[i]    <-  stats.me[4, 2]
  sim.resUni$SDMS[i]    <- stats.me[9, 2]
  sim.resUni$SDML[i]    <-  stats.me[5, 2]
  
 if(sim.resUni$Type[i] == '1-km Uniform'){
  myconsts <- c(sim_list_uniformSpacing[k_vals[i]][[1]]$Nf_true, sim_list_uniformSpacing[k_vals[i]][[1]]$Nm_true)
 }
  if(sim.resUni$Type[i] == '0.5-km Uniform'){
    myconsts <- c(sim_list_uniformSmall[k_vals[i]][[1]]$Nf_true, sim_list_uniformSmall[k_vals[i]][[1]]$Nm_true)
  }
  if(sim.resUni$Type[i] == '0.25-km Clustered'){
    myconsts <- c(sim_list_clustered[k_vals[i]][[1]]$Nf_true, sim_list_clustered[k_vals[i]][[1]]$Nm_true)
  }

  sim.resUni$N1CI[i] <- sum.me[2, 1] <= myconsts[1]  & sum.me[2, 5] >=  myconsts[1] #CI performance
  sim.resUni$N2CI[i] <- sum.me[3, 1] <=  myconsts[2] & sum.me[3, 5] >=  myconsts[2]
  sim.resUni$Lam01CI[i] <- sum.me[4, 1] <= .07  & sum.me[4, 5] >=  .07 #CI performance
  sim.resUni$Lam02CI[i] <- sum.me[5, 1] <=  .03 & sum.me[5, 5] >=  .03
  sim.resUni$Sig1CI[i] <- sum.me[8, 1] <= 0.6  & sum.me[8, 5] >=  0.6
  sim.resUni$Sig2CI[i] <- sum.me[9, 1] <=  .9 & sum.me[9, 5] >=  .9

  sim.resUni$TrueNF[i] <- myconsts[1]
  sim.resUni$TrueNM[i] <- myconsts[2]
}

### This is the object produced by code above:
# sim.res <- read.csv('./Simulation/simulation.results.csv')
# sim.res$Type <- 'Case Study Array'

sim.plot <- rbind(sim.res[,-1], sim.resUni)
sim.plot$CV.F <- sim.plot$SDF/sim.plot$NF*100
sim.plot$CV.M <- sim.plot$SDM/sim.plot$NM*100

simCIs <- sim.plot %>% 
  group_by(Model, Type) %>%
  summarise(NCI.F = mean(N1CI),
            NCI.M = mean(N2CI),
            SigCIF = mean(Sig1CI),
            SigCIM = mean(Sig2CI),
            LamCIF = mean(Lam01CI),
            LamCIM = mean(Lam02CI),
            SD.F = median(SDF),
            SD.M = median(SDM),
            SD.sigF = median(SDFS),
            SD.lamF = median(SDFL),
            SD.sigM = median(SDMS),
            SD.lamM = median(SDML),
            CV.F = median(CV.F),
            CV.M = median(CV.M),
            n = n())

#CI performance:
simCIs



mysummary <- function(x, y){
  a <- median(x-y)
  b <- sd(x-y)
  c <- median((x-y)/y)
  d <- sd((x-y)/y)
  return(round(c(a,b,c,d), 2))
}

BandRes <- subset(sim.plot, sim.plot$Model == "BandSMR" & sim.plot$Type == 'Case Study Array')
TelemRes <- subset(sim.plot, sim.plot$Model == 'TelemetrySMR' & sim.plot$Type == 'Case Study Array')
BandRes_U1 <- subset(sim.plot, sim.plot$Model == "BandSMR" & sim.plot$Type != "1-km Uniform" )
TelemRes_U1 <- subset(sim.plot, sim.plot$Model == 'TelemetrySMR' & sim.plot$Type != "1-km Uniform" )
BandRes_U2 <- subset(sim.plot, sim.plot$Model == "BandSMR" & sim.plot$Type != "0.5-km Uniform" )
TelemRes_U2 <- subset(sim.plot, sim.plot$Model == 'TelemetrySMR' & sim.plot$Type != "0.5-km Uniform" )
BandRes_clust <- subset(sim.plot, sim.plot$Model == "BandSMR" & sim.plot$Type != "0.25-km Clustered" )
TelemRes_clust <- subset(sim.plot, sim.plot$Model == 'TelemetrySMR' & sim.plot$Type != "0.25-km Clustered" )


quick_summary <- function(xx){
  NF_a <- mysummary(x = xx$NF, y = xx$TrueNF) #median bias/sd; relative bias/sd
  NM_a <- mysummary(x = xx$NM, y = xx$TrueNM) 
  SigM_a <- mysummary(x = xx$sig2, y = xx$truesig2) 
  SigF_a <- mysummary(x = xx$sig1, y = xx$truesig1) 
  lamF_a <- mysummary(x = xx$lam01, y = xx$truelam01) 
  lamM_a <- mysummary(x = xx$lam02, y = xx$truelam02) 
  aa <- rbind(NF_a, NM_a, lamF_a, lamM_a, SigF_a, SigM_a)
  colnames(aa) <- c('Median Bias', 'Median Bias SD', 'Relative Bias', 'Relative Bias SD')
  return(aa)
}
quick_summary(BandRes)
quick_summary(BandRes_U1)
quick_summary(BandRes_U2)
quick_summary(BandRes_clust)
quick_summary(rbind(BandRes_U1, BandRes_U2, BandRes_clust))
quick_summary(TelemRes)
quick_summary(TelemRes_U1)
quick_summary(TelemRes_U2)
quick_summary(TelemRes_clust)
quick_summary(rbind(TelemRes, TelemRes_U1, TelemRes_U2, TelemRes_clust))

sim.visual <- sim.plot %>% 
  mutate(
    bias_NF = (NF - TrueNF)/TrueNF*100,
    bias_NM = (NM - TrueNM)/TrueNM*100
  ) %>%
  select(Model, Type, bias_NF, bias_NM) %>%
  pivot_longer(
    cols = starts_with("bias_"),
    names_to = "Group",
    values_to = "Bias"
  ) %>%
  mutate(Group = recode(Group,
                        bias_NF = "Female",
                        bias_NM = "Male"))


aa <- ggplot(sim.visual, aes(x= Group, y = Bias))+
  geom_boxplot(aes(fill = Model), alpha =.5)+
  theme_bw()+
  facet_wrap(~Type, nrow = 1)+
  ylab('% Relative bias')+
  xlab('')+
  scale_y_continuous(
    labels = label_percent(scale = 1)
  )+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_manual(values = c('grey30', '#092050'))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = 'top')

sim.visual2 <- sim.plot %>%
  select(Model, Type, CV.F, CV.M) %>%
  pivot_longer(
    cols = starts_with("CV"),
    names_to = "Group",
    values_to = "Uncertainty"
  ) %>%
  mutate(Group = recode(Group,
                        CV.F = "Female",
                        CV.M = "Male"))

bb <- ggplot(sim.visual2, aes(x= Group, y = Uncertainty))+
  geom_boxplot(aes(fill = Model), alpha =.5, outliers = F)+
  theme_bw()+
  facet_wrap(~Type, nrow = 1)+
  ylab('CV')+
  xlab('')+
  scale_y_continuous(
    labels = label_percent(scale = 1)
  )+
  scale_fill_manual(values = c('grey30', '#092050'))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = 'top')

library(ggpubr)
ggarrange(aa,bb, nrow = 2, common.legend = TRUE, align = "v")
ggsave('Figure_S1.jpeg', width = 11, height = 8)

### Could look at the trap arrays at the top too
#These objects are created in the simulation functions: 
traps.plot <- data.frame(x = c(OG_array$traps[,1],s1_array$traps[,1],s2_array$traps[,1],s3_array$traps[,1]),
                         y = c(OG_array$traps[,2],s1_array$traps[,2],s2_array$traps[,2],s3_array$traps[,2]),
                         array = c(rep('Case Study Array', 98), rep( c("1-km Uniform", "0.5-km Uniform", "0.25-km Clustered"), each = 100)))

cc <- ggplot(traps.plot, aes(x = x, y= y))+
  geom_point()+
  theme_bw()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = 'top')+
  facet_wrap(~array, nrow  = 1)

ggarrange(aa,bb, nrow = 2, common.legend = TRUE)
#ggsave('Figure_S1.jpeg', width = 11, height = 8)


### for manuscript, summarize detections:
det_info <- array(NA, c(100, 4))
uniqueTraps <- function(arr) {
  apply(arr, 1, function(mat) {
    vals <- as.vector(mat)
    vals <- vals[vals != -1]
    length(unique(vals))
  })
}
for(i in 1:100){
  me <- sim_list[[i]]
  medat <- me$nimdat$detNums[1:45,] #all the tagged ones
  seen <- sum(rowSums(medat)>0) #seen at least once by a camera
  unique_cams <- mean(uniqueTraps(me$nimdat$detIndices[1:45,,]))
  avgF <- mean(rowSums(me$nimdat$nsum[,,1]))
  avgM <- mean(rowSums(me$nimdat$nsum[,,2]))
  det_info[i, ] <- c(seen, unique_cams, avgF, avgM)
}

summary(det_info)
