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
WinterCams <- read.csv('WinterCamera_Locations.csv')
BFG <- vect('bf grant boundary.shp')
CC <- vect('cedarcreek_boundary2014.shp')

WinterCams <- subset(WinterCams, WinterCams$CamerID %notin% c('EB01','EB02', 'CC71', 'SA02',  'SG01',  'SG02'))
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

dmax = 25
trapLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = scaledObjects$coordsDataScaled,
                             dmax = dmax,
                             resizeFactor = 1,
                             plot.check = T)

#Simulation function #####
#simpler than real function, assumes population is the same for all 4 fortnights 
simpop <- function(sig = c(.6, .9), #female and male sigmas (scaled for nimble)
                   lam0 = c(.07, .03) #female and male lam0
){ 
  
  totZ_M <- 200 # true M pop
  totZ_F <- 300 # true F pop
  
  tagged_F <- 30 #count of tagged
  tagged_M <- 15 
  
  xbounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,1])
  ybounds <- range(scaledObjects$coordsHabitatGridCenterScaled[,2])
  
  
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


### Evaluate results ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)
# sim_list <- readRDS('simulation_forSMR.rds')
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
#   sim.res$NF[i]    <- sum.me[1, 3] #total abund
#   sim.res$NF[i]    <- sum.me[2, 3] #female abund
#   sim.res$SDF[i]    <- stats.me[2, 2] #SD of abund
#   sim.res$SDM[i]    <-  stats.me[3, 2]  #female abund
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

### This is the object produced by code above:
sim.res <- read.csv('simulation.results.csv')

simCIs <- sim.res %>% 
  group_by(Model) %>%
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

BandRes <- subset(sim.res, sim.res$Model == "BandSMR")
TelemRes <- subset(sim.res, sim.res$Model == 'TelemetrySMR')


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
quick_summary(TelemRes)

sim.visual <- sim.res %>%
  mutate(
    bias_NF = NF - TrueNF,
    bias_NM = NM - TrueNM
  ) %>%
  select(Model, bias_NF, bias_NM) %>%
  pivot_longer(
    cols = starts_with("bias_"),
    names_to = "Group",
    values_to = "Bias"
  ) %>%
  mutate(Group = recode(Group,
                        bias_NF = "Female",
                        bias_NM = "Male"))

ggplot(sim.visual, aes(x= Group, y = Bias))+
  geom_boxplot(aes(fill = Model), alpha =.5)+
  theme_bw()+
  ylab('Bias in abundance')+
  geom_hline(yintercept = 0, lty = 2)+
  scale_fill_manual(values = c('grey30', '#092050'))+
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = 'top')
