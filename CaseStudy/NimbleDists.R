library(nimble)
deregisterDistributions('dbinomLocal_normal_p') #just in case
### Function to do local evaluation, assuming normal movement from activity center when y is binary and p is Pr(at least one det)

### this function is almost identical to the one in the nimbleSCR package, with the exception of line 143, which turns detections into a binary format
dbinomLocal_normal_p <- nimbleFunction(
  run = function( x = double(1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(1),
                  p0 = double(0, default = -999),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0),
                  lengthYCombined = double(0, default = 0),
                  log = integer(0, default = 0)
  ) {
    ## Specify return type
    returnType(double(0))
    
    ## Deal with cases where detection info is combined in one vector 
    if(detNums==-999){
      detNums <- x[1]
      nMaxDetectors <- (lengthYCombined-1)/2
      detIndices1 <- x[(nMaxDetectors+2):lengthYCombined]
      x1 <- x[2:(nMaxDetectors+1)]
    }else{
      x1 <- x
      detIndices1 <- detIndices
    }
    
    ## Shortcut if the current individual is not available for detection
    if(indicator == 0){
      if(detNums == 0){
        if(log == 0) return(1.0)
        else return(0.0)
      } else {
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
    
    ## Retrieve the index of the habitat cell where the current AC is
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## Retrieve the indices of the local traps surrounding the selected habita grid cell
    theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
    
    ## CHECK IF DETECTIONS ARE WITHIN THE LIST OF  LOCAL TRAPS
    if(detNums > 0){
      for(r in 1:detNums){
        if(sum(detIndices1[r] == theseLocalTraps) == 0){
          if(log == 0) return(0.0)
          else return(-Inf)
        }
      }
    }
    
    ## Calculate the log-probability of the vector of detections
    alpha <- -1.0 / (2.0 * sigma * sigma)
    logProb <- 0.0 
    detIndices1 <- c(detIndices1, 0)
    count <- 1 
    
    # when p0 is provide through p0
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          lambda <- p0 * exp(alpha * d2)
          psi <- 1 - exp(-lambda)
          logProb <-  logProb + dbinom(x1[count], prob = psi, size = size[theseLocalTraps[r]], log = TRUE)
          count <- count + 1
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          lambda <- p0 * exp(alpha * d2)
          psi <- 1 - exp(-lambda)
          logProb <- logProb + dbinom(0, prob = psi, size = size[theseLocalTraps[r]], log = TRUE)
        }
      }
    
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })

rbinomLocal_normal_p <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(1),
                  p0 = double(0, default = -999),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0),
                  lengthYCombined = double(0, default = 0)
  ) {
    ## Specify return type
    returnType(double(1))
    if(detNums >= 0) stop("Random generation for the rbinomLocal_normal distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
    
    #========================================================
    # RETURN TYPE DECLARATION
    if(n!=1){print("rbinomLocal_normal only allows n = 1; using n = 1")}
    # returnType(double(3))
    # len <- 2*MAX + 1
    ## GET NECESSARY INFO
    alpha <- -1.0 / (2.0 * sigma * sigma)
    # n.detectors <- dim(detector.xy)[1]
    # nMAxDetections <- length(detIndices)
    nMAxDetections <- (lengthYCombined-1)/2
    ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
    #if(indicator == 0){return(rep(0.0, 2*nMAxDetections + 1))}
    if(indicator == 0){return(rep(0.0, lengthYCombined))}
    
    ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
    theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]
    
    ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
    detectOut <- rep(0, localTrapsNum[sID])
    ys <- rep(-1, nMAxDetections)
    dets <- rep(-1, nMAxDetections)
    count <- 1
    
    ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
    ## when p0 is provided through p0
      for(r in 1:localTrapsNum[sID]){
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        lambda <- p0 * exp(alpha * d2)
        psi <- 1 - exp(-lambda)
        # Draw the observation at detector j from a binomial distribution with probability psi
        detectOut[r] <- rbinom(1, size[theseLocalTraps[r]], psi)
        if(detectOut[r] >0){
          if(nMAxDetections<count){stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                          You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
          ys[count] <- detectOut[r]
          dets[count] <- theseLocalTraps[r]
          count <- count + 1
        }#if
      }#r 
      
    count <- count - 1
    
    
    # out <- rep(-1, 2*nMAxDetections + 1)
    out <- rep(-1, lengthYCombined)
    
    out[1] <- count
    if(count >= 1){
      out[2:(count+1)] <- ys[1:count]
      out[(nMAxDetections+2):(nMAxDetections+count+1)] <- dets[1:count]
    }
    ## OUTPUT
    return(out)
  })

### Stage 2

### Function for getting the p needed for the binomial draw of counts at each camera
#
## This is just calculating detection probability in a function, instead of inside the nimble code, to save computation time
## Get_p ####
get_p <- nimbleFunction(
  run = function(s = double(2), 
                 M = double(0),
                 lam0=double(0), 
                 sigma=double(0), 
                 trapCoords = double(2),
                 ntraps=double(0), 
                 z = double(1)){
                 #localTrapsIndices = double(2)){
                 #resizeFactor = double(0, default = 1), 
                 #habitatGrid = double(2)){ 
    returnType(double(1))
    
    #sID <- array(NA, M)
    lam <- array(0, c(M, ntraps))
    Lambda <- array(0, ntraps)
    p <- array(0, ntraps)
    
    # for(i in 1:M){
    # sID[i] <- habitatGrid[trunc(s[i,2]/resizeFactor)+1, trunc(s[i,1]/resizeFactor)+1]
    # }
    
    for(j in 1:ntraps){
      for(i in 1:M){
        if(z[i] == 0){#if not real, p = 0
          lam[i,j] <- 0
        } else{ #if z == 1, calculate p
          d2 <- sqrt((s[i,1]-trapCoords[j,1])^2 + (s[i,2]-trapCoords[j,2])^2)
          lam[i,j] <- lam0*exp(-d2^2/(2*sigma^2))
        }
      } #end i
      Lambda[j] <- sum(lam[1:M,j])
      p[j] <- 1 - exp(-Lambda[j]) #p(at least 1) = 1-exp(-lambda) #from poisson math
    } #end j
    
    return(p)
  }
)


### Function for getting the p needed for the binomial draw of counts at each camera, 2 sexes 
### get p_2sex ####
get_p2sex <- nimbleFunction(
  run = function(s = double(2),  #matrix
                 M = double(0),
                 k = double(0),
                 sex = double(1), #vector
                 lam0=double(1), #vector 
                 sigma=double(1),  #vector
                 trapCoords = double(2), #matrix
                 ntraps=double(0), #number
                 z = double(1)){ #vector
    returnType(double(1))
    
    lam <- array(0, c(M, ntraps))
    Lambda <- array(0, ntraps)
    p <- array(0, ntraps)
    
    for(j in 1:ntraps){
      for(i in 1:M){
        if(z[i] == 0 | sex[i] != (k-1)){#if not real, p = 0
          lam[i,j] <- 0
        } else{ #if z == 1, calculate p
          d2 <- sqrt((s[i,1]-trapCoords[j,1])^2 + (s[i,2]-trapCoords[j,2])^2)
          lam[i,j] <- lam0[k]*exp(-d2^2/(2*sigma[k]^2))
        }
      } #end i
      Lambda[j] <- sum(lam[1:M,j])
      p[j] <- 1 - exp(-Lambda[j]) #p(at least 1) = 1-exp(-lambda) #from poisson math
    } #end j
    
    return(p)
  }
)
