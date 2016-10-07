

# this is the function to compute the matern covariance. This shouldn't
# be used directly by the user except in special cases. There are two 
# wrappers, maternIsotropic and maternSphereTime below that take in
# parameters and locations, rather than a distance matrix

# there is also a simulation function at the end


maternFun <- function(distmat, covparms, returnD1 = FALSE, eee = 1e-8 ){
    
    # covariance parameters are...
    # (variance (including nugget), range, smoothness, signal to signal + noise ratio)
    # range of values are...
    # (0,Inf), (0,Inf) , (0,Inf), (0,1) 
    
    n1 <- dim(distmat)[1]
    n2 <- dim(distmat)[2]
    covmat      <-  matrix(0,n1,n2)
    scaledist   <-  distmat/covparms[2]

    # to avoid bad spots for smoothness parameter
    covparms[3] <- min(max(covparms[3],1e-12),20)
    normcon     <-  2^(covparms[3]-1)*gamma(covparms[3])

    # special cases for nu = 1/2 and 3/2
    if( covparms[3] == 1/2 ){ maternpart <- exp(-scaledist) }
    else if( covparms[3] == 3/2 ){ maternpart <- (1 + scaledist)*exp(-scaledist) }
    else {
        besspart <- besselK(scaledist,covparms[3])
        maternpart  <-  1/normcon*(scaledist)^covparms[3]*besspart
    }
    covmat[distmat != 0] <- covparms[1]*covparms[4]*maternpart[distmat!=0]
    covmat[distmat == 0] <- covparms[1]

    
    # if specified, return vector of first derivatives
    if(returnD1){
        
        if(!exists("besspart")) besspart <- besselK(scaledist,covparms[3])
        
        d1 <- covmat/covparms[1]
        
        d2 <- matrix(0,n1,n2)
        bessm1 <- besselK(scaledist,covparms[3]-1)
        bessp1 <- besselK(scaledist,covparms[3]+1)
        d2 <- -1/normcon*covparms[1]*covparms[4]*
            ( covparms[3]*scaledist^(covparms[3]-1)*besspart -
                  scaledist^(covparms[3])*(bessm1 + bessp1)/2 )*
            scaledist/covparms[2]
        d2[distmat==0] <- 0
        
        dcovparms <- covparms
        dcovparms[3] <- covparms[3]+eee
        d3 <- 1/eee*(maternFun(distmat,dcovparms)
                     -maternFun(distmat,covparms))
        
        d4 <- covmat/covparms[4]
        d4[distmat==0] = 0
        
        # not sure if I like this, maybe return a list instead
        attr(covmat,"dmats") <- list( d1=d1, d2=d2, d3=d3, d4=d4 ) 
    }
    covmat
}



maternIsotropic <- function(locs,covparms,returnD1=FALSE){
    
    # use maternFun to compute isotropic matern covariance
    # locs are locations in a Euclidean space.
    distmat <- rdist(locs)
    covmat <- maternFun(distmat,covparms,returnD1)
    
}


maternSphereTime <- function(locstime,covparms,returnD1=FALSE){

    # isotropic matern in space, matern in time
    # locs <- (lon,lat,time)
    # nonseparable
    # covparms = (variance,spatial range,temporal range,smoothness,Sig2Noise)
    
    locs <- locstime[,1:2,drop=FALSE]
    time <- locstime[,3,drop=FALSE]
    chordaldist <- 2*sin(rdist.earth(locs,R=1)/2)
    timedist <- rdist(time)
    distmat <- sqrt( chordaldist^2/covparms[2]^2 + timedist^2/covparms[3]^2 )
    covmat <- maternFun(distmat,c(covparms[1],1,covparms[4],covparms[5]))
}





simulateData <- function(locs, covparms, covfun){
    n <- dim(locs)[1]
    covmat <- covfun( locs, covparms )  
    cholmat <- chol(covmat) 
    y <- c(t(cholmat) %*% rnorm(n))
}


