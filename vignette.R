

# a short vignette demonstrating how to use the functions

# required packages
require(fields)          # for the rdist function
require(matrixStats)     # for the colMins function

# R source files
sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source)
#set.seed(1)

# grid size for data locations
gsize <- 40
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
locs <- simulateGrid(nvec,jittersize=0)
plot(locs[,1],locs[,2])

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, sig2noise = 1)

# simulate some data and plot a map of it
# uses Cholesky method so don't try this for large n
y <- simulateData(locs,covparms,covfun)
image( matrix(y,nvec) )

# generate an ordering and plot the first n/8
ord <- orderMaxMinLocal(locs)
n0 <- round(n/8)
plot( locs[ord[1:n0],1],locs[ord[1:n0],2] )

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 30
NNarray <- findOrderedNNfast(locsord,m)

# automatically group the observations
NNlist <- groupNN(NNarray)
NNlist[1:4]  # list elements are subsets of rows of NNarray

# compute exact loglik, ungrouped, and grouped ordered composite logliks
# system.time(  ll0 <- mvnMargLik(covparms,covfun,y,locs) # only do this if n is small 
system.time(  ll1 <- orderedCompLik(covparms,covfun,yord,locsord,NNarray)      )
system.time(  ll2 <- orderedGroupCompLik(covparms,covfun,yord,locsord,NNlist)  )





# an attempt to write a wrapper function to do all of this stuff:
# ordering, finding neighbors, maximizing parameters
# only covariance function implemented is isotropic matern 
# interesting thing: block independent likelihood is ridiculously
# good at finding parameter estimates that nearly maximize
# vecchia's likelihood approximation
# This function uses block independent likelihood to quickly get starting values
# for an optimization of Vecchia's approximation

fitmodel <- function(y,locs,covfun, numneighbors = 30, orderfun = "maxmindist", 
                     grouped = TRUE, fixedparameters = NA ){
    
    # need to add capability for including a linear mean
    n <- length(y)
    
    # check to see if the number of data points matches the 
    # number of spatial locations
    if( nrow(locs) != n ) stop("length of y not equal to number of locations")
    
    # order the points according to orderfun argument
    # only maxmindist implemented. Others are easy too
    if( orderfun == "maxmindist" ){
        ord <- orderMaxMinLocal(locs)
    } else {
        stop("Unrecognized ordering method specified in orderfun argument")
    }
    
    # define link functions for each parameter, i.e. take logs of
    # positive-valued parameters, logit of parameters in (0,1)
    if( identical(covfun,maternIsotropic, ignore.environment = TRUE) ){
        linkfun <- list( function(x) log(x), function(x) log(x), function(x) log(x), function(x) log(x)/log(1-x) )
        invlinkfun <- list(function(x) exp(x),function(x) exp(x),function(x) exp(x), function(x) exp(x)/(1+exp(x)))
        if(identical(NA,fixedparameters)) fixedparameters <- rep(NA,4)
        startvals <- rep(0,4)
    }
    
    # apply the ordering
    yord <- y[ord]
    locsord <- locs[ord,]
    
    # get nearest neighbors and do grouping
    NNarray <- findOrderedNNfast(locsord,numneighbors)
    NNlist <- groupNN(NNarray)
    
    # fixedparameters allows us to fix some of the parameters. any parameter
    # taking on NA in fixed parameters gets estimated. Those with specified
    # values get fixed
    covparms <- fixedparameters
    notfixedinds <- which(is.na(fixedparameters))  # indices of parms to estimate
    
    # block independent approximation to get starting values
    nblocks <- round(n/50)
    nside <- ceiling( sqrt(nblocks) )
    blocklist <- getBlockList(locs,nside,nside)
    
    # independent blocks negative loglikelihood
    f0 <- function(x){
        for(j in 1:length(notfixedinds)){
            covparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](x[j])
        }
        ll <- mvnIndepBlocks(covparms,covfun,y,locs,blocklist)
        return(-ll)
    }
    
    # Vecchia's approximation
    f1 <- function(x){
        for(j in 1:length(notfixedinds)){
            covparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](x[j])
        }
        if(grouped){
            ll <- orderedGroupCompLik(covparms,covfun,yord,locsord,NNlist)
        } else {
            ll <- orderedCompLik(covparms,covfun,yord,locsord,NNarray)
        }
        return(-ll)
    }
    
    # use nelder mead (which seems more stable) to move towards maximum
    # of independent blocks likelihood approximation
    result0 <- optim(startvals[notfixedinds],f0,method="Nelder-Mead",control=list(maxit=50,trace=1))
    # use BFGS to get to the optimum once we are close
    result0 <- optim(result0$par,f0,method="BFGS",control=list(maxit=100,trace=1))
    # use the maximum independent blocks likelihood estimates to start 
    # an optimization of Vecchia's likelihood
    result1 <- optim(result0$par,f1,method="BFGS",control=list(maxit=10,trace=5))
    outparms <- fixedparameters
    # transform back to original parameter domain with inverse link
    for(j in 1:length(notfixedinds)){
        outparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](result1$par[j])
    }
    result1$par <- outparms
    return(result1)
}       
    
# get some paramter estimates with the new function
system.time( result <- fitmodel(y,locs,maternIsotropic,numneighbors=30,fixedparameters=c(NA,NA,NA,1))  )


# see how long it takes to get estimates with exact likelihood    
fixedparameters <- c(NA,NA,NA,1)
notfixedinds <- which(is.na(fixedparameters))  # indices of parms to estimate
linkfun <- list( function(x) log(x), function(x) log(x), function(x) log(x), function(x) log(x)/log(1-x) )
invlinkfun <- list(function(x) exp(x),function(x) exp(x),function(x) exp(x), function(x) exp(x)/(1+exp(x)))
startvals <- rep(0,4)

f2 <- function(x){
    for(j in 1:length(notfixedinds)){
        covparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](x[j])
    }
    ll <- mvnMargLik(covparms,covfun,y,locs)
    return(-ll)
}
system.time({
    result2 <- optim(startvals[notfixedinds],f2,method="Nelder-Mead",control=list(maxit=500,trace=1))
})
outparms <- fixedparameters
# transform back to original parameter domain with inverse link
for(j in 1:length(notfixedinds)){
    outparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](result2$par[j])
}
result2$par <- outparms


# compare the estimates
# approximate methods
result$par
# exact likelihood
result2$par











