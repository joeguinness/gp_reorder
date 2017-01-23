

# All of these functions assume that the arguments y and locs, which
# hold the data and their locations, have already been ordered, so we
# simply operate on them in the canonical ordering (i.e. there is no
# need to reorder them again after they have been ordered and input
# into these functions)


# ordered composite likelihood. This is Vecchia's approximation
# NNarray is an array whose first column is 1 through n, representing
# indices for the observations. The other columns give the indices of 
# the observations we condition on. 
# Thus, NNarray[i,j] = i for j = 1, and NNarray[i,j] < i for j > 1.
# this particular function just does apply() on the function mvnCondLikInds
# and sums the loglikelihoods.
orderedCompLik <- function( covparms, covfun, y, locs, NNarray, returnD1 = FALSE ){
    ocl <- apply(NNarray,1,mvnCondLikInds,covparms=covparms,covfun=covfun,y=y,locs=locs,lastk=1,returnD1=returnD1)
    if(returnD1) return(rowSums(ocl))
    else return(sum(ocl))
}


# We can probably use mapply() or Map() to absorb this intermediate
# step into the orderedCompLik function. But here it is
mvnCondLikInds <- function( inds, covparms, covfun, y, locs, lastk = 1, returnD1 ){
    # inds in NNarray are ordered nearest current point to furthest
    # so if we want current point last, need to reverse
    # there are also NAs if current point early in the ordering
    inds0 <- rev(inds[!is.na(inds)])
    mvnCondLik(covparms,covfun,y[inds0],locs[inds0,,drop=FALSE],lastk,returnD1)
}


# Same thing for the grouped version of the likelihood.
# NNlist is a list whose elements are subarrays of NNarray, representating
# the observations that get grouped together.
orderedGroupCompLik <- function( covparms, covfun, y, locs, NNlist ){

    ocl <- lapply(NNlist,mvnMultiCondLikInds,covparms=covparms,covfun=covfun,y=y,locs=locs)
    return(sum(unlist(ocl)))
}


# again, this could probably be absorbed with smart use of mapply() of Map()
mvnMultiCondLikInds <- function( indarray, covparms, covfun, y, locs, returnz = FALSE ){
    # inds in NNarray are ordered nearest current point to furthest
    # so if we want current point last, need to reverse
    # there are also NAs if current point early in the ordering
    
    inds1 <- indarray[,1]
    allinds <- sort(unique(c(indarray)), na.last = NA)
    likinds <- which( allinds %in% inds1 )
    mvnMultiCondLik(covparms,covfun,y[allinds],locs[allinds,,drop=FALSE],likinds,returnz)
}


# this is the conditional likelihood function used for 
# the grouped version. We compute the conditional loglikelihood
# for the elements of y specified by likinds, and condition on 
# all previous points y's ordering
mvnMultiCondLik <- function( covparms, covfun, y, locs, likinds, returnz = FALSE ){
    
    # computes the ordered conditional multivariate normal density for the
    # data in the likinds position of y

    # for example, if y has length 4, and likinds = c(2,4), 
    # this returns the log density for the second observation 
    # given the first plus the log density for the fourth 
    # observation given the first three
    
    leny    <- length(y)
    distmat <- rdist(locs,locs)
    covmat  <- covfun( locs, covparms, returnD1 = FALSE )
    cholmat <- tryCatch(t(chol(covmat)) , error = function(a) numeric(0) )
    if( length(cholmat) == 0 ){
        # sometimes the covmat is not numerically positive definite.
        # we probably need a better solution for this.
        return(-9999999)
        print("One of the Choleskys failed")
    }

    # decorrelating transform
    z       <- forwardsolve(cholmat,y)
    # sum of log conditional variances
    logsig  <- 2*sum(log(diag(cholmat)[likinds]))
    # simply add them up to get conditional loglikelihood
    loglik  <- -length(likinds)/2*log(2*pi) - 1/2*logsig - 1/2*sum(z[likinds]^2)
}


# this is the non-grouped version. I wrote this to compute the 
# conditional loglikelihood for the 'lastk' observations in y
# given the previous ones. This allows us to use this function
# to get the exact likelihood, by letting y and locs be the whole
# dataset rather than subsets and setting lastk = length(n). It also
# allows us to use this function to compute block independent 
# loglikelihoods by inputting y and locs as subsets (blocks) of the
# whole dataset and setting lastk = length(y)
mvnCondLik <- function( covparms, covfun, y, locs, lastk = 1, returnD1 = FALSE ){
    
    # computes the conditional multivariate normal density for the
    # 'lastk' elements of y given the first
    leny    <- length(y)
    covmat  <- covfun( locs, covparms, returnD1 )
    cholmat <- tryCatch(t(chol(covmat)) , error = function(a) numeric(0) )
    if( length(cholmat) == 0 ) return(-9999999)
    
    z       <- forwardsolve(cholmat,y)
    inds2    <- (leny-lastk+1):leny
    logsig <- 2*sum(log(diag(cholmat)[inds2]))
    loglik  <- -lastk/2*log(2*pi) - 1/2*logsig - 1/2*sum(z[inds2]^2)

    dlogf <- c()
    d2 <- c()
    
    # derivatives
    if( returnD1 ){
        
        dlogf <- rep(0,length(covparms))
        inds1 <- seq(length.out=leny-lastk)
        
        if(length(inds1)==0){
            for(j in 1:length(covparms)){
                dcov <- attr(covmat,"dmats")[[j]]
                dsig <- dcov[inds2,inds2] 
                dlogf[j] <- -1/2*dsig*exp(-logsig) + 1/2*exp(-1*logsig)*dsig*sum(z[inds2]^2)
                d2[j] <- dsig
            }
        } else {
            v1 <- backsolve( t(cholmat[inds1,inds1]), cholmat[inds2,inds1] )
            v2 <- backsolve( t(cholmat[inds1,inds1]), z[inds1] )
            for(j in 1:length(covparms)){
                dcov <- attr(covmat,"dmats")[[j]]
                dsig <- dcov[inds2,inds2] + 
                    crossprod(v1,dcov[inds1,inds1] %*% v1) -
                    2*crossprod(v1,dcov[inds1,inds2])
                x1 <- -1/2*dsig*exp(-logsig)
                x2 <- 1/2*exp(-logsig)*dsig*sum(z[inds2]^2)
                x3 <- z[inds2]*exp(-1/2*logsig)*( crossprod(dcov[inds1,inds2],v2) - 
                                                      crossprod(v1, dcov[inds1,inds1] %*% v2 ) )
                dlogf[j] <- x1 + x2 + x3
                d2[j] <- dsig
            }
        }
    }
    return(c(loglik,dlogf))
}


# using mvnCondLik to get exact loglikelihood
mvnMargLik <- function(covparms,covfun,y,locs){
    lastk = length(y)
    mvnCondLik(covparms,covfun,y,locs,lastk)
}

# auxiliary function that takes in the whole dataset
# as well as a set of indices. Function returns the 
# exact loglikelihood for the data with indices 'inds'
# used with lapply in the next function
mvnMargLikInds <- function(inds,covparms,covfun,y,locs){
    mvnMargLik(covparms,covfun,y[inds],locs[inds,])
}

# blocklist is a list whose elements are vectors of indices
# that define the blocks
mvnIndepBlocks <- function(covparms,covfun,y,locs,blocklist){
    
    # blocklist is a list. each list element is a vector 
    # of indices, defining the blocks.
    liklist <- lapply(blocklist,mvnMargLikInds,
                      covparms=covparms,covfun=covfun,y=y,locs=locs)
    return(sum(unlist(liklist)))
}


# function n1*n2 recangular blocks and returning
# a list of indices for locations that fall inside
# each of the rectangular blocks. 
# Used as input to mvnIndepBlocks
getBlockList <- function(locs,n1,n2){
    
    # blocks the rectangle enclosing locs into a 
    # n1 x n2 grid, and returns a list of the indices
    # of locs that fall in each block
    n <- nrow(locs)
    r1 <- range(locs[,1])
    r2 <- range(locs[,2])
    locsround <- matrix(0,n,2)
    locsround[,1] <- floor( (locs[,1] - r1[1])/diff(r1)*(n1-1e-6) + 1)
    locsround[,2] <- floor( (locs[,2] - r2[1])/diff(r2)*(n2-1e-6) + 1)
    blocklist <- vector("list",length=n1*n2)
    for(j1 in 1:n1){
        for(j2 in 1:n2){
            blocklist[[(j1-1)*n2 + j2]] <- which(locsround[,1]==j1 & locsround[,2]==j2)     
        }
    }
    npoints <- sapply( blocklist,length )
    blocklist[npoints==0] <- NULL
    return(blocklist)
}


# for spde comparisons. This was for the original paper. We can ignore this for now.
spdefun <- function( locs,range0,tau0,cutoff=1e-6,offset=c(0.2,0.2),max.edge=c(0.2,0.2), getKL=FALSE, covmat=NULL, cholmat=NULL ){

    # construct spde mesh
    mesh = inla.mesh.2d(
        loc=locs,
        cutoff = cutoff,
        offset=offset,
        max.edge=max.edge)
    
    kappa0 = sqrt(8)/range0

    ## Construct spde object
    spde=inla.spde2.matern(mesh, alpha=2,
                           B.tau=cbind(log(tau0),1,0),
                           B.kappa=cbind(log(kappa0),0,1),
                           theta.prior.mean=c(0,0),
                           theta.prior.prec=c(0.1,1))

    # get precision matrix Q 
    Q=inla.spde2.precision(
        spde,theta=c(0,0))

    # compute KLdivergence if indicated
    if(getKL){
        allind <- 1:dim(Q)[1]
        idx1 <- mesh$idx$loc   # observation indices
        idx2 <- allind[-idx1]
        Q11 <- Q[idx1,idx1]
        Q12 <- Q[idx1,idx2]
        Q21 <- Q[idx2,idx1]
        Q22 <- Q[idx2,idx2]
        siginv <- Q11 - Q12 %*% solve(Q22,Q21)
        kldiv <- -sum(log(diag(cholmat)))-sum(log(diag(chol(siginv)))) + 1/2*sum(diag(siginv %*% covmat)) - n/2 
        return(kldiv)
    } else { # otherwise just return the precision
        return(Q)
    }
}                 



# these functions are used to compute L^{-1}X, where
# L^{-1} is the approximation to the inverse Cholesky
# implied by the Vecchia approximation, and X
# is an n x p matrix
sparseInvCholSolveRows <- function( indarray, covparms, covfun, X, locs ){

    inds1 <- indarray[,1]
    allinds <- sort(unique(c(indarray)), na.last = NA)
    likinds <- which( allinds %in% inds1 )
    
    locs0 <- locs[allinds,,drop=FALSE]    
    
    X0 <- X[allinds,,drop=FALSE]
    #distmat <- rdist(locs0,locs0)
    covmat  <- covfun( locs0, covparms, returnD1 = FALSE )
    cholmat <- tryCatch(t(chol(covmat)) , error = function(a) numeric(0) )
    if( length(cholmat) == 0 ){
        # sometimes the covmat is not numerically positive definite.
        # we probably need a better solution for this.
        print("One of the Choleskys failed")
        return(cbind(allinds[likinds],X0[likinds,,drop=FALSE]))    
    }
    
    # decorrelating transform
    Z0 <- forwardsolve(cholmat,X0)
    return(cbind(allinds[likinds],Z0[likinds,,drop=FALSE]))    

}


# this is the wrapper that calls lapply and cleans up
sparseInvCholSolve <- function( NNlist, covparms, covfun, X, locs ){
    
    Zlist <- lapply(NNlist, sparseInvCholSolveRows, covparms,covfun,X,locs)
    dimX <- dim(X)
    Z <- array(NA, dimX )
    for(j in 1:length(NNlist)){
        #print(Zlist[[j]])
        Z[ Zlist[[j]][,1], ] <- Zlist[[j]][,2:(dimX[2]+1)]
    }
    return( Z )

}





orderedGroupCompProfLik <- function( covparms, covfun, y, X, locs, NNlist ){
    
    # profile out regression coefficients
    # profile out overall variance
    dimX <- dim(X)
    n <- dimX[1]
    p <- dimX[2]
    yandX <- cbind(y,X)
    covparms[1] <- 1
    Z <- sparseInvCholSolve( NNlist, covparms, covfun, yandX, locs )
    Zcross <- crossprod( Z[,2:(p+1)] )
    if( min( eigen(Zcross)$values ) < 1e-6 ||
        max( eigen(Zcross)$values ) > 1e20 ) return(list(loglik=-999999))
    #print( crossprod( Z[,2:(p+1)] ) )
    betacovmat <- solve( Zcross )
    betahat <- betacovmat %*% crossprod( Z[,2:(p+1)],Z[,1] ) 
    y0 <- y - X %*% betahat
    Z <- sparseInvCholSolve(NNlist,covparms,covfun,y0,locs)
    sigmasq <- c(crossprod( Z )/n)
    covparms[1] <- sigmasq
    print(covparms)
    loglik <- orderedGroupCompLik( covparms, covfun, y0, locs, NNlist )
    return(list(loglik=loglik, sigmasq = sigmasq, 
                betahat = betahat, betacovmat = sigmasq*betacovmat))
}






# wrapper function for model fitting 
# uses grouped profile likelihood, maxmindist ordering implemented,
# right now, cannot fix variance parameter, always profiled out
fitmodel <- function(y, X, locs, covfun, numneighbors = 30, orderfun = "maxmindist", 
                     fixedparameters = NA ){
    
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
    Xord <- X[ord,]
    
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

    # independent blocks negative loglikelihood with Least Squares betahat
    betahat <- solve( crossprod(Xord), crossprod(Xord,yord) )
    f0 <- function(x){
        for(j in 1:length(notfixedinds)){
            covparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](x[j])
        }
        y0 <- y - X %*% betahat
        covparms[1] <- crossprod( y0 )/n
        ll <- mvnIndepBlocks(covparms,covfun,y0,locs,blocklist)
        return(-ll)
    }
    
    # Vecchia's approximation
    f1 <- function(x){
        for(j in 1:length(notfixedinds)){
            covparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](x[j])
        }
        ll <- orderedGroupCompProfLik(covparms,covfun,yord,Xord,locsord,NNlist)
        return(-ll$loglik)
    }
    
    # use nelder mead (which seems more stable) to move towards maximum
    # of independent blocks likelihood approximation
    result0 <- optim(startvals[notfixedinds],f0,method="Nelder-Mead",control=list(maxit=50,trace=1))
    # use BFGS to get to the optimum once we are close
    result0 <- optim(result0$par,f0,method="BFGS",control=list(maxit=100,trace=1))
    # use the maximum independent blocks likelihood estimates to start 
    # an optimization of Vecchia's likelihood
    result1 <- optim(result0$par,f1,method="BFGS",control=list(maxit=20,trace=5))
    outparms <- fixedparameters
    # transform back to original parameter domain with inverse link
    for(j in 1:length(notfixedinds)){
        outparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](result1$par[j])
    }
    proflik <- orderedGroupCompProfLik(outparms,covfun,yord,Xord,locsord,NNlist)
    outparms[1] <- proflik$sigmasq
    names(outparms) <- c("variance","range","smoothness","signal2noise")
    
    # get variance and mean parameters
    returnobj <- list( covparms = outparms, betahat = as.vector(proflik$betahat), 
                       betacovmat = proflik$betacovmat, loglik = proflik$loglik )
    return(returnobj)
}       








