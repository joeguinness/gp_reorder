

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
    blocklist
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

