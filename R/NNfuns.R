


# functions to find nearest neighbors and do the grouping operations


# naive nearest neighbor finder

findOrderedNN <- function( locs, m ){
    # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
    # by convention, this includes locs[j,], which is distance 0
    n <- dim(locs)[1]
    NNarray <- matrix(NA,n,m+1)
    for(j in 1:n ){
        distvec <- c(rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
        NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
    }    
    NNarray
}


# include some distance neighbors, with proportion distant 1-pnear

findOrderedNNdistant <- function( locs, m, pnear = 1 ){
    # find the pnear(m+1) nearest neighbors to locs[j,] in locs[1:j,]
    # by convention, this includes locs[j,], which is distance 0

    # also adds (1-pnear)(m+1) points from the distant past
    # for a total of m+1 neighbors
    
    n <- dim(locs)[1]
    NNarray <- matrix(NA,n,m+1)
    mnear <- floor(pnear*(m+1))
    mfar <- m+1 - mnear
    
    for(j in 1:n ){
        distvec <- c(rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
        odist <- order(distvec)
        if( j <= m+1) NNarray[j,1:j] <- odist[1:j]        
        if( j > m+1){
            NNarray[j,1:mnear] <- odist[1:mnear]
            if( mfar > 0 ){
                # have to double check to make sure this has no duplicates
                # (and that it's doing what I think it is
                leftinds <- (mnear+1):j
                inds <- floor( seq(mnear+1,j,length.out=mfar) ) 
                NNarray[j,(mnear+1):(m+1)] <- odist[inds]
            }
        }
    }    
    NNarray
}



# take in an array of nearest neighbors, and automatically group
# the observations into groups that share neighbors
# this is helpful to speed the computations and improve their accuracy

groupNN <- function(NNarray){
    n <- nrow(NNarray)
    m <- ncol(NNarray)-1
    
    clust <- vector("list",n)
    for(j in 1:n) clust[[j]] <- j
    for( ell in 2:(m+1) ){  # 2:(m+1)?
        sv <- which( NNarray[,1] - NNarray[,ell] < n )
        for(j in sv){
            k <- NNarray[j,ell]
            if( length(clust[[k]]) > 0){
                nunique <- length(unique(c(NNarray[c(clust[[j]],clust[[k]]),])))

                # this is the rule for deciding whether two groups
                # should be combined
                if( nunique^2 <= length(unique(c(NNarray[clust[[j]],])))^2 + length(unique(c(NNarray[clust[[k]],])))^2 ) {
                    clust[[j]] <- c(clust[[j]],clust[[k]])
                    clust[[k]] <- numeric(0)
                }
            }
        }
    }
    zeroinds <- unlist(lapply(clust,length)==0)
    clust[zeroinds] <- NULL
    NNlist <- lapply(clust,function(inds) NNarray[inds,,drop=FALSE])
    return(NNlist)
}



# faster algorithm to find nearest neighbors. This one splits the
# observation domain into grid boxes and searches neighboring
# grid boxes to find neighbors

findOrderedNNfast <- function( locs, m ){

    n <- dim(locs)[1]
    
    # number of grid boxes in each dimension
    nside <- ceiling( 1/10*sqrt(n) ) # perhaps change to 2*sqrt(n)

    eps <- sqrt(.Machine$double.eps)
    
    # this is a 3D array that will hold the indices of the points
    # within each grid box. This is a bit tricky because we don't
    # know ahead of time how many points are in each grid box. Not sure
    # of the best way to deal with that problem. For now I just pick
    # something big (10*n/nside^2)
    indcube <- array(0,c(nside,nside,ceiling(10*n/nside^2))) 
    
    # rectangular bounds of observation domain
    lims <- matrix( c(apply(locs,2,min)-eps, apply(locs,2,max)+eps ),2,2 )

    # simply round the coordinates to assign them to grid boxes    
    locround <- ceiling( cbind( nside*(locs[,1]-lims[1,1])/(lims[1,2]-lims[1,1]),
                                nside*(locs[,2]-lims[2,1])/(lims[2,2]-lims[2,1]) ) )

    # could be faster, but this is not the bottle neck.
    # just puts the indices into indcube, the object that 
    # maps grid boxes to observation indices
    for(j in 1:n){
        inds <- locround[j,]
        k <- which( indcube[inds[1],inds[2],] == 0 )[1]
        indcube[inds[1],inds[2],k] <- j
    }

    # initialize    
    NNarray <- matrix(NA,n,m+1)

    # loop over all of the grid boxes
    for(i1 in 1:nside){
        for(i2 in 1:nside){
            
            # number of observations in current grid box
            # maybe better to initialize indcube with NAs
            nk <- sum( indcube[i1,i2,] > 0 )    
            
            # loop over observation indices
            for(k in seq(length.out=nk) ){
                j <- indcube[i1,i2,k]
                
                # s gives the number of neighboring grid boxes to search
                # s=1 means search current and adjacent 8 grid boxes 
                s <- max(1,floor( 1/2*nside^2/j ))
                subind = c()
                nlocal <- 0

                # increase s until we have found "enough" PREVIOUS points
                # in neighboring grid boxes. Must be PREVIOUS in ordering
                while( nlocal < min(j,m*2) ){  # could be m*3?
                    l1 <- max(1,i1-s):min(nside,i1+s)
                    l2 <- max(1,i2-s):min(nside,i2+s)
                    subcube <- indcube[l1,l2,]
                    subind <- subcube[ subcube > 0 & subcube <= j ] # PREVIOUS!
                    nlocal <- length(subind)
                    s <- max(s*2,s+1)
                }
                # pick the closest m of them
                subdistmat <- rdist(locs[j,,drop=FALSE],locs[subind,,drop=FALSE])
                orderind <- order(subdistmat)
                NNarray[j,1:min(j,m+1)] <- as.integer(subind[orderind[1:min(j,m+1)]])
            }
        }
    }
    return(NNarray)
}
