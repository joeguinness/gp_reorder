

# I think this is what expand.grid does. So maybe this should be updated
# to simply call that function.
createGrid <- function(nvec){
    
    if(missing(nvec) || length(nvec) == 0) stop("Must supply grid dimensions")
    d <- length(nvec)
    n <- prod(nvec)
    gridlocs <- matrix(0,n,d)
    for(j in 1:d){
        perm <- (1:d)[-j]
        perm <- c(j,perm)
        tempvec <- nvec[-j]
        tempvec <- c(nvec[j],tempvec)
        a1 <- ((1:nvec[j])-1/2)/nvec[j]
        a2 <- rep(a1,n/nvec[j])
        a3 <- array(a2,tempvec)
        a4 <- aperm(a3,perm)
        gridlocs[,j] <- c(a4)
    }
    gridlocs
}


# was hoping to use this as a wrapper for different
# simulation methods
simulateLocations <- function(n,d,method="uniform"){
    if( method == "uniform" )  locs <- matrix(runif(n*d),n,d)   
}

# simulate from a jittered grid.
# jittersize = 0 has no jittering (i.e. regular grid).
# I believe jittersize = 1 is the maximum amount of jittering
# that avoids points overlapping each other
simulateGrid <- function(nvec,jittersize=0){
    n <- prod(nvec)
    d <- length(nvec)
    gridlocs <- createGrid(nvec)
    u <- matrix(runif(n*d)-1/2,n,d) %*% diag(jittersize/nvec)
    locs <- gridlocs + u
}
