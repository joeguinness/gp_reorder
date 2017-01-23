

# analyze microscale soil elemental data

# clean up workspace and load in functions
rm(list=ls())
library(fields)
library(matrixStats)
sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source)

# load in the data (you'll have to change the path)
load( file = "/Users/jsguinne/Documents/misc/BeforeAndAfterData_AsSandProject_EAlves.RData")

# flag the repeated measurements
badinds <- which( abs( diff(as.matrix(y_pos_After)) ) > 10 ,arr.ind= TRUE)
det2_As_K_After[41:80,1:40] <- NA

# define the element names
elnames <- c("As","Ca","Cr","Cu","Fe","Mn","Ni","Si","Ti","Zn")

# these lines tell us how we are going to subset the data
# right now, it says we take every fourth pixel, if we decrease
# to "by=2" it uses more data (every other pixel) but runs slower
dim_elements <- dim(det2_As_K_After)
inds1 <- seq(1,dim_elements[1],by=4)
inds2 <- seq(1,dim_elements[2],by=4)

# create a matrix to hold the data and fill it in
concentration_all <- matrix(NA, length(inds1)*length(inds2), length(elnames) )
colnames(concentration_all) <- elnames
for(j in 1:length(elnames)){
    objname <- paste0( "det2_",elnames[j],"_K_After" )
    a <- as.matrix(get(objname))
    concentration_all[,j] <- as.vector( a[inds1,inds2] )
}

# get the data locations
locs_all <- cbind(as.vector(as.matrix(x_pos_After[inds1,inds2])), 
                  as.vector(as.matrix(y_pos_After[inds1,inds2])))

# figure out which locations have missing data
ismissing <- which( is.na(concentration_all) | is.nan(concentration_all), arr.ind = TRUE )[,1]

# remove those from the dataset
concentration <- concentration_all[-ismissing,]
locs <- locs_all[-ismissing,]

# response is log arsenic
y <- log( concentration[,1] )

# predictors are log other elements
X <- cbind( rep(1,length(y)), log( concentration[,2:length(elnames)] ) )

# choose covariance function
covfun <- maternIsotropic

# this is the function that fits the model
result <- fitmodel(y, X, locs, covfun, numneighbors = 30, fixedparameters = c(NA,NA,NA,NA) )
# and this gives the Z-scores
zscores <- result$betahat/sqrt(diag(result$betacovmat))


# you can look at the z-scores and decide which columns of X
# you want to drop, and then fit it again with fewer columns

