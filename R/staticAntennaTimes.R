##
## read antenna pointing directions from LPI output directory and
## return a vector of unix times, such that when used as integratin period
## period limits in plasma parameter estimation, the time periods during
## which both TX and RX antennas were static are integrated together,
## and the time resolution of the LPI analysis is used when either of
## the antennas is moving. 
##
##
## I. Virtanen 2013
##
##


staticAntennaTimes <- function( ddir='.' , recursive=TRUE )
    {
        # list LPI output files
        dfiles <- dir( ddir , pattern="[[:digit:]]LP.Rdata" , full.names=TRUE , recursive=recursive )

        # number of files
        nfiles <- length(dfiles)

        # a vector for timestamps
        times <- rep( NA , length=nfiles)

        # a matrix for pointing directions
        pointings <- matrix(nrow=nfiles,ncol=4)

        for( k in seq( nfiles ) ){
            # load a data file
            load( dfiles[k] )
            # read the timestamp
            times[k] <- ACF[["time.s"]]
            # read TX and RX pointings
            pointings[k,] <- c(ACF[["azelT"]],ACF[["azelR"]])
        }

        # order everything according to the timestamps
        torder <- order( times )
        times <- times[torder]
        pointings <- pointings[ torder , ]
        dfiles <- dfiles[ torder ]

        # differences in pointing directions
        dpoint <- abs( diff( pointings ) )

        # a difference larger than 0.1 degrees is considered as a movement
        movingInds <-  apply( dpoint , function(x){ any(x>.1)} , MARGIN=1 ) 

        # repeat the first and last index
        movingInds <- c( movingInds[1] , movingInds , movingInds[nfiles-1] )

        # select the integration periods whose pointing directions differ from
        # pointings in both adjacent measurements
        iperInds <- which( movingInds[1:nfiles] & movingInds[2:(nfiles+1)] )

        # load the first data file to check its integratio time
        load( dfiles[1] )
        
        times <- unique( c( times[1]-ACF[["timeRes.s"]]  , times[iperInds] , times[nfiles] ) )

        return(times)



    }
