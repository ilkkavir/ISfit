postInteg.time <- function( parlist , timeRes.s ){

    # time range of the measurement
    rtime <- range( parlist[["time_sec"]] , na.rm=TRUE )

    # new integration time limits
    tlimsNew <- seq( rtime[1] , rtime[2] , by=timeRes.s )

    # number of integration periods after post integration
    nIperNew <- length(tlimsNew) - 1

    # dimensions of the original param array
    dimOld <- dim( parlist[["param"]] )

    # allocate new arrays
    param     <- array(NA,dim=c(dimOld[1],dimOld[2],nIperNew))
    std       <- array(NA,dim=c(dimOld[1],dimOld[2],nIperNew))
    model     <- array(NA,dim=c(dimOld[1],dimOld[2],nIperNew))
    height    <- array(NA,dim=c(dimOld[1],nIperNew))
    status    <- array(NA,dim=c(dimOld[1],nIperNew))
    time_sec  <- vector(length=nIperNew,mode='numeric')
    POSIXtime <- vector(length=nIperNew,mode='list')

    statinds <- which( parlist[["status"]]!=0 , arr.ind=TRUE )
    parlist[["param"]][statinds[1],,statinds[2]] <- NA
    parlist[["std"]][statinds[1],,statinds[2]] <- Inf

    # the integration
    for( k in seq( nIperNew ) ){
        indk <- which( ( parlist[["time_sec"]] > tlimsNew[k] ) & (parlist[["time_sec"]] <= (tlimsNew[k+1]+.01)) )
        if( length(indk) > 1 ){
            std[,,k] <- sqrt( 1 / apply( 1 / parlist[["std"]][,,indk]**2 , FUN=sum , MARGIN=c(1,2) ) )
            param[,,k] <- apply( parlist[["param"]][,,indk] / parlist[["std"]][,,indk]**2 , FUN=sum , MARGIN=c(1,2) , na.rm=TRUE ) * std[,,k]**2
            model[,,k] <- apply( parlist[["model"]][,,indk] , FUN=mean , MARGIN=c(1,2) )
            height[,k] <- rowMeans( parlist[["height"]][,indk] )
            status[,k] <- pmin( apply( parlist[["status"]][,indk] , FUN=max , MARGIN=1 ) , apply( parlist[["status"]][,indk] , FUN=min , MARGIN=1 ) )
        }
        time_sec[k] <- tlimsNew[k+1]        
        POSIXtime[[k]] <- as.POSIXlt(time_sec[k],origin='1970-01-01',tz='ut')
    }


    oldnames <- dimnames(parlist[["param"]])

    dimnames( param ) <- list( oldnames[[1]] , oldnames[[2]] , seq(nIperNew) )
    dimnames( std ) <- list( oldnames[[1]] , oldnames[[2]] , seq(nIperNew) )
    dimnames( model ) <- list( oldnames[[1]] , oldnames[[2]] , seq(nIperNew) )
    

    return(list(param=param,std=std,model=model,status=status,height=height,time_sec=time_sec,POSIXtime=POSIXtime,n=nIperNew,nHeight=dimOld[1],mIon=parlist[["mIon"]]))


}
