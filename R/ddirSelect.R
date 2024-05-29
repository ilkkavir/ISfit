ddirSelect <- function(mlatlims,mlonlims,hlims,ddirs){

    nbinmlat <- length(mlatlims)-1
    nbinmlon <- length(mlonlims)-1
    nbinlatlon <- nbinmlat*nbinmlon
    hmaxmin <- range(hlims)

    ndirs <- length(ddirs)

    keepdirs <- list()
    for(ibin in seq(nbinlatlon)){
        keepdirs[[ibin]] <- rep(FALSE,ndirs)
    }

    for(idir in seq(ndirs)){
        df <- dir(ddirs[idir],pattern='LP.Rdata',full.names=T)
        load(df[1])
        t <- as.POSIXlt( ACF$time.s , origin='1970-01-01' , tz='utc')
        date <- c(t$year+1900,t$mon+1,t$mday,t$hour,t$min,t$sec)
        nr <- length(ACF$range.km)

        llhr <- matrix(NA,nrow=nr,ncol=3)
        for(ir in seq(nr)){
            llhr[ir,] <- range2llh( r=ACF$range.km[ir]*1000 , llhT=ACF$llhT , azelT=ACF$azelT , llhR=ACF$llhR )
        }
        
        llhm <- aacgmv2(llhr[,1],llhr[,2],llhr[,3]/1000,date,'G2A')

        for(ibin in seq(nbinlatlon)){
            imlat <- ibin%%nbinmlat
            if(imlat==0) imlat <- nbinmlat
            imlon <- floor((ibin-1)/nbinmlat)+1

            keepdirs[[ibin]][idir]  <-  any(llhm[[1]]>=mlatlims[imlat] & llhm[[1]]<=mlatlims[imlat+1]&llhm[[2]]>=mlonlims[imlon]&llhm[[2]]<=mlonlims[imlon+1]&llhr[,3]>=1000*hmaxmin[1]&llhr[,3]<=1000*hmaxmin[2])
        }
        
    }

    return(keepdirs)
#    return(ddirs[keepdirs])

}
