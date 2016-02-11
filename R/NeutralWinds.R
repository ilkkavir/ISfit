NeutralWinds <- function(dpath,timeRes=600,hlimsE=c(80,150),hlimsF=c(200,400),recursive=FALSE){
    #
    # Neutral winds from multistatic velocity measurements. 
    #
    # Electric field is first estimated from F region observations with time resolution of
    # the original ISfit run.
    # 
    # Neutral winds from the E-region gates are then estimated with the same time resolution,
    # and postintegrated to resolution timeRes (in seconds)
    #
    #
    # INPUT:
    #   dpath     data path(s)
    #   timeRes   time resolution for the neutral wind estimation [s]
    #   hlimsE    minimum and maximum height [km] considered as E region 
    #   hlimsF    minimum and maximum height [km] considered as F region
    #   recursive logical, if TRUE, dpath is searched recursively
    #
    # OUTPUT
    #  a list with elements
    #   nWind     a matrix of neutral wind estimates
    #   std       a matrix of neutral wind standard deviations
    #   time      timestamps [s]
    #   height    height gate centres [km]
    #   
    #  I. Virtanen 2016
    #


    # list data files and count them
    if(is.null(dpath)) return(NULL)
    if(length(dpath)==0) return(NULL)
    fpath <- dpath[!file.info(dpath)$isdir]
    fpath <- fpath[length(grep('PP.Rdata',fpath))>0]
    dpath <- dpath[file.info(dpath)$isdir]
    dfiles <- dir(dpath[1],pattern='PP.Rdata',recursive=recursive,full.names=TRUE)
    if(length(dpath)>1){
        for(k in seq(2,length(dpath))){
            dfiles <- c(dfiles,dir(dpath[1],pattern='PP.Rdata',recursive=recursive,full.names=TRUE))
        }
    }
    dfiles <- c(dfiles,fpath)
    if(length(dfiles)==0) return(NULL)
    nf <- length(dfiles)


    # First analyse with the original time resolution, collect results in arrays
    E <- matrix(ncol=2,nrow=nf)
    Ecov <- array(dim=c(nf,2,2))
    time <- rep(NA,nf)
    # load the first data file to count the E-region gates, we will assume that the gates are fixed
    load(dfiles[1])
    Egates <- which((PP[["height"]]>=hlimsE[1])&(PP[["height"]]<=hlimsE[2]))
    nEgates <- length(Egates)
    nWind <- array(dim=c(nf,nEgates,3))
    nWcov <- array(dim=c(nf,nEgates,3,3))
    
    for(k in seq(nf)){

        # load the file
        load(dfiles[k])

        # Electric field from F region
        tmpE <- EfieldFPP(PP,hmin=hlimsF[1],hmax=hlimsF[2])
        E[k,] <- tmpE[["E"]]
        Ecov[k,,] <- tmpE[["cov"]]
        time[k] <- as.numeric(PP[["POSIXtime"]])

        # Neutral wind in each E region gate
        viE <- PP[["param"]][Egates,c('Vix','Viy','Viz')]
        vicovE <- lapply(PP[["covar"]][Egates],FUN=function(x){return(x[c('Vix','Viy','Viz'),c('Vix','Viy','Viz')])})
        hE <- PP[["height"]][Egates]
        latE <- PP[["latitude"]][Egates]
        lonE <- PP[["longitude"]][Egates]
        BE <- PP[["B"]][Egates,]
        tmpN <- Nwind(E[k,],Ecov[k,,],viE,vicovE,hE,BE,PP[["POSIXtime"]],latE,lonE)
        nWind[k,,] <- tmpN[["nWind"]]
        nWcov[k,,,] <- tmpN[["cov"]]
        cat('\r',k,'/',nf)
    }


    # filter nWind and nWcov
    nWindF <- array(dim=c(nf,nEgates,3))
    nWcovF <- array(dim=c(nf,nEgates,3,3))
    for(k in seq(nf)){
        for(hh in seq(nEgates)){
            Q <- matrix(0,ncol=3,nrow=3)
            Mtmp <- c(0,0,0)
            for(n in seq(k)){
                if( (time[k]-time[n]) <= timeRes){
                    Qn <- tryCatch( solve(nWcov[n,hh,,]) , error=function(e){matrix(0,ncol=3,nrow=3)})
                    if((!(any(is.na(Qn))))&(!any(is.na(nWind[n,hh,])))){
                        Q <- Q + Qn
                        Mtmp <- Mtmp + Qn%*%nWind[n,hh,]
                    }
                }
            }
            nWcovF[k,hh,,] <- tryCatch( solve(Q) , error=function(e){matrix(0,ncol=3,nrow=3)})
            nWindF[k,hh,] <- nWcovF[k,hh,,]%*%Mtmp
        }
    }
    

    return(list(nWind=nWindF,cov=nWcovF,E=E,Ecov=Ecov,time=time,nWindHres=nWind,covHres=nWcov,height=hE))


}
