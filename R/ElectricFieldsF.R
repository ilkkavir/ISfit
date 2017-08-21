ElectricFieldsF <- function(dpath,hmin=200,hmax=400,vipar0=FALSE,recursive=FALSE,chisqrLim=10,...){
    #
    # Electric field from multistatic F-region velocity measurements.
    #
    # I. Virtanen 2016
    #

    if(is.null(dpath)) return(NULL)
    if(length(dpath)==0) return(NULL)

    fpath <- dpath[!file.info(dpath)$isdir]
    fpath <- fpath[length(grep('PP.Rdata',fpath))>0]
    dpath <- dpath[file.info(dpath)$isdir]

    # list the PPI result files
    dfiles <- dir(dpath[1],pattern='PP.Rdata',recursive=recursive,full.names=TRUE)
    if(length(dpath)>1){
        for(k in seq(2,length(dpath))){
            dfiles <- c(dfiles,dir(dpath[1],pattern='PP.Rdata',recursive=recursive,full.names=TRUE))
        }
    }
    dfiles <- c(dfiles,fpath)

    if(length(dfiles)==0) return(NULL)


    nf <- length(dfiles)

    E <- matrix(ncol=2,nrow=nf)
    Ecov <- array(dim=c(nf,2,2))
    time <- rep(NA,nf)
    tres <- rep(NA,nf)
    latitudes <- list()
    longitudes <- list()
    heights <- list()
    B <- list()
    date <- list()
    sites <- list()
    POSIXtime <- list()

    for(k in seq(nf)){
        load(dfiles[k])

        # cut off all failed iterations and large residuals
        rminds <- PP[["status"]]!=0 | PP[["chisqr"]]>chisqrLim

        if (any(rminds)){
            for (irm in which(rminds)){
                PP[["covar"]][[irm]][,] <- Inf
            }
        }

        #PP[["param"]][PP[["status"]]!=0,] <- NA

        # cut off very large chi-squared values
        #PP[["param"]][PP[["chisqr"]]>chisqrLim,] <- NA

        tmp <- EfieldFPP(PP,hmin=hmin,hmax=hmax,vipar0=vipar0,...)
        E[k,] <- tmp[["E"]]
        Ecov[k,,] <- tmp[["cov"]]
        time[k] <- tmp[["time"]]
        tres[k] <- tmp[["tres"]]
        latitudes[[k]] <- tmp[["latitudes"]]
        longitudes[[k]] <- tmp[["longitudes"]]
        heights[[k]] <- tmp[["heights"]]
        B[[k]] <- tmp[["B"]]
        date[[k]] <- tmp[["date"]]
        sites[[k]] <- tmp[["sites"]]
        POSIXtime[[k]] <- tmp[["POSIXtime"]]
    }

    latmin <- sapply(latitudes,FUN=min)
    latmax <- sapply(latitudes,FUN=max)
    lonmin <- sapply(longitudes,FUN=min)
    lonmax <- sapply(longitudes,FUN=max)
    hmin <- sapply(heights,FUN=min)
    hmax <- sapply(heights,FUN=max)



    return(list(E=E,Ecov=Ecov,time=time,tres=tres,
                latitudes=latitudes , longitudes=longitudes ,
                heights=heights, B=B , date=date , sites=sites ,
                latmin=latmin , latmax=latmax ,
                lonmin=lonmin , lonmax=lonmax ,
                hmin=hmin , hmax=hmax , POSIXtime=POSIXtime)
           )

}
