ElectricFieldsF <- function(dpath,hmin=200,hmax=400,vipar0=FALSE,recursive=FALSE){
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

    for(k in seq(nf)){
        load(dfiles[k])
        tmp <- EfieldFPP(PP,hmin=hmin,hmax=hmax,vipar0=vipar0)
        E[k,] <- tmp[["E"]]
        Ecov[k,,] <- tmp[["cov"]]
        time[k] <- as.numeric(PP[["POSIXtime"]])
    }

    return(list(E=E,Ecov=Ecov,time=time))

}
