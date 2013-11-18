readPP.3D <- function(dpath,recursive=F){
# 
#
# read plasma parameters from files
# 
# I. Virtanen 2010, 2013
# 

    if(is.null(dpath))   return(NULL)
    if(length(dpath)==0) return(NULL)
    
    fpath <- dpath[!file.info(dpath)$isdir]
    fpath <- fpath[length(grep('PP.Rdata',fpath))>0]
    dpath <- dpath[file.info(dpath)$isdir]
    
    # list the PPI result files
    flist <- dir(dpath[1],pattern='PP.Rdata',recursive=recursive,full.names=TRUE)
    if(length(dpath)>1){
        for(k in seq(2,length(dpath))){
            flist <- c(flist,dir(dpath[1],pattern='PP.Rdata',recursive=recursive,full.names=TRUE))
        }
    }
    flist <- c(flist,fpath)
    
    if(length(flist)==0) return(NULL)

    # read first of the data files
    load(flist[1])

    # number of data files
    nFile  <- length(flist)
    # number of height gates
    nHeight <- length(PP$height)
    # number of plasma parameters
    nPar   <- dim(PP$param)[2]
    mi     <- PPI_param$mi
    # number of receiver sites
    nSites <- dim(PP$sites)[1]

    # allocate the necessary arrays
    param     <- array(NA,dim=c(nHeight,nPar+4*nSites+5,nFile))
    std       <- array(NA,dim=c(nHeight,nPar+4*nSites+5,nFile))
    model     <- array(NA,dim=c(nHeight,nPar+4*nSites+5,nFile))
    height    <- array(NA,dim=c(nHeight,nFile))
    status    <- array(NA,dim=c(nHeight,nFile))
    chisqr    <- array(NA,dim=c(nHeight,nFile))
    sites     <- array(NA,dim=c(dim(PP$sites),nFile))
    time_sec  <- vector(length=nFile,mode='numeric')
    date      <- vector(length=nFile,mode='list')
    POSIXtime <- vector(length=nFile,mode='list')
    llhT      <- PP$llhT
    llhR      <- PP$llhR
    covar     <- array(NA,dim=c(nHeight,nPar+4*nSites+5,nPar+4*nSites+5,nFile))
    
    # read the data from files
    for (k in seq(nFile)){
        load(flist[k])
        # check the number of unknowns, if it is larger than the current nPar we must reallocate
        nSitesk <- dim(PP$sites)[1]
        nPark <- dim(PP$param)[2]
        if(nPark>nPar){
            partmp <- param
            stdtmp <- std
            modeltmp <- model
            covartmp <- covar
            sitestmp <- sites
            
            param     <- array(NA,dim=c(nHeight,nPark+4*nSitesk+5,nFile))
            std       <- array(NA,dim=c(nHeight,nPark+4*nSitesk+5,nFile))
            model     <- array(NA,dim=c(nHeight,nPark+4*nSitesk+5,nFile))
            covar     <- array(NA,dim=c(nHeight,nPark+4*nSitesk+5,nPark+4*nSitesk+5,nFile))
            sites     <- array(NA,dim=c(dim(PP$sites),nFile))
            
            param[,1:nPar,] <- partmp[,1:nPar,]
            param[,(nPark+1):(nPark+5),] <- partmp[,(nPar+1):(nPar+5),]
            param[,(nPark+6):(nPark+5+nSites),] <- partmp[,(nPar+6):(nPar+5+nSites),]
            param[,(nPark+nSitesk+6):(nPark+nSitesk+5+nSites),] <- partmp[,(nPar+nSites+6):(nPar+5+2*nSites),]
            param[,(nPark+2*nSitesk+6):(nPark+2*nSitesk+5+nSites),] <- partmp[,(nPar+2*nSites+6):(nPar+5+3*nSites),]
            param[,(nPark+3*nSitesk+6):(nPark+3*nSitesk+5+nSites),] <- partmp[,(nPar+3*nSites+6):(nPar+5+4*nSites),]
            
            std[,1:nPar,] <- stdtmp[,1:nPar,]
            std[,(nPark+1):(nPark+5),] <- stdtmp[,(nPar+1):(nPar+5),]
            std[,(nPark+6):(nPark+5+nSites),] <- stdtmp[,(nPar+6):(nPar+5+nSites),]
            std[,(nPark+nSitesk+6):(nPark+nSitesk+5+nSites),] <- stdtmp[,(nPar+nSites+6):(nPar+5+2*nSites),]
            std[,(nPark+2*nSitesk+6):(nPark+2*nSitesk+5+nSites),] <- stdtmp[,(nPar+2*nSites+6):(nPar+5+3*nSites),]
            std[,(nPark+3*nSitesk+6):(nPark+3*nSitesk+5+nSites),] <- stdtmp[,(nPar+3*nSites+6):(nPar+5+4*nSites),]
            
            model[,1:nPar,] <- modeltmp[,1:nPar,]
            model[,(nPark+1):(nPark+5),] <- modeltmp[,(nPar+1):(nPar+5),]
            model[,(nPark+6):(nPark+5+nSites),] <- modeltmp[,(nPar+6):(nPar+5+nSites),]
            model[,(nPark+nSitesk+6):(nPark+nSitesk+5+nSites),] <- modeltmp[,(nPar+nSites+6):(nPar+5+2*nSites),]
            model[,(nPark+2*nSitesk+6):(nPark+2*nSitesk+5+nSites),] <- modeltmp[,(nPar+2*nSites+6):(nPar+5+3*nSites),]
            model[,(nPark+3*nSitesk+6):(nPark+3*nSitesk+5+nSites),] <- modeltmp[,(nPar+3*nSites+6):(nPar+5+4*nSites),]
            
            covar[,1:nPar,1:nPar,] <- covartmp[,1:nPar,1:nPar,]
            covar[,(nPark+1):(nPark+5),(nPark+1):(nPark+5),] <- covartmp[,(nPar+1):(nPar+5),(nPar+1):(nPar+5),]
            covar[,(nPark+6):(nPark+5+nSites),(nPark+6):(nPark+5+nSites),] <- covartmp[,(nPar+6):(nPar+5+nSites),(nPar+6):(nPar+5+nSites),]
            covar[,(nPark+nSitesk+6):(nPark+nSitesk+5+nSites),(nPark+nSitesk+6):(nPark+nSitesk+5+nSites),] <- covartmp[,(nPar+nSites+6):(nPar+5+2*nSites),(nPar+nSites+6):(nPar+5+2*nSites),]
            covar[,(nPark+2*nSitesk+6):(nPark+2*nSitesk+5+nSites),(nPark+2*nSitesk+6):(nPark+2*nSitesk+5+nSites),] <- covartmp[,(nPar+2*nSites+6):(nPar+5+3*nSites),(nPar+2*nSites+6):(nPar+5+3*nSites),]
            covar[,(nPark+3*nSitesk+6):(nPark+3*nSitesk+5+nSites),(nPark+3*nSitesk+6):(nPark+3*nSitesk+5+nSites),] <- covartmp[,(nPar+3*nSites+6):(nPar+5+4*nSites),(nPar+3*nSites+6):(nPar+5+4*nSites),]

            nPar <- nPark
            nSites <- nSitesk
            
            if(!all(sitestmp[,,(k-1)]==unique(PP[["sites"]][1:dim(sitestmp)[1],1:dim(sitestmp)[2]]))) warning("Site indexing is not identical in all integration periods")
            
            rm(partmp,stdtmp,modeltmp,covartmp,sitestmp)
        }
        param[,1:nPark,k]   <- PP$param
        std[,1:nPark,k]     <- PP$std
        model[,1:nPark,k]   <- PP$model
        chisqr[,k]   <- PP$chisqr
        status[,k]   <- PP$status
        time_sec[k]  <- PP$time_sec
        date[[k]]    <- PP$date
        POSIXtime    <- PP$POSIXtime
        height[,k]   <- PP$height
        sites[1:nSitesk,,k]   <- PP$sites

        for(r in seq(nHeight)){
            # if there will be only NA's there is no need to proceed
            if(any(!is.na(PP$covar[[r]]))){
                covar[r,1:nPark,1:nPark,k] <- PP$covar[[r]]
                # ion temperature (Ti_par + 2*Ti_perp)/3
                param[r,nPar+1,k] <- (PP$param[r,2] + 2*PP$param[r,3] ) / 3 
                # electron temperature (Te_par + 2*Te_perp)/3
                param[r,nPar+2,k] <- (PP$param[r,4] + 2*PP$param[r,5] ) / 3
                # ion velocity along magnetic field
                param[r,nPar+3,k] <- PP$param[r,7:9]%*%PP$B[r,]/sqrt(sum(PP$B[r,]**2))
                # ion perpendicular/parallel temperature ratio
                # approximation of the ratio of two correlated normal random variables
                # from Hayya et al., A note on the ratio of two normally distributed variables,
                # Management Science, 21 (11), 1338-1341, 1975
                param[r,nPar+4,k] <- PP$param[r,3]/PP$param[r,2] +
                    PP$covar[[r]][2,2]*PP$param[r,3]/PP$param[r,2]**3 -
                        PP$covar[[r]][2,3]/PP$param[r,2]**2
                # electron perpendicular/parallel temperature ratio
                param[r,nPar+5,k] <- PP$param[r,5]/PP$param[r,4] +
                    PP$covar[[r]][4,4]*PP$param[r,5]/PP$param[r,4]**3 -
                        PP$covar[[r]][4,5]/PP$param[r,4]**2
                std[r,nPar+1,k] <- sqrt(c(1,2)%*%PP$covar[[r]][2:3,2:3]%*%c(1,2)/sum(c(1,2)**2))
                std[r,nPar+2,k] <- sqrt(c(1,2)%*%PP$covar[[r]][4:5,4:5]%*%c(1,2)/sum(c(1,2)**2))
                std[r,nPar+3,k] <- sqrt(PP$B[r,]%*%PP$covar[[r]][7:9,7:9]%*%PP$B[r,]/sum(PP$B[r,]**2))
                std[r,nPar+4,k] <- sqrt( PP$covar[[r]][2,2]*PP$param[r,3]**2/PP$param[r,2]**4 +
                                        PP$covar[[r]][3,3]/PP$param[r,2]**2 -
                                        2*PP$covar[[r]][2,3]*PP$param[r,3]/PP$param[r,2]**3
                                        )
                std[r,nPar+5,k] <- sqrt( PP$covar[[r]][4,4]*PP$param[r,5]**2/PP$param[r,4]**4 +
                                        PP$covar[[r]][5,5]/PP$param[r,4]**2 -
                                        2*PP$covar[[r]][4,5]*PP$param[r,5]/PP$param[r,4]**3
                                        )
                for( s in seq(nSitesk)){
                    
                    # ion velocity seen at site s
                    param[r,nPar+2*s+4,k] <- PP$param[r,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sqrt(sum(PP$intersect[[r]][[s]]$k.ENU**2))
                    std[r,nPar+2*s+4,k] <- sqrt(PP$intersect[[r]][[s]]$k.ENU%*%PP$covar[[r]][7:9,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sum(PP$intersect[[r]][[s]]$k.ENU**2))
                    # horizontal ion velocity component seen at site s
                    khor <- c(PP$intersect[[r]][[s]]$k.ENU[c(1,2)],0)
                    khor <- khor / sqrt(sum(khor**2))
                    param[r,nPar+2*s+5,k] <- PP$param[r,7:9]%*%khor
                    std[r,nPar+2*s+5,k] <- sqrt(khor%*%PP$covar[[r]][7:9,7:9]%*%khor)
                    # ion temperature at site s
                    phi <- radarPointings:::vectorAngle.cartesian(PP$intersect[[r]][[s]]$k.ENU,PP$B[r,],degrees=FALSE)
                    prvec <- c( cos(phi)**2 , sin(phi)**2 )
                    param[r,nPar+2*nSites+s+5,k] <- PP$param[r,2:3]%*%prvec
                    std[r,nPar+2*nSites+s+5,k] <- sqrt(prvec%*%PP$covar[[r]][2:3,2:3]%*%prvec)
                    # electron temperature at site s
                    param[r,nPar+3*nSites+s+5,k] <- PP$param[r,4:5]%*%prvec
                    std[r,nPar+3*nSites+s+5,k] <- sqrt(prvec%*%PP$covar[[r]][4:5,4:5]%*%prvec)
                }
                
                # remove parameters that must be based solely on the prior model
                if(!is.null(PP$contribSites)){
                    # if there is only one site, remove all but the projections to that receiver
                    if(length(PP$contribSites[[r]])==1){
                        param[r,c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5),k] <- NA
                        std[r,c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5),k] <- NA
                        for(s in setdiff(PP$contribSites[[r]],seq(nSitesk))){
                            param[r,c(nPar+2*s+4,nPar+2*s+5,nPar+2*nSites+s+5,nPar+3*nSites+s+5),k] <- NA
                            std[r,c(nPar+2*s+4,nPar+2*s+5,nPar+2*nSites+s+5,nPar+3*nSites+s+5),k] <- NA
                        }
                    }
                    # if there are two sites we can accept also the temperature anisotropy estimates
                    if(length(PP$contribSites[[r]])==1){
                        param[r,c(7,8,9,nPar+3),k] <- NA
                        std[r,c(7,8,9,nPar+3),k] <- NA
                        for(s in setdiff(PP$contribSites[[r]],seq(nSitesk))){
                            param[r,c(nPar+2*s+4,nPar+2*s+5,nPar+2*nSites+s+5,nPar+3*nSites+s+5),k] <- NA
                            std[r,c(nPar+2*s+4,nPar+2*s+5,nPar+2*nSites+s+5,nPar+3+nSites+s+5),k] <- NA
                        }
                    }
                }
            }
        }
    }
  
    dimnames(param) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
    dimnames(std) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
    dimnames(model) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
    dimnames(covar) <- c(list(paste(seq(nHeight))),lapply(dimnames(PP[["covar"]][[1]]),FUN=function(x,nSites){c(x[1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep=''))},nSites=nSites),list(paste(seq(nFile))))

    return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,height=height,time_sec=time_sec,date=date,POSIXtime=POSIXtime,sites=sites,n=nFile,nPar=nPar,nHeight=nHeight,mi=mi,covar=covar))
    
} # readPP



