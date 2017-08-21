readPP.3D <- function(dpath,measuredOnly=T,recursive=F,...){
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
    #
    mIon   <- PP[["mIon"]]
    if(is.null(mIon)){
        if(exists('PPI_param')) mIon <- PPI_param$mi
    }
    # number of receiver sites
    nSites <- dim(PP$sites)[1]

    # allocate the necessary arrays
    param     <- array(NA,dim=c(nHeight,nPar+4*nSites+7,nFile))
    std       <- array(NA,dim=c(nHeight,nPar+4*nSites+7,nFile))
    model     <- array(NA,dim=c(nHeight,nPar+4*nSites+7,nFile))
    height    <- array(NA,dim=c(nHeight,nFile))
    status    <- array(NA,dim=c(nHeight,nFile))
    chisqr    <- array(NA,dim=c(nHeight,nFile))
    sites     <- array(NA,dim=c(dim(PP$sites),nFile))
    time_sec  <- vector(length=nFile,mode='numeric')
    timeLimits <- array(NA,dim=c(2,nFile))
    date      <- vector(length=nFile,mode='list')
    POSIXtime <- vector(length=nFile,mode='list')
    llhT      <- PP$llhT
    llhR      <- PP$llhR
    covar     <- array(NA,dim=c(nHeight,nPar+4*nSites+7,nPar+4*nSites+7,nFile))

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

            param     <- array(NA,dim=c(nHeight,nPark+4*nSitesk+7,nFile))
            std       <- array(NA,dim=c(nHeight,nPark+4*nSitesk+7,nFile))
            model     <- array(NA,dim=c(nHeight,nPark+4*nSitesk+7,nFile))
            covar     <- array(NA,dim=c(nHeight,nPark+4*nSitesk+7,nPark+4*nSitesk+7,nFile))
            sites     <- array(NA,dim=c(dim(PP$sites),nFile))

            param[,1:nPar,] <- partmp[,1:nPar,]
            param[,(nPark+1):(nPark+7),] <- partmp[,(nPar+1):(nPar+7),]
            param[,(nPark+8):(nPark+7+nSites),] <- partmp[,(nPar+8):(nPar+7+nSites),]
            param[,(nPark+nSitesk+8):(nPark+nSitesk+7+nSites),] <- partmp[,(nPar+nSites+8):(nPar+7+2*nSites),]
            param[,(nPark+2*nSitesk+8):(nPark+2*nSitesk+7+nSites),] <- partmp[,(nPar+2*nSites+8):(nPar+7+3*nSites),]
            param[,(nPark+3*nSitesk+8):(nPark+3*nSitesk+7+nSites),] <- partmp[,(nPar+3*nSites+8):(nPar+7+4*nSites),]

            std[,1:nPar,] <- stdtmp[,1:nPar,]
            std[,(nPark+1):(nPark+7),] <- stdtmp[,(nPar+1):(nPar+7),]
            std[,(nPark+8):(nPark+7+nSites),] <- stdtmp[,(nPar+8):(nPar+7+nSites),]
            std[,(nPark+nSitesk+8):(nPark+nSitesk+7+nSites),] <- stdtmp[,(nPar+nSites+8):(nPar+7+2*nSites),]
            std[,(nPark+2*nSitesk+8):(nPark+2*nSitesk+7+nSites),] <- stdtmp[,(nPar+2*nSites+8):(nPar+7+3*nSites),]
            std[,(nPark+3*nSitesk+8):(nPark+3*nSitesk+7+nSites),] <- stdtmp[,(nPar+3*nSites+8):(nPar+7+4*nSites),]

            model[,1:nPar,] <- modeltmp[,1:nPar,]
            model[,(nPark+1):(nPark+7),] <- modeltmp[,(nPar+1):(nPar+7),]
            model[,(nPark+8):(nPark+7+nSites),] <- modeltmp[,(nPar+8):(nPar+7+nSites),]
            model[,(nPark+nSitesk+8):(nPark+nSitesk+7+nSites),] <- modeltmp[,(nPar+nSites+8):(nPar+7+2*nSites),]
            model[,(nPark+2*nSitesk+8):(nPark+2*nSitesk+7+nSites),] <- modeltmp[,(nPar+2*nSites+8):(nPar+7+3*nSites),]
            model[,(nPark+3*nSitesk+8):(nPark+3*nSitesk+7+nSites),] <- modeltmp[,(nPar+3*nSites+8):(nPar+7+4*nSites),]

            covar[,1:nPar,1:nPar,] <- covartmp[,1:nPar,1:nPar,]
            covar[,(nPark+1):(nPark+7),(nPark+1):(nPark+7),] <- covartmp[,(nPar+1):(nPar+7),(nPar+1):(nPar+7),]
            covar[,(nPark+8):(nPark+7+nSites),(nPark+8):(nPark+7+nSites),] <- covartmp[,(nPar+8):(nPar+7+nSites),(nPar+8):(nPar+7+nSites),]
            covar[,(nPark+nSitesk+8):(nPark+nSitesk+7+nSites),(nPark+nSitesk+8):(nPark+nSitesk+7+nSites),] <- covartmp[,(nPar+nSites+8):(nPar+7+2*nSites),(nPar+nSites+8):(nPar+7+2*nSites),]
            covar[,(nPark+2*nSitesk+8):(nPark+2*nSitesk+7+nSites),(nPark+2*nSitesk+8):(nPark+2*nSitesk+7+nSites),] <- covartmp[,(nPar+2*nSites+8):(nPar+7+3*nSites),(nPar+2*nSites+8):(nPar+7+3*nSites),]
            covar[,(nPark+3*nSitesk+8):(nPark+3*nSitesk+7+nSites),(nPark+3*nSitesk+8):(nPark+3*nSitesk+7+nSites),] <- covartmp[,(nPar+3*nSites+8):(nPar+7+4*nSites),(nPar+3*nSites+8):(nPar+7+4*nSites),]

            nPar <- nPark
            nSites <- nSitesk

            if(!all(sitestmp[,,(k-1)]==unique(PP[["sites"]][1:dim(sitestmp)[1],1:dim(sitestmp)[2]]),na.rm=T)) warning("Site indexing is not identical in all integration periods")

            rm(partmp,stdtmp,modeltmp,covartmp,sitestmp)
        }
        param[,1:nPark,k]   <- PP$param
        std[,1:nPark,k]     <- PP$std
        model[,1:nPark,k]   <- PP$model
        chisqr[,k]          <- PP$chisqr
        status[,k]          <- PP$status
        time_sec[k]         <- PP$time_sec
        timeLimits[,k]      <- PP$timeLimits.s
        date[[k]]           <- PP$date
        POSIXtime[[k]]      <- PP$POSIXtime
        height[,k]          <- PP$height
        sites[1:nSitesk,,k] <- PP$sites

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
                # ..ignored the additional terms from the above reference...
                # could perhaps use the full Geary-Hinkley transformation but will leave it like this
                # now. The different approximations could be validated by means of MCMC fits
#                param[r,nPar+4,k] <- PP$param[r,3]/PP$param[r,2] #+
##                    PP$covar[[r]][2,2]*PP$param[r,3]/PP$param[r,2]**3 -
##                        PP$covar[[r]][2,3]/PP$param[r,2]**2
                # replaced the ratio with temperature difference, 2016-01-28
                param[r,nPar+4,k] <- PP$param[r,3] - PP$param[r,2]

#                # electron perpendicular/parallel temperature ratio
#                param[r,nPar+5,k] <- PP$param[r,5]/PP$param[r,4] #+
##                    PP$covar[[r]][4,4]*PP$param[r,5]/PP$param[r,4]**3 -
##                        PP$covar[[r]][4,5]/PP$param[r,4]**2

                # electron temperature difference
                param[r,nPar+5,k] <- PP$param[r,5] - PP$param[r,4]


                std[r,nPar+1,k] <- sqrt(c(1,2)%*%PP$covar[[r]][2:3,2:3]%*%c(1,2)/9)
                std[r,nPar+2,k] <- sqrt(c(1,2)%*%PP$covar[[r]][4:5,4:5]%*%c(1,2)/9)
                std[r,nPar+3,k] <- sqrt(PP$B[r,]%*%PP$covar[[r]][7:9,7:9]%*%PP$B[r,]/sum(PP$B[r,]**2))
#                std[r,nPar+4,k] <- sqrt( PP$covar[[r]][2,2]*PP$param[r,3]**2/PP$param[r,2]**4 +
#                                        PP$covar[[r]][3,3]/PP$param[r,2]**2 #-
##                                        2*PP$covar[[r]][2,3]*PP$param[r,3]/PP$param[r,2]**3
#                                        )
                std[r,nPar+4,k] <- sqrt(PP$covar[[r]][2,2]+PP$covar[[r]][3,3]+2*PP$covar[[r]][2,3])#
#                std[r,nPar+5,k] <- sqrt( PP$covar[[r]][4,4]*PP$param[r,5]**2/PP$param[r,4]**4 +
##                                        PP$covar[[r]][5,5]/PP$param[r,4]**2 #-
#                                        2*PP$covar[[r]][4,5]*PP$param[r,5]/PP$param[r,4]**3
#                                        )
                std[r,nPar+5,k] <- sqrt(PP$covar[[r]][4,4]+PP$covar[[r]][5,5]+2*PP$covar[[r]][4,5])
                # ion velocity components perpendicular to the magnetic field..
                # geomagnetic north
                Bhor <- c(PP$B[r,1:2],0)
                # geomagnetic east is perpendicular both to B and Bhor
                Bx <- radarPointings:::vectorProduct.cartesian(PP$B[r,],Bhor)
                Bx <- Bx / sqrt(sum(Bx**2))

                By <- radarPointings:::vectorProduct.cartesian(Bx,PP$B[r,])
                By <- By / sqrt(sum(By**2))

                param[r,nPar+6,k] <- PP$param[r,7:9]%*%Bx
                param[r,nPar+7,k] <- PP$param[r,7:9]%*%By
                std[r,nPar+6,k] <- sqrt(Bx%*%PP$covar[[r]][7:9,7:9]%*%Bx)
                std[r,nPar+7,k] <- sqrt(By%*%PP$covar[[r]][7:9,7:9]%*%By)

                for( s in PP$contribSites[[r]]){

                    # ion velocity seen at site s
                    param[r,nPar+2*s+6,k] <- PP$param[r,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sqrt(sum(PP$intersect[[r]][[s]]$k.ENU**2))
                    std[r,nPar+2*s+6,k] <- sqrt(PP$intersect[[r]][[s]]$k.ENU%*%PP$covar[[r]][7:9,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sum(PP$intersect[[r]][[s]]$k.ENU**2))
                    # horizontal ion velocity component seen at site s
                    khor <- c(PP$intersect[[r]][[s]]$k.ENU[c(1,2)],0)
                    khor <- khor / sqrt(sum(khor**2))
                    param[r,nPar+2*s+7,k] <- PP$param[r,7:9]%*%khor
                    std[r,nPar+2*s+7,k] <- sqrt(khor%*%PP$covar[[r]][7:9,7:9]%*%khor)
                    # copy the horizontal velocity projections to all beams of a multibeam receiver
                    for( ss in seq( nSitesk ) ){
                        if( !any(is.na(PP[["sites"]][ss,2:10])) ){
                            if( sum( PP[["sites"]][ss,2:10] - PP[["sites"]][s,2:10] ) == 0 ){
                                if(is.na(std[r,nPar+2*ss+7,k])){
                                    param[r,nPar+2*ss+7,k] <-  param[r,nPar+2*s+7,k]
                                    std[r,nPar+2*ss+7,k] <-  std[r,nPar+2*s+7,k]
                                }else if(std[r,nPar+2*ss+7,k]>std[r,nPar+2*s+7,k]){
                                    param[r,nPar+2*ss+7,k] <-  param[r,nPar+2*s+7,k]
                                    std[r,nPar+2*ss+7,k] <-  std[r,nPar+2*s+7,k]
                                }
                            }
                        }
                    }
                    # ion temperature at site s
                    phi <- radarPointings:::vectorAngle.cartesian(PP$intersect[[r]][[s]]$k.ENU,PP$B[r,],degrees=FALSE)
                    prvec <- c( cos(phi)**2 , sin(phi)**2 )
                    param[r,nPar+2*nSites+s+7,k] <- PP$param[r,2:3]%*%prvec
                    std[r,nPar+2*nSites+s+7,k] <- sqrt(prvec%*%PP$covar[[r]][2:3,2:3]%*%prvec)
                    # electron temperature at site s
                    param[r,nPar+3*nSites+s+7,k] <- PP$param[r,4:5]%*%prvec
                    std[r,nPar+3*nSites+s+7,k] <- sqrt(prvec%*%PP$covar[[r]][4:5,4:5]%*%prvec)
                }


                # remove parameters that must be based solely on the prior model
                # sometimes there are minor antenna movements that are treated as separate
                # "sites" by the analysis, strip these off before counting the sites
                if(measuredOnly){
                    nsitesr <- length(PP$contribSites[[r]])
                    #if(nsitesr>1) nsitesr <- dim(unique(floor(t(t(PP$sites[PP$contribSites[[r]],])*c(0,100,10,10,1,1,1,10,10,1,1,1)))))[1]
                    if(nsitesr>1) nsitesr <- dim(unique(floor(t(t(PP$sites[PP$contribSites[[r]],])*c(0,0,10,10,0,1,1,10,10,0,0,0)))))[1]
                    # if there is only one site, remove all but the projections to that receiver
                    if(nsitesr==1){
                        param[r,c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5,nPar+6,nPar+7),k] <- NA
                        std[r,c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5,nPar+6,nPar+7),k] <- NA
                        covar[r,c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5,nPar+6,nPar+7),c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5,nPar+6,nPar+7),k] <- NA
                    }
                    # if there are two sites we can accept also the temperature anisotropy estimates
                    if(nsitesr==2){
                        param[r,c(7,8,9,nPar+3,nPar+6,nPar+7),k] <- NA
                        std[r,c(7,8,9,nPar+3,nPar+6,nPar+7),k] <- NA
                        covar[r,c(7,8,9,nPar+3,nPar+6,nPar+7),c(7,8,9,nPar+3,nPar+6,nPar+7),k] <- NA
                    }
                }
            }
        }
    }

    dimnames(param) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tidiff','Tediff','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
    dimnames(std) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tidiff','Tediff','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
    dimnames(model) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tidiff','Tediff','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
    dimnames(covar) <- c(list(paste(seq(nHeight))),lapply(dimnames(PP[["covar"]][[1]]),FUN=function(x,nSites){c(x[1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tidiff','Tediff','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep=''))},nSites=nSites),list(paste(seq(nFile))))

    if(!measuredOnly) warning("Returning also parameters that are based solely on the prior model.")

    return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,height=height,time_sec=time_sec,timeLimits=timeLimits,date=date,POSIXtime=POSIXtime,sites=sites,n=nFile,nPar=nPar,nHeight=nHeight,mIon=mIon,covar=covar))

}



