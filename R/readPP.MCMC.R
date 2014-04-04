readPP.MCMC <- function( dpath , recursive=TRUE , MCMClist=FALSE)
    {
        #
        # read plasma parameters from the MCMC fit results
        # This function uses the MCMC part of the
        # files, the iteration results can be read with
        # function readPP.3D
        #
        # I. Virtanen 2014
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
        date      <- vector(length=nFile,mode='list')
        POSIXtime <- vector(length=nFile,mode='list')
        llhT      <- PP$llhT
        llhR      <- PP$llhR
        covar     <- array(NA,dim=c(nHeight,nPar+4*nSites+7,nPar+4*nSites+7,nFile))
        MCMC  <- list()


        # read the data from files
        for (k in seq(nFile)){
            MCMC[[k]] <- vector(mode='list',length=nHeight)
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
            model[,1:nPark,k]   <- PP$model
            time_sec[k]         <- PP$time_sec
            date[[k]]           <- PP$date
            POSIXtime[[k]]      <- PP$POSIXtime
            height[,k]          <- PP$height
            sites[1:nSitesk,,k] <- PP$sites

            for(r in seq(nHeight)){

                if(!is.null(PP[["MCMC"]][[r]])){
                    # allocate a new matrix to the MCMC list, it will be deleted afterwards if MCMClist==FALSE
                    nMCMC <- dim(PP[["MCMC"]][[r]][["pars"]])
                    MCMC[[k]][[r]] <- matrix(nrow=nMCMC[1],ncol=(nPar+4*nSites+7))
                    MCMC[[k]][[r]][,1:nMCMC[2]] <- PP[["MCMC"]][[r]][["pars"]]
                    # ion temperature (Ti_par + 2*Ti_perp)/3
                    MCMC[[k]][[r]][,nPar+1] <- ( MCMC[[k]][[r]][,2] + 2* MCMC[[k]][[r]][,3]) / 3
                    # electron temperature (Te_par + 2*Te_perp)/3
                    MCMC[[k]][[r]][,nPar+2] <- ( MCMC[[k]][[r]][,4] + 2* MCMC[[k]][[r]][,5]) / 3
                    # ion velocity along magnetic field
                    MCMC[[k]][[r]][,nPar+3] <- MCMC[[k]][[r]][,7:9]%*%PP$B[r,]/sqrt(sum(PP$B[r,]**2))
                    # ion perpendicular/parallel temperature ratio
                    MCMC[[k]][[r]][,nPar+4] <- MCMC[[k]][[r]][,3] / MCMC[[k]][[r]][,2]
                    # electron perpendicular/parallel temperature ratio
                    MCMC[[k]][[r]][,nPar+5] <- MCMC[[k]][[r]][,5] / MCMC[[k]][[r]][,4]
                    # ion velocity components perpendicular to the magnetic field..
                    # geomagnetic north
                    Bhor <- c(PP$B[r,1:2],0)
                    # geomagnetic east is perpendicular both to B and Bhor
                    Bx <- radarPointings:::vectorProduct.cartesian(PP$B[r,],Bhor)
                    Bx <- Bx / sqrt(sum(Bx**2))

                    By <- radarPointings:::vectorProduct.cartesian(Bx,PP$B[r,])
                    By <- By / sqrt(sum(By**2))

                    MCMC[[k]][[r]][,nPar+6] <- MCMC[[k]][[r]][,7:9]%*%Bx
                    MCMC[[k]][[r]][,nPar+7] <- MCMC[[k]][[r]][,7:9]%*%By

                    for( s in PP$contribSites[[r]]){

                        # ion velocity seen at site s
                        MCMC[[k]][[r]][,nPar+2*s+6] <- MCMC[[k]][[r]][,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sqrt(sum(PP$intersect[[r]][[s]]$k.ENU**2))
                        # horizontal ion velocity component seen at site s
                        khor <- c(PP$intersect[[r]][[s]]$k.ENU[c(1,2)],0)
                        khor <- khor / sqrt(sum(khor**2))
                        MCMC[[k]][[r]][,nPar+2*s+7] <- MCMC[[k]][[r]][,7:9]%*%khor
                        # copy the horizontal velocity projections to all beams of a multibeam receiver
                        for( ss in seq( nSitesk ) ){
                            if( !any(is.na(PP[["sites"]][ss,2:10])) ){
                                if( sum( PP[["sites"]][ss,2:10] - PP[["sites"]][s,2:10] ) == 0 ){
                                    if(any(is.na(MCMC[[k]][[r]][,nPar+2*ss+7]))){
                                        MCMC[[k]][[r]][,nPar+2*ss+7] <- MCMC[[k]][[r]][,nPar+2*s+7]
                                    }else if(sd(MCMC[[k]][[r]][,nPar+2*ss+7])>sd(MCMC[[k]][[r]][,nPar+2*s+7])){
                                        MCMC[[k]][[r]][,nPar+2*ss+7] <- MCMC[[k]][[r]][,nPar+2*s+7]
                                    }
                                }
                            }
                        }
                        # ion temperature at site s
                        phi <- radarPointings:::vectorAngle.cartesian(PP$intersect[[r]][[s]]$k.ENU,PP$B[r,],degrees=FALSE)
                        prvec <- c( cos(phi)**2 , sin(phi)**2 )
                        MCMC[[k]][[r]][,nPar+2*nSites+s+7] <- MCMC[[k]][[r]][,2:3]%*%prvec
                        # electron temperature at site s
                        MCMC[[k]][[r]][,nPar+3*nSites+s+7] <- MCMC[[k]][[r]][,4:5]%*%prvec
                    }


                    # remove parameters that must be based solely on the prior model
                    # sometimes there are minor antenna movements that are treated as separate
                    # "sites" by the analysis, strip these off before counting the sites
                    nsitesr <- length(PP$contribSites[[r]])
                    if(nsitesr>1) nsitesr <- dim(unique(floor(t(t(PP$sites[PP$contribSites[[r]],])*c(0,100,10,10,1,1,1,10,10,1,1,1)))))[1]
                    # if there is only one site, remove all but the projections to that receiver
                    if(nsitesr==1){
                        MCMC[[k]][[r]][,c(2,3,4,5,7,8,9,nPar+1,nPar+2,nPar+3,nPar+4,nPar+5,nPar+6,nPar+7)] <- NA
                    }
                    # if there are two sites we can accept also the temperature anisotropy estimates
                    if(nsitesr==2){
                        MCMC[[k]][[r]][,c(7,8,9,nPar+3,nPar+6,nPar+7)] <- NA
                    }

                    # then finally extract the parameters and their covariances from the MCMC data
                    param[r,,k] <- colMeans(MCMC[[k]][[r]])
                    std[r,,k] <- apply(MCMC[[k]][[r]],FUN=sd,MARGIN=2)
                    covar[r,,,k] <- cov(MCMC[[k]][[r]])


                    # if MCMClist==TRUE add correct dimnames, otherwise remove the list element
                    if(MCMClist){
                        dimnames(MCMC[[k]][[r]]) <- list(c(),c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')))
                    }else{
                        MCMC[[k]][[r]] <- NULL
                    }

                }
            }

        }

        dimnames(param) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
        dimnames(std) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
        dimnames(model) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]][1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep='')),paste(seq(nFile)))
        dimnames(covar) <- c(list(paste(seq(nHeight))),lapply(dimnames(PP[["covar"]][[1]]),FUN=function(x,nSites){c(x[1:12],paste('Site',seq(nSites),sep=''),'Ti','Te','ViB','Tiratio','Teratio','ViBx','ViBy', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''),paste('TiR',seq(nSites),sep=''),paste('TeR',seq(nSites),sep=''))},nSites=nSites),list(paste(seq(nFile))))

        return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,height=height,time_sec=time_sec,date=date,POSIXtime=POSIXtime,sites=sites,n=nFile,nPar=nPar,nHeight=nHeight,mIon=mIon,covar=covar,MCMClist=MCMC,functionCall=PP[["functionCall"]]))



    }
