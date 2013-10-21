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
  param     <- array(NA,dim=c(nHeight,nPar+2*nSites+3,nFile))
  std       <- array(NA,dim=c(nHeight,nPar+2*nSites+3,nFile))
  model     <- array(NA,dim=c(nHeight,nPar+2*nSites+3,nFile))
  height    <- array(NA,dim=c(nHeight,nFile))
  status    <- array(NA,dim=c(nHeight,nFile))
  chisqr    <- array(NA,dim=c(nHeight,nFile))
  sites     <- array(NA,dim=c(dim(PP$sites),nFile))
  time_sec  <- vector(length=nFile,mode='numeric')
  date      <- vector(length=nFile,mode='list')
  POSIXtime <- vector(length=nFile,mode='list')
  llhT      <- PP$llhT
  llhR      <- PP$llhR
  covar     <- array(NA,dim=c(nHeight,nPar+2*nSites+3,nPar+2*nSites+3,nFile))

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

        param     <- array(NA,dim=c(nHeight,nPark+2*nSitesk+3,nFile))
        std       <- array(NA,dim=c(nHeight,nPark+2*nSitesk+3,nFile))
        model     <- array(NA,dim=c(nHeight,nPark+2*nSitesk+3,nFile))
        covar     <- array(NA,dim=c(nHeight,nPark+2*nSitesk+3,nPark+2*nSites+3,nFile))
        sites     <- array(NA,dim=c(dim(PP$sites),nFile))

        param[,1:(nPar+2*nSites+3),] <- partmp
        std[,1:(nPar+2*nSites+3),] <- stdtmp
        model[,1:(nPar+2*nSites+3),] <- modeltmp
        covar[,1:(nPar+2*nSites+3),1:(nPar+nSites+3),] <- covartmp
        sites[1:dim(sitestmp)[1],1:dim(sitestmp)[2],] <- sitestmp


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
    sites[,,k]   <- PP$sites
    
    for(r in seq(nHeight)){
      covar[r,1:nPark,1:nPark,k] <- PP$covar[[r]]
      # electron temperature (Te_par + 2*Te_perp)/3
      param[r,nPark+1,k] <- (PP$param[r,2] + 2*PP$param[r,3] ) / 3 
      # ion temperature (Ti_par + 2*Ti_perp)/3
      param[r,nPark+2,k] <- (PP$param[r,4] + 2*PP$param[r,5] ) / 3
      # ion velocity along magnetic field
      param[r,nPark+3,k] <- PP$param[r,7:9]%*%PP$B[r,]/sqrt(sum(PP$B[r,]**2))
      std[r,nPark+1,k] <- sqrt(c(1,2)%*%PP$covar[[r]][2:3,2:3]%*%c(1,2)/sum(c(1,2)**2))
      std[r,nPark+2,k] <- sqrt(c(1,2)%*%PP$covar[[r]][4:5,4:5]%*%c(1,2)/sum(c(1,2)**2))
      std[r,nPark+3,k] <- sqrt(PP$B[r,]%*%PP$covar[[r]][7:9,7:9]%*%PP$B[r,]/sum(PP$B[r,]**2))
      for( s in seq(nSitesk)){
          # ion velocity seen at site k
          param[r,nPark+2*(s-1)+4,k] <- PP$param[r,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sqrt(sum(PP$intersect[[r]][[s]]$k.ENU**2))
          std[r,nPark+2*(s-1)+4,k] <- sqrt(PP$intersect[[r]][[s]]$k.ENU%*%PP$covar[[r]][7:9,7:9]%*%PP$intersect[[r]][[s]]$k.ENU/sum(PP$intersect[[r]]$k.ENU**2))
          # horizontal ion velocity component seen at site k
          khor <- c(PP$intersect[[r]][[s]]$k.ENU[c(1,2)],0)
          khor <- khor / sqrt(sum(khor**2))
          param[r,nPark+2*s+3,k] <- PP$param[r,7:9]%*%khor
          std[r,nPark+2*s+3,k] <- sqrt(khor%*%PP$covar[[r]][7:9,7:9]%*%khor)
      }
    }
  }
  
  dimnames(param) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]],'Ti','Te','ViB', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep='')),paste(seq(nFile)))
  dimnames(std) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]],'Ti','Te','ViB', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep='')),paste(seq(nFile)))
  dimnames(model) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]],'Ti','Te','ViB', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep='')),paste(seq(nFile)))
  dimnames(covar) <- c(list(paste(seq(nHeight))),lapply(dimnames(PP[["covar"]][[1]]),FUN=function(x,nSites){c(x,'Ti','Te','ViB', paste('ViR',paste(rep(seq(nSites),each=2),c('','hor'),sep=''),sep=''))},nSites=nSites),list(paste(seq(nFile))))

  return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,height=height,time_sec=time_sec,date=date,POSIXtime=POSIXtime,sites=sites,n=nFile,nPar=nPar,nHeight=nHeight,mi=mi,covar=covar))

} # readPP



