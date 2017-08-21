readPP <- function(dpath,recursive=F,...){
#
#
# read plasma parameters from files
# This is a very simple version that does not calculate
# projections to sites etc. Use readPP.3D for more advanced
# treatment of parameters.
#
# I. Virtanen 2010, 2013
#

  if(is.null(dpath))   return(NULL)
  if(length(dpath)==0) return(NULL)

  fpath <- dpath[!file.info(dpath)$isdir]
  fpath <- fpath[length(grep('PP.Rdata',fpath))>0]
  dpath <- dpath[file.info(dpath)$isdir]

  # list the plasma parameter files
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

  # allocate the necessary arrays
  param     <- array(NA,dim=c(nHeight,nPar,nFile))
  std       <- array(NA,dim=c(nHeight,nPar,nFile))
  model     <- array(NA,dim=c(nHeight,nPar,nFile))
  height    <- array(NA,dim=c(nHeight,nFile))
  status    <- array(NA,dim=c(nHeight,nFile))
  chisqr    <- array(NA,dim=c(nHeight,nFile))
  time_sec  <- vector(length=nFile,mode='numeric')
  timeLimits <- array(NA,dim=c(2,nFile))
  date      <- vector(length=nFile,mode='list')
  POSIXtime <- vector(length=nFile,mode='list')

  # read the data from files
  for (k in seq(nFile)){
    load(flist[k])
    nhk <- length(PP$height)
    if( nhk > nHeight ){
        partmp <- param
        stdtmp <- std
        modeltmp <- model
        heighttmp <- height
        statustmp <- status
        chisqrtmp <- chisqr
        param <- std <- model <- array(NA,dim=c(nhk,nPar,nFile))
        height <- status <- chisqr <- array(NA,dim=c(nhk,nFile))
        par[1:nHeight,,] <- partmp
        std[1:nHeight,,] <- stdtmp
        model[1:nHeight,,] <- modeltmp
        height[1:nHeight] <- heighttmp
        status[1:nHeight] <- statustmp
        chisqr[1:nHeight] <- chisqrtmp
        nHeight <- nhk
    }
    param[1:nhk,,k]   <- PP$param
    std[1:nhk,,k]     <- PP$std
    model[1:nhk,,k]   <- PP$model
    chisqr[1:nhk,k]   <- PP$chisqr
    status[1:nhk,k]   <- PP$status
    time_sec[k]  <- PP$time_sec
    timeLimits[,k]      <- PP$time_limits
    date[[k]]    <- PP$date
    POSIXtime[[k]]<- PP$POSIXtime
    height[1:nhk,k]   <- PP$height
  }

  dimnames(param) <- list(paste("gate",seq(nHeight),sep=''),dimnames(PP[["param"]])[[2]],seq(nFile))
  dimnames(std) <- list(paste("gate",seq(nHeight),sep=''),dimnames(PP[["param"]])[[2]],seq(nFile))
  dimnames(model) <- list(paste("gate",seq(nHeight),sep=''),dimnames(PP[["param"]])[[2]],seq(nFile))

  # warn the user about dangers of this function..
  warning("This function does not make any checks about data reliability, use readPP.3D instead if you are not sure what you are doing.")

  return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,height=height,time_sec=time_sec,timeLimits=timeLimits,date=date,POSIXtime=POSIXtime, n=nFile,nPar=nPar,nHeight=nHeight))

} # readPP



