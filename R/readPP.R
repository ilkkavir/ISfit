readPP <- function(dpath,recursive=F){
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
  # number of range gates
  nRange <- length(PP$range)
  # number of plasma parameters
  nPar   <- dim(PP$param)[2]
  mi     <- PPI_param$mi

  # allocate the necessary arrays
  param     <- array(NA,dim=c(nRange,nPar+4,nFile))
  std       <- array(NA,dim=c(nRange,nPar+4,nFile))
  model     <- array(NA,dim=c(nRange,nPar+4,nFile))
  range     <- PP$range
  height    <- array(NA,dim=c(nRange,nFile))
  status    <- array(NA,dim=c(nRange,nFile))
  chisqr    <- array(NA,dim=c(nRange,nFile))
  time_sec  <- vector(length=nFile,mode='numeric')
  date      <- vector(length=nFile,mode='list')
  POSIXtime <- vector(length=nFile,mode='list')
  llhT      <- PP$llhT
  llhR      <- PP$llhR
  azelT     <- matrix(NA,nrow=2,ncol=nFile)
  covar <- array(NA,dim=c(nRange,nPar+4,nPar+4,nFile))

  # read the data from files
  for (k in seq(nFile)){
    load(flist[k])
    param[,1:nPar,k]   <- PP$param
    std[,1:nPar,k]     <- PP$std
    model[,1:nPar,k]   <- PP$model
    chisqr[,k]   <- PP$chisqr
    status[,k]   <- PP$status
    time_sec[k]  <- PP$time_sec
    date[[k]]    <- PP$date
    POSIXtime    <- PP$POSIXtime
    height[,k]   <- PP$height
    if(!all(is.null(PP$azelT))) azelT[,k]    <- PP$azelT
    for(r in seq(nRange)){
      covar[r,1:nPar,1:nPar,k] <- PP$covar[[r]]
      param[r,nPar+1,k] <- (PP$param[r,2] + 2*PP$param[r,3] ) / 3 
      param[r,nPar+2,k] <- (PP$param[r,4] + 2*PP$param[r,5] ) / 3 
      param[r,nPar+3,k] <- PP$param[r,7:9]%*%PP$B[r,]/sqrt(sum(PP$B[r,]**2))
      param[r,nPar+4,k] <- PP$param[r,7:9]%*%PP$intersect[[r]]$k.ENU/sqrt(sum(PP$intersect[[r]]$k.ENU**2))
      std[r,nPar+1,k] <- sqrt(c(1,2)%*%PP$covar[[r]][2:3,2:3]%*%c(1,2)/sum(c(1,2)**2))
      std[r,nPar+2,k] <- sqrt(c(1,2)%*%PP$covar[[r]][4:5,4:5]%*%c(1,2)/sum(c(1,2)**2))
      std[r,nPar+3,k] <- sqrt(PP$B[r,]%*%PP$covar[[r]][7:9,7:9]%*%PP$B[r,]/sum(PP$B[r,]**2))
      std[r,nPar+4,k] <- sqrt(PP$intersect[[r]]$k.ENU%*%PP$covar[[r]][7:9,7:9]%*%PP$intersect[[r]]$k.ENU/sum(PP$intersect[[r]]$k.ENU**2))
    }
  }
  
  dimnames(param) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]],'Ti','Te','ViB','ViR1'),paste(seq(nFile)))
  dimnames(std) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]],'Ti','Te','ViB','ViR1'),paste(seq(nFile)))
  dimnames(model) <- list(dimnames(PP[["param"]])[[1]],c(dimnames(PP[["param"]])[[2]],'Ti','Te','ViB','ViR1'),paste(seq(nFile)))
  dimnames(covar) <- c(list(paste(seq(nRange))),lapply(dimnames(PP[["covar"]][[1]]),FUN=function(x){c(x,'Ti','Te','ViB','ViR1')}),list(paste(seq(nFile))))

  return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,range=range,height=height,time_sec=time_sec,date=date,POSIXtime=POSIXtime,
              llhT=llhT,llhR=llhR,azelT=azelT,n=nFile,nPar=nPar,nRange=nRange,mi=mi,covar=covar))

} # readPP



