readPPIdir <- function(dpath,recursive=F){
# 
# read all PPI result files in the given directory(ies) (and file(s))
# if recursive, search recursively in sub-directories
#
# returns a list containing three-dimensional arrays of plasma parameters and their errors,
# two-dimensional array status and chisqr,
# plus vectors, range, time_sec, date, POSIXtime
# 
# I. Virtanen 2010
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
  param     <- array(NA,dim=c(nRange,nPar,nFile))
  std       <- array(NA,dim=c(nRange,nPar,nFile))
  model     <- array(NA,dim=c(nRange,nPar,nFile))
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

  # read the data from files
  for (k in seq(nFile)){
    load(flist[k])
    param[,,k]   <- PP$param
    std[,,k]     <- PP$std
    model[,,k]   <- PP$model
    chisqr[,k]   <- PP$chisqr
    status[,k]   <- PP$status
    time_sec[k]  <- PP$time_sec
    date[[k]]    <- PP$date
    POSIXtime    <- PP$POSIXtime
    height[,k]   <- PP$height
    if(!all(is.null(PP$azelT))) azelT[,k]    <- PP$azelT
  }


  return(list(param=param,std=std,model=model,chisqr=chisqr,status=status,range=range,height=height,time_sec=time_sec,date=date,POSIXtime=POSIXtime,
              llhT=llhT,llhR=llhR,azelT=azelT,n=nFile,nPar=nPar,nRange=nRange,mi=mi))

} # readPPIdir



