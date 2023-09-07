## file:readPlasmaLines.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
##
## An average of all plasma line frequencies in directories/files listed in dpath
##
## Arguments:
##  dpath     data path(s), a directory, a vector of file names
##            or a combination of these
##  ranges    range gate selection
##  stdThrsh  standard deviation threshold, data points with standard
##            deviation larger than stdThrsh are not used in the average
##
## Returns:
##  PL a PL list of the averaged data
## 

readPlasmaLines <- function( dpath , ranges=NULL , stdThrsh=Inf )
  {
    
    if(is.null(dpath))   return(NULL)
    if(length(dpath)==0) return(NULL)
    
    fpath <- dpath[!file.info(dpath)$isdir]
    fpath <- fpath[length(grep('LP.Rdata',fpath))>0]
    dpath <- dpath[file.info(dpath)$isdir]
    
    # list the result files
    flist <- dir(dpath[1],pattern='PL.Rdata',full.names=TRUE)
    if(length(dpath)>1){
      for(k in seq(2,length(dpath))){
        flist <- c(flist,dir(dpath[k],pattern='PL.Rdata',full.names=TRUE))
      }
    }
    flist <- c(flist,fpath)
    
    if(length(flist)==0) return(NULL)
    
    # number of data files
    nFile  <- length(flist)
    
    # read the last data file
    load(flist[nFile])

    # number of range-gates
    if(is.null(ranges)){
      nRange   <- length(PL$range)
      ranges   <- seq(nRange)
    }else{
      ranges <- ranges[ranges<length(PL$range)]
      ranges <- ranges[ranges>0]
      nRange <- length(ranges)
    }
    

    # allocate the necessary arrays
    fupave        <- rep(0,nRange)
    varupave      <- rep(0,nRange)
    fdownave      <- rep(0,nRange)
    vardownave    <- rep(0,nRange)
    azelT         <- c(0,0)
    azelR         <- c(0,0)
    range.km      <- PL[["range.km"]] * 0
    
    # read the data and average
    for( fn in flist){
      load(fn)
      indup <- ( sqrt(PL$varup[ranges]) < stdThrsh ) & (!is.na(PL$fup[ranges,lags]))
      indup[is.na(indup)] <- FALSE
      fupave[indup] <- fupave[indup] + (PL$fup[ranges] / PL$varup[ranges])[indup]
      varupave[indup] <- varupave[indup] + (1/PL$varup[ranges])[indup]

      inddown <- ( sqrt(PL$vardown[ranges]) < stdThrsh ) & (!is.na(PL$fdown[ranges,lags]))
      inddown[is.na(inddown)] <- FALSE
      fdownave[inddown] <- fdownave[inddown] + (PL$fdown[ranges] / PL$vardown[ranges])[inddown]
      vardownave[inddown] <- vardownave[inddown] + (1/PL$vardown[ranges])[inddown]
      
      azelT <- azelT + ACF[["azelT"]]
      azelR <- azelR + ACF[["azelR"]]
      range <- range + ACF[["range"]]
      range.km <- range.km + ACF[["range.km"]]
    }
    varupave <- 1/varupave
    fupave <- fupave * varupave
    vardownave <- 1/vardownave
    fdownave <- fdownave * vardownave
    azelT <- azelT / nFile
    azelR <- azelR / nFile
    range <- range / nFile
    range.km <- range.km / nFile
    
    # replace the values in the last data file with the averaged ones and return
    PL[["fup"]] <- fupave
    PL[["varup"]] <- varupave
    PL[["fdown"]] <- fdownave
    PL[["vardown"]] <- vardownave
    PL[["range.km"]] <- range.km[ranges]
    PL[["nGates"]] <- nRange
    PL[["azelT"]] <- azelT
    PL[["azelR"]] <- azelR
    
    return(PL)
    
  }
