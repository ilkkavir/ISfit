calibrate.eiscat.dynasonde <- function( dataDir , beginTime=c(1970,1,1,0,0,0) , endTime=c(2500,1,1,0,0,0) , region='F' , neMin=0 , multistatic=TRUE , timeplot=FALSE){
# electron density calibration with the EISCAT dynasonde


  btime <- as.double(ISOdate( beginTime[1], beginTime[2], beginTime[3], beginTime[4], beginTime[5], beginTime[6])) + beginTime[6]%%1
  etime <- as.double(ISOdate( endTime[1], endTime[2], endTime[3], endTime[4], endTime[5], endTime[6])) + endTime[6]%%1

  if(multistatic){
      PP <- readPP.3D( dataDir )
  }else{
      PP <- readPP( dataDir )
  }
      
  if(is.null(PP)) stop('No radar data')
  
  PPtimes <- range(PP$time_sec)

  trange <- c(btime,etime)
  trange[1] <- max(PPtimes[1],btime)
  trange[2] <- min(PPtimes[2],etime)

  if(trange[1]>etime) stop('No radar data from the given time interval')
  if(trange[2]<btime) stop('No radar data from the given time interval')


  # count number of days (will be two if the experiment continues over midnight)
  ndays <- as.POSIXlt(trange[2],origin='1970-01-01',tz='ut')[["yday"]]-as.POSIXlt(trange[1],origin='1970-01-01',tz='ut')[["yday"]] + 1


  datavec <- c()
  for( k in seq(ndays) ){
    ttmp <- as.POSIXlt(trange[1] + (k-1)*24*60*60 , origin='1970-01-01' , tz='ut' )
    tstr <- paste( ttmp$year+1900 , ttmp$mon+1 , ttmp$mday , sep='_')
    if( (ttmp$mon<10) & (ttmp$mday<10) )   tstr <- paste( ttmp$year+1900 , paste('0' ,ttmp$mon+1,'_0' , ttmp$mday,sep='') , sep='_')
    if( (ttmp$mon<10) & (ttmp$mday>=10) )  tstr <- paste( ttmp$year+1900 , paste('0' ,ttmp$mon+1 ,sep='') , ttmp$mday , sep='_')
    if( (ttmp$mon>=10) & (ttmp$mday<10) )  tstr <- paste( ttmp$year+1900 , ttmp$mon+1, paste( '0' , ttmp$mday,sep='') , sep='_')
    ofile <- paste(tstr,'.dat',sep='')
    if(!file.exists(ofile)) download.file(paste('http://dynserv.eiscat.uit.no/DD/myque.php?q=select%20dDay,foE,foF2%20from%20tromso.resul',tstr,'%20where%20foE%3E0%20or%20foF2%3E0',sep=''),destfile=ofile)
    datalines <- readLines(ofile)
    datalines <- datalines[10:(length(datalines)-10)]
    datatmp <- as.numeric(unlist(sapply(strsplit(datalines,'  '),function(x){return(x[1:3])})))
    datatmp[seq(1,length(datatmp),by=3)] <- (datatmp[seq(1,length(datatmp),by=3)]-1)*24*3600 + as.double(ISOdate(ttmp$year+1900,1,1,0,0,0))
    datavec <- c(datavec,datatmp)
  }
  dynamat <- matrix(datavec,ncol=3,byrow=TRUE)
  dynamat <- dynamat[which( (dynamat[,1]<=trange[2]) & (dynamat[,1]>=trange[1]) ) , ]

  if(region=='F'){
    dynamat <- dynamat[dynamat[,3]>0,c(1,3)]
  }else if(region=='E'){
    dynamat <- dynamat[dynamat[,2]>0,c(1,2)]
  }else{
    stop('"region" must be either "E" or "F"')
  }

  # the dynasonde data are measured with 120s temporal resolution, take 2 min of radar data around each dynasonde data point
  ndyna <- dim(dynamat)[1]

  if(ndyna==0){
      if(region=='F'){
          stop("No F-region dynasonde data found")
      }else{
          stop("No E-region dynasonda data found")
      }
  }
  radarheight <- radarne <- radarvar<- matrix(ncol=ndyna,nrow=dim(PP$height)[1])

  for(k in seq(ndyna)){
    radarinds <- which( ( PP$time_sec>= (dynamat[k,1]-59) ) & ( PP$time_sec<= dynamat[k,1]+59) )
    if(length(radarinds)>1){
      radarvar[,k] <- 1/rowSums(1/PP$std[,1,radarinds]**2,na.rm=T)
      radarne[,k] <- rowSums(PP$param[,1,radarinds]/PP$std[,1,radarinds]**2,na.rm=T) * radarvar[,k]
      radarheight[,k] <- rowMeans(PP$height[,radarinds])
    }else if (length(radarinds)==1){
      radarvar[,k] <- 1/c(1/PP$std[,1,radarinds]**2)
      radarne[,k] <- c(PP$param[,1,radarinds]/PP$std[,1,radarinds]**2)*radarvar[,k]
      radarheight[,k] <- c(PP$height[,radarinds])
    }
  }


  necompar <- matrix(ncol=4,nrow=ndyna)
  necompar[,1] <- dynamat[,1]
  necompar[,2] <-  (dynamat[,2]/8980)**2*1e18
  
  if(region=='F'){
    for(k in seq(ndyna)){
       necompar[k,3] <- max( radarne[ which( ( radarheight[,k] < 500 ) & ( radarheight[,k] > 200 ) ) , k ] )
       if(is.finite(necompar[k,3])&(!is.na(necompar[k,3]))){
         necompar[k,4] <- radarvar[which(radarne[,k]==necompar[k,3]),k]
       }else{
         necompar[k,4] <- Inf
       }
     }
  }else if( region=='E'){
    for(k in seq(ndyna)){
      necompar[k,3] <- max( radarne[ which( ( radarheight[,k] > 90 ) & ( radarheight[,k] < 200 ) ) , k ] )
       if(is.finite(necompar[k,3])&(!is.na(necompar[k,3]))){
         necompar[k,4] <- radarvar[which(radarne[,k]==necompar[k,3]),k]
       }else{
         necompar[k,4] <- Inf
       }
    }
  }else{
    stop('"region" must be either "E" of "F"')
  }

  necompar <- necompar[necompar[,2]>neMin,]

  if(timeplot){
      plot(necompar[,1],necompar[,2],ylim=c(0,max(necompar[,2:3],na.rm=TRUE,finite=TRUE)))
      points(necompar[,1],necompar[,3],col='red')
      dev.new()
  }

  maxne <- max(necompar[,2:3])

  fitvar <- 1/sum(necompar[,2]**2/necompar[,4],na.rm=T)
  fitres <- sum(necompar[,2]*necompar[,3]/necompar[,4],na.rm=T)*fitvar

  d <- (necompar[,3]-necompar[,2]*fitres)**2/necompar[,4]
  d <- d[!is.na(d)]
  chisqr <- sum(d)/length(d)

  x <- c(0,1e14)

  plot(necompar[,2],necompar[,3],xlim=c(0,maxne)*1.2,ylim=c(0,maxne)*1.2,xlab='Dynasonde',ylab='IS radar',main=paste(sprintf("Ratio = %6.3f, std = %6.3f, chisqr = %6.3f",fitres,sqrt(fitvar),chisqr)))
  lines(x,x*fitres,col='red',lwd=2)
  lines(x,x*(fitres+2*sqrt(fitvar)),col='blue')
  lines(x,x*(fitres-2*sqrt(fitvar)),col='blue')


  
  return(c(fitres,fitvar,chisqr))
  
}
