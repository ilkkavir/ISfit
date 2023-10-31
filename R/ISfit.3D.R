ISfit.3D <- function( ddirs='.' , odir='.' ,  heightLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , fitFun=leastSquare.lvmrq , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=10 , absCalib=FALSE , TiIsotropic=TRUE , TeIsotropic=TRUE , recursive=TRUE , aprioriFunction=ISaprioriH , scaleFun=acfscales , siteScales=NULL, calScale=1, MCMCsettings=list( niter=10000 , updatecov=100 , burninlength=5000 , outputlength=5000 ) , maxdev=2 , trueHessian=FALSE , nCores=1 , reverseTime=FALSE , burnin.s=0 ,  ... )
  {

      # 3D incoherent scatter plasma parameter fit using LPI output files in ddirs
      #
      # This is a GUISDAP-style fit, which uses a single ion temperature and collision frequency. Currently 3 ions ( 30.5 , 16.0 , 1.0 )
      #
      #
      #
      # INPUT:
      #   ddirs           a vector of data directory paths, each directory must contain data from exactly one site
      #   odir            Output directory
      #   heightLimits.km analysis height-gate limits. If NA, range gates of the reference site are used
      #   timeRes.s       time resolution (integration time)
      #   beginTime       c(year,month,day,hour,minute,seconds) analysis start time
      #   endTime         c(year,month,day,hour,minute,seconds) analysis end time
      #   absLimit        limit for absolute value of the residual.
      #                   The iteration will not be stopped (unles maxIter is reached) before the residual is below absLimit
      #   diffLimit       Upper limit for fractional change in residual in an iteration step.
      #   maxLambda       maximum Lambda value in Levenberg-Marquardt iteration
      #   maxIter         maximum number of iterations
      #   absCalib        TRUE if the remotes are absolutely calibrated, FALSE to allow for scaling of their calibration coefficients
      #   TiIsotropic     TRUE if ion thermal velocity distribution is modeled as isotropic, FALSE if bi-maxwellian
      #   recursive       logical, should the data directories be searched recursively
      #   scaleFun        function that returns acf scaling factors for each site
      #   siteScales      ACF scales for each site as returned by siteCalib. (Run first with siteScales=NULL, then run siteCalib
      #                   and use its output as siteScales in a second analysis run). This scaling affects only the relative site scales
      #                   actual electron density calibration is done wiht calScale
      #   calScale        additional scaling factor from ionosonde calibration applied to ALL ACF samples
      #   MCMCsettings    a list of input arguments for the modMCMC function
      #   maxdev          maximum angular deviation from the beam centre intersection
      #   trueHessian     logical, calculate the Hessian from finite differences of cost function instead of the direct theory approximation?
      #   nCores          number of parallel processes
      #   reverseTime     logical, should the integration periods be analysed from the last to the first one? (needed for the BAFIM analysis)
      #   burnin.s        duration of a burnin period. If reverseTime=FALSE, te analysis is started burnin.s seconds after the actual start time, runs backwards in time until the start time, and then runs forward until end of the analysis period. If reverseTime=TRUE, the corresponding thing is done at end of the analysis period. (This is needed for the BAFM analysis)
      #
      # OUTPUT:
      #   None, the results are written to files in odir.
      #



      
      
      # create the output directory
      dir.create( odir , recursive=TRUE , showWarnings=FALSE)

      # copy the original function call, it will be stored in each data files
# should also explicitly store every single input argument, as those with default values are not returned by  match.call
      functionCall <- match.call(expand.dots=TRUE)





      # Initialize an empty list for the plasma parameters
      PP <- list()

      # list all data files
      dfiles <- list()
      ddirs <- unique(ddirs)
      nd <- length(ddirs)
      for( n in seq(nd) ){
          dfiles[[n]] <- dir( ddirs[n] , full.names=T , recursive=recursive , pattern="[[:digit:]]*LP.Rdata")
      }
      
      # stop if there is no data
      if ( sum( sapply( dfiles , length ) ) == 0 ) stop( "Could not find any data files." )

      # read timestamps from file names. This is the unix time at end of integration period
      tstamps <- lapply(dfiles,function(x){sapply(x,function(x){as.numeric(substr(rev(unlist(strsplit(x,.Platform[["file.sep"]])))[1],1,13))/1000})})

      # convert beginTime and endTime into unix time format
      bTime <- as.double(ISOdate(beginTime[1],beginTime[2], beginTime[3],beginTime[4], beginTime[5],beginTime[6])) + beginTime[6]%%1
      eTime <- as.double(ISOdate(endTime[1],endTime[2], endTime[3],endTime[4], endTime[5],endTime[6])) + endTime[6]%%1

      if(length(timeRes.s)==1){
          # first integration period from which we actually have data
          iperFirst <- max( 1, floor( ( min( sapply( tstamps , min ) ) - bTime ) / timeRes.s ) )
          
          # last integration period
          iperLast <- ceiling( ( min( max( sapply( tstamps , max ) ) , eTime ) - bTime ) / timeRes.s )

          # stop if analysis end time is before its start time
          if( iperLast < iperFirst ) stop("Last integration period is before the first one.")

          # integration period limits
          iperLimits <- seq( iperFirst - 1 , iperLast ) * timeRes.s + bTime

      }else{
          iperLimits <- timeRes.s
      }














      

      # number of integration periods
      nIper <- length(iperLimits) - 1

      # logical for permanent site selection
      siteselprev <- 0

      # walk through all integration periods, taking into accout reverseTime and burnin.s inputs
      ipers <- seq(nIper)
      if(reverseTime){
          ipers <- rev(ipers)
      }
      nburnin <- 0
      if(burnin.s>0){
          for(k in seq(nIper)){
              if(abs(iperLimits[ipers[k]]-iperLimits[ipers[1]])>burnin.s){
                  break
              }
          }
          nburnin <- k-1
          if(k>1){
              ipers <- c(rev(ipers[2:k]),ipers)
          }
      }

      nnn <- 0
#      for( k in seq( nIper ) ){
      for( k in ipers ){

          nnn <- nnn + 1
          
          # output file name
          resFile <- paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "PP.Rdata" , sep='')
          
          # look for data files from the current integration period
          iperFiles <- lapply( tstamps , function(x,l1,l2){ which( ( x > l1 ) & ( x <= l2 ) ) } , l1=iperLimits[k] , l2=iperLimits[k+1] )

          # load all data and collect it in a list
          nFiles <- sum( sapply( iperFiles , length ) )
          if( nFiles > 0 ){
              
              dlist <- vector( mode='list', length=nd)
              
              for( n in seq( nd ) ){
                  if( length(iperFiles[[n]]) > 0 ){
                      dlist[[n]] <- readACF( dfiles[[n]][ iperFiles[[n]] ] )
                      if(length(dlist[[n]][["nGates"]])!=length(dlist[[n]][["lag.us"]])){
                          dlist[[n]][["nGates"]] <- rep(length(dlist[[n]][["range.km"]]),length(dlist[[n]][["lag.us"]]))
                      }
                      # some sanity checks
                      if( (dlist[[n]][['azelT']][2]<0) | (dlist[[n]][['azelR']][2]<0) | (dlist[[n]][['radarFreq']]<0) ){
                          dlist[[n]] <- list()
                          dlist[[n]][["ACF"]] <- matrix(NA)
                          dlist[[n]][["var"]] <- matrix(NA)
                          dlist[[n]][["lag.us"]] <- NA
                          dlist[[n]][["nGates"]] <- 1
                          dlist[[n]][["range.km"]] <- NA
                          dlist[[n]][["azelT"]] <- c(NA,NA)
                          dlist[[n]][["azelR"]] <- c(NA,NA)
                          dlist[[n]][["llhT"]] <- c(NA,NA,NA)
                          dlist[[n]][["llhR"]] <- c(NA,NA,NA)
                          dlist[[n]][["radarFreq"]] <- NA
                      }
                  }else{
                      dlist[[n]] <- list()
                      dlist[[n]][["ACF"]] <- matrix(NA)
                      dlist[[n]][["var"]] <- matrix(NA)
                      dlist[[n]][["lag.us"]] <- NA
                      dlist[[n]][["nGates"]] <- 1
                      dlist[[n]][["range.km"]] <- NA
                      dlist[[n]][["azelT"]] <- c(NA,NA)
                      dlist[[n]][["azelR"]] <- c(NA,NA)
                      dlist[[n]][["llhT"]] <- c(NA,NA,NA)
                      dlist[[n]][["llhR"]] <- c(NA,NA,NA)
                      dlist[[n]][["radarFreq"]] <- NA
                  }
              }


              # read acf, variance, lag, range, pointing directions, and TX / RX location of each data point
              acf   <- calScale * unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["ACF"]] , n=x[["nGates"]] ) ) ) } ) )
              var   <- calScale**2 * unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["var"]] , n=x[["nGates"]] ) ) ) } ) )
              lag   <- 1e-6*unlist( lapply( dlist , function(x){ return( rep( x[["lag.us"]] , times=x[["nGates"]] ) ) } ) )
              ran   <- 1000*unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] ] ) } , x=x[["range.km"]] , n=x[["nGates"]] ) ) ) } ) )
              sinds <- unlist( lapply( seq(nd) , function( i , x ){ return( rep( i , sum( x[[i]][["nGates"]] ) ) ) } , x=dlist ) )
              
              # number of sites is equal to number of data directories, create a matrix with position and pointing direction information for each site

              sites <- matrix( nrow=nd , ncol=12 )
              for(n in seq(nd)){
                  sites[n,1]     <- n
                  sites[n,2]     <- dlist[[n]][["radarFreq"]]
                  sites[n,3:5]   <- dlist[[n]][["llhT"]]
                  sites[n,6:7]   <- dlist[[n]][["azelT"]]
                  sites[n,8:10]  <- dlist[[n]][["llhR"]]
                  sites[n,11:12] <- dlist[[n]][["azelR"]]
              }


              # if siteScales is null, we will give NULL to ISapriori that will
              # then allow scaling for all others but the reference site
              if(is.null(siteScales)){
                  sScales <- NULL
              }else{
                  # if siteScales was given, we need to pick correct rows from it
                  sScales <- matrix(nrow=nd,ncol=dim(siteScales)[2])
                  for(n in seq(nd)){
                      # select the site that is closest to site n
                      # it must be reasonably close, since the calibration
                      # was done with ISfit.3D results from identical data.
                      # There are small differences because siteCalib rounds the
                      # site information to reasonably accuracy before site
                      # identification
                      sdiffs <- apply( siteScales[,2:12] , FUN=function(x,y){sum(abs(x-y)**2)},MARGIN=1,y=sites[n,2:12])
                      if(!all(is.na(sdiffs))){
                          sScales[n,] <- siteScales[which(sdiffs==min(sdiffs)),]
                          sScales[n,1] <- n
                      }
                  }
              }

              # select the reference site
              # if there is data for only one site this is trivial
              # if absCalib==TRUE we can just pick any site
              if( nd==1 | absCalib | sum(!is.na(rowSums(sites)))==1 ){
                  refsite <- which(!is.na(rowSums(sites)))[1]
              }else{ # otherwise look for a monostatic site

                  refsite <- which( apply( sites[,3:7] == sites[,8:12] , FUN=all , MARGIN=1 ) )

                  # if there are more than one candidate we will first check if the antenna was
                  # just moving slightly during the integration
                  # (this was left from an earlier, slightly different, version. Should never happen the present version...)
                  if(length(refsite)>1){
                      # largest variation in azimuth / elevation
                      maxPointDiff <- max( apply( sites[refsite,c(6,7,11,12)] ,
                                                 FUN=function(x){diff(range(x))} ,
                                                 MARGIN=2) )
                      # this should not be too common in well-designed experiments
                      # simply proceed by picking the first candidate
                      if(maxPointDiff < .1 ){
                          refsite <- refsite[1]
                      }
                      
                  }

                 # just pick the first available site and warn the user
                  if(length(refsite) != 1){
                      if(length(refsite)==0){
                          cat("Could not find a monostatic site, using the first bistatic as a reference\n")
                          refsite <- which(!is.na(rowSums(sites)))[1] #sometimes eiscat data contains NaNs as coordinates...
                      }else{
                          cat("Found two or more monostatic sites, using the first one as a reference\n")
                          refsite <- refsite[1]
                      }
                      
                  }
                  
              }

              
              # scale all data with values returned by scaleFun
              dscales <- scaleFun( sites )
              for( n in seq( nd ) ){
                  spoints <- which(sinds==n)
                  acf[spoints] <- acf[spoints] * dscales[n,6]
                  if(!is.na(dscales[n,5])){
                      if( dscales[n,5] ) acf[spoints] <- Conj(acf[spoints])
                  }
                  var[spoints] <- var[spoints] * dscales[n,6]**2
              }
              

              # a time vector converted from iperLimits
              t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='utc')
              date <- c(t$year+1900,t$mon+1,t$mday,t$hour,t$min,t$sec)
              tprev <- as.POSIXlt( iperLimits[k] , origin='1970-01-01' , tz='utc')
              dateprev <- c(tprev$year+1900,tprev$mon+1,tprev$mday,tprev$hour,tprev$min,tprev$sec)
              
              # height gate limits
              if( all( is.na( heightLimits.km ) ) ){
                  rlims <- sort( unique( ran[ which( sinds==refsite ) ] ) )
                  rlims <- c( rlims , max(rlims) + 1)
                  hlims <- rlims
                  for( h in seq(length(hlims))) hlims[h] <- range2llh(r=rlims[h],llhT=sites[refsite,3:5],llhR=sites[refsite,8:10],azelT=sites[refsite,6:7])['h']
              }else{
                  hlims <- unique( heightLimits.km )*1000
              }
              

              # IRI does not work above 70 km
              hlims <- hlims[hlims>=70000]
              
              
              nh <- length( hlims ) - 1
              
              covar <- intersect <- list()
              for(h in seq(nh)) covar[[h]] <- matrix(ncol=12+nd,nrow=12+nd)
              
              model <- std <- param <- matrix(ncol=12+nd,nrow=nh)
              latitude <- longitude <- height <- status <- chisqr <- rep(-1,nh)
              B <- B2 <- matrix(ncol=3,nrow=nh)
              MCMC <- vector(length=nh,mode='list')
              
              # convert all ranges to latitude, longitude, height
              llh <- matrix(nrow=length(ran),ncol=3)
              llhlist <- mclapply(seq(length(ran)) , FUN=range2llhParFun ,  ran=ran , sites=sites , sinds=sinds , mc.cores=nCores)
              for( dind in seq(length(ran))){
#                  if(!is.na(ran[dind])) llh[dind,] <- range2llh( r=ran[dind] , llhT=sites[sinds[dind],3:5] , llhR=sites[sinds[dind],8:10] , azelT=sites[sinds[dind],6:7])
                  if(!is.na(ran[dind])) llh[dind,] <- llhlist[[dind]]
              }

              # a list for site indices contributing at each height
              contribSites <- apriori <- vector(mode='list',length=nh)#list()
              
              
              










              acf.gate <- var.gate <- lag.gate <- ran.gate <- llh.gate <- sites.gate <- list()
              acf.site <- var.site <- lag.site <-  nlags.site <- ind.site <- aSite  <- kSite  <- kBsite <- Brotmat <- list()
              initParam <- limitParam  <- list()
              fitGate <- rep(T,nh)
              
              for( h in seq( nh ) ){

                  # Initial values, these will be immediately updated if any data is found
                  height[h]      <- sum(hlims[h:(h+1)])/2000
                  latitude[h]    <- sites[refsite,3]
                  longitude[h]   <- sites[refsite,4]
                  Btmp           <- igrf(date=date[1:3],lat=latitude[h],lon=longitude[h],height=height[h],isv=0,itype=1)
                  B[h,]          <- c(Btmp$y,Btmp$x,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards east,
                  dimnames(covar[[h]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
                  
                  # if there are actual data, all initializations will be updated


                  gateinds       <- which( ( llh[,3] >= hlims[h] ) & (llh[,3] < hlims[h+1]) )
                  
                  acf.gate[[h]]       <- acf[ gateinds ]
                  var.gate[[h]]       <- var[ gateinds ]
                  lag.gate[[h]]       <- lag[ gateinds ]
                  ran.gate[[h]]       <- ran[ gateinds ]
                  llh.gate[[h]]       <- llh[ gateinds , ]
                  sites.gate[[h]]     <- sites[ sinds[gateinds] , ]
                  
                  if(!is.matrix(llh.gate[[h]])) llh.gate[[h]] <- matrix(llh.gate[[h]],nrow=1)
                  if(!is.matrix(sites.gate[[h]])) sites.gate[[h]] <- matrix(sites.gate[[h]],nrow=1)
                  
                  # remove NA values
                  nainds         <- is.na(acf.gate[[h]]) | is.na(var.gate[[h]])

                  if(any(!nainds)){

                      acf.gate[[h]]       <- acf.gate[[h]][ !nainds ]
                      var.gate[[h]]       <- var.gate[[h]][ !nainds ]
                      lag.gate[[h]]       <- lag.gate[[h]][ !nainds ]
                      llh.gate[[h]]       <- llh.gate[[h]][ !nainds , ]
                      sites.gate[[h]]     <- sites.gate[[h]][ !nainds , ]
                      
                      if(!is.matrix(llh.gate[[h]])) llh.gate[[h]] <- matrix(llh.gate[[h]],nrow=1)
                      if(!is.matrix(sites.gate[[h]])) sites.gate[[h]] <- matrix(sites.gate[[h]],nrow=1)

                      # coordinates of the measurement volume
                      llhTarget      <- colMeans( llh.gate[[h]] )
                      height[h]      <- llhTarget[3] / 1000
                      latitude[h]    <- llhTarget[1]
                      longitude[h]   <- llhTarget[2]
                      
                      
                      
                      # beam intersections
                      intersect[[h]] <- list()
                      gainR <- aSite[[h]] <- nlags.site[[h]] <-  rep(NA,nd)
                      kSite[[h]] <- list()
                      lag.site[[h]] <- list()
                      acf.site[[h]] <- list()
                      var.site[[h]] <- list()
                      ind.site[[h]] <- list()
                      for( s in seq(nd) ){
                          
                          if(any(is.na(sites[s,]))){
                              gainR[s] <- NA
                              aSite[[h]][s] <- NA
                              kSite[[h]][[s]] <- c(NA,NA,NA)
                          }else{
                              # the beam widths and antenna types are stored in dscales
                              intersect[[h]][[s]] <- beamIntersection( llhT=sites[s,3:5] , llhR=sites[s,8:10] , azelT=sites[s,6:7] , azelR=sites[s,11:12] , fwhmT=dscales[s,1] , fwhmR=dscales[s,3] , phArrT=dscales[s,2]>0 , phArrR=dscales[s,4]>0 , freq.Hz=sites[s,2] )
                              
                              # conversion from lat, lon, height to range in this gate
                              rs1 <- height2range( llhT=sites[s,3:5] , azelT=sites[s,6:7] , llhR=sites[s,8:10] , h=hlims[h] )
                              rs2 <- height2range( llhT=sites[s,3:5] , azelT=sites[s,6:7] , llhR=sites[s,8:10] , h=hlims[h+1] )
                              
                              # gain integral
                              gainR[s] <- gategain( intersect[[h]][[s]] , c(rs1,rs2) , maxdev=maxdev)
                              
                              # scattering angle
                              aSite[[h]][s] <- intersect[[h]][[s]][["phi"]]
                              
                              # scattering wave vector
                              kSite[[h]][[s]] <- intersect[[h]][[s]][["k.ENU"]]
                          }

                          if(is.na(gainR[s])){
                              lag.site[[h]][[s]] <- c()
                              nlags.site[[h]][[s]] <- 0
                              acf.site[[h]][[s]] <- c()
                              var.site[[h]][[s]] <- c()
                              ind.site[[h]][[s]] <- c()
                          }else{
                              
                              # data points from this site
                              data.site <- which(sites.gate[[h]][,1]==s)
                              
                              # lags measured at this site
                              lag.site[[h]][[s]] <- unique(lag.gate[[h]][data.site])
                              
                              nlags.site[[h]][s] <- length(lag.site[[h]][[s]])
                              
                              # average data points from each lag value
                              acf.site[[h]][[s]] <- rep(0+0i,nlags.site[[h]][s])
                              var.site[[h]][[s]] <- rep(0,nlags.site[[h]][s])
                              for( l in seq(nlags.site[[h]][s])){
                                  lagind <- which(lag.site[[h]][[s]][l]==lag.gate[[h]][data.site])
                                  for( lind in lagind ){
                                      acf.site[[h]][[s]][l] <- acf.site[[h]][[s]][l] + acf.gate[[h]][data.site][lind]/var.gate[[h]][data.site][lind]
                                      var.site[[h]][[s]][l] <- var.site[[h]][[s]][l] + 1/var.gate[[h]][data.site][lind]
                                  }
                              }
                              var.site[[h]][[s]] <- 1/var.site[[h]][[s]]
                              acf.site[[h]][[s]] <- acf.site[[h]][[s]] * var.site[[h]][[s]]
                              
                              acf.site[[h]][[s]] <- acf.site[[h]][[s]] / gainR[s]
                              var.site[[h]][[s]] <- var.site[[h]][[s]] / gainR[s]**2
                              
                              ind.site[[h]][[s]] <- rep(s,nlags.site[[h]][s])
                          }
                          
                      }


                  }else{
                      fitGate[h] <- FALSE
                  }
              }



              # magnetic field from igrf
              for(h in seq(nh)){
                  if(fitGate[h]){
                  
                      if(sum(nlags.site[[h]])>0){
                          
                          # magnetic field direction
                          Btmp           <- igrf(date=date[1:3],lat=latitude[h],lon=longitude[h],height=height[h],isv=0,itype=1)
                          # the model has y-axis to east and z-axis downwards, we have x towards east and z upwards
                          B[h,]          <- c(Btmp$y,Btmp$x,-Btmp$z)
                      }
                  }

                  # rotate all k-vectors to magnetic coordinates
                  # the z-axis is along B but upwards
                  Bz <- -B[h,]/sqrt(sum(B[h,]**2))
                  # horizontal component of B
                  Bhor <- c(B[h,1:2],0)
                  # geomagnetic east is perpendicular both to B and Bhor
                  Bx <- radarPointings:::vectorProduct.cartesian(B[h,],Bhor)
                  Bx <- Bx / sqrt(sum(Bx**2))
                  # y completes the right-handed coordinate system
                  By <- radarPointings:::vectorProduct.cartesian(Bx,B[h,])
                  By <- By / sqrt(sum(By**2))

                  Brotmat[[h]] <- matrix(c(Bx,By,Bz),ncol=3,byrow=F)

                  
                  # checking that the conversion was done correctly
                  B2[h,] <- B[h,]%*%Brotmat[[h]]

                  kBsite[[h]] <- list()

                  if(fitGate[h]){
                      for(ss in seq(length(kSite[[h]]))){
                          kBsite[[h]][[ss]] <- kSite[[h]][[ss]]%*%Brotmat[[h]]
#                          print(kBsite[[h]][[s]])
                      }
                  }
              }

              # the prior model
              apriori <- aprioriFunction( PP=PP , date=date , dateprev=dateprev , latitude=latitude , longitude=longitude , height=height , nSite=nd , nIon=3 , absCalib=absCalib , TiIsotropic=TiIsotropic , TeIsotropic=TeIsotropic , refSite=refsite , siteScales=sScales , B=B2 ,  nCores=nCores , resFile=resFile , updateFile=ifelse(nnn>nburnin,TRUE,FALSE) , ... )

              
              # copy the model/initial values in a matrix for backward compatibility
              for (h in seq(nh)){
                  model[h,] <- apriori[[h]][["aprioriParam"]]
              }

              ## # initialize flipchem if it will be used
              ## fc <- NULL
              ## chisqrFun <- chiSquare
              ## if(!is.null(flipchemStd)){
              ##     library(reticulate)
              ##     Sys.setenv(RETICULATE_AUTOCREATE_PACKAGE_VENV="no")
              ##     fcfile <- system.file('python','startflipchem.py',package='ISfit')
              ##     source_python(fcfile)
              ##     idate <- as.integer(date)
              ##     fc <- startflipchem(idate[1],idate[2],idate[3],idate[4],idate[5],idate[6])
              ##     chisqrFun <- chisqrFlipchem
              ##     trueHessian <- TRUE # chisqrFlipchem works only when the Hessian matrix is calculated from finite differences of chi square
              ## }
              
              # run the actual iterative fit in parallel              
              fitpar <- mclapply(seq(nh),FUN=ISparamfitParallel,
                                 acf             = acf.site,
                                 var             = var.site,
                                 lags            = lag.site,
                                 nData           = nlags.site,
                                 fSite           = sites[,2],
                                 aSite           = aSite,
                                 kSite           = kBsite,
                                 iSite           = ind.site,
                                 B               = B2,
                                 apriori         = apriori,
                                 directTheory    = ISdirectTheory,
                                 absLimit        = absLimit,
                                 diffLimit       = diffLimit,
                                 scaleFun        = scaleParams,
                                 maxLambda       = maxLambda,
                                 maxIter         = maxIter,
                                 fitFun          = fitFun,
                                 MCMCsettings    = MCMCsettings,
                                 trueHessian     = trueHessian,
                                 heights         = height,
                                 latitude        = latitude,
                                 longitude       = longitude,
                                 fitGate         = fitGate,
                                 mc.cores=nCores
                                 )
              
              ## fitpar <- list()
              ## for(h in seq(nh)){
              ##     fitpar[[h]] <- ISparamfitParallel(h,
              ##                    acf             = acf.site,
              ##                    var             = var.site,
              ##                    lags            = lag.site,
              ##                    nData           = nlags.site,
              ##                    fSite           = sites[,2],
              ##                    aSite           = aSite,
              ##                    kSite           = kBsite,
              ##                    iSite           = ind.site,
              ##                    B               = B2,
              ##                    apriori         = apriori,
              ##                    directTheory    = ISdirectTheory,
              ##                    absLimit        = absLimit,
              ##                    diffLimit       = diffLimit,
              ##                    scaleFun        = scaleParams,
              ##                    maxLambda       = maxLambda,
              ##                    maxIter         = maxIter,
              ##                    fitFun          = fitFun,
              ##                    MCMCsettings    = MCMCsettings,
              ##                    trueHessian     = trueHessian,
              ##                    heights         = height,
              ##                    latitude        = latitude,
              ##                    longitude       = longitude,
              ##                    fitGate         = fitGate)#,
              ##     }
              
              # scale the results back to physical units and add some metadata
              for(h in seq(nh)){
                  if(fitGate[h]){
                      param[h,] <- scaleParams( fitpar[[h]]$param , scale=apriori[[h]]$parScales , inverse=TRUE )
                      covar[[h]] <- scaleCovar( fitpar[[h]]$covar , scale=apriori[[h]]$parScales , inverse=TRUE)
                      std[h,]   <- sqrt(diag(covar[[h]]))
                      chisqr[h] <- fitpar[[h]][["chisqr"]]
                      status[h] <- fitpar[[h]][["fitStatus"]]
                      dimnames(covar[[h]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
                      
                      contribSites[[h]] <- unique(unlist(ind.site[[h]]))
                      if(!is.null(fitpar[[h]]$MCMC)){
                          MCMC[[h]] <- fitpar[[h]]$MCMC
                          for( sr in seq(dim(MCMC[[h]][["pars"]][1]))) MCMC[[h]][["pars"]][k,] <-  scaleParams( MCMC[[h]][["pars"]][k,] , parScales[[h]] , inverse=T )
                          MCMC[[h]][["bestpar"]] <- scaleParams( MCMC[[h]][["bestpar"]] , apriori[[h]]$parScales , inverse=TRUE )
                          MCMC[[h]][["pars"]] <- t( apply( MCMC[[h]][["pars"]] , FUN=scaleParams , MARGIN=1  , scale=apriori[[h]]$parScales , inverse=TRUE) )
                          dimnames(MCMC[[h]][["pars"]]) <- list( NULL , c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')) )
                          names(MCMC[[h]][["bestpar"]]) <- c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep=''))
                      }
                  }
              }
              



              time_sec <- iperLimits[k+1]
              POSIXtime <- as.POSIXlt(time_sec,origin='1970-01-01',tz='utc')
              std[is.na(std)] <- Inf
              
              
              dimnames(param) <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
              dimnames(std)   <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
              dimnames(model)   <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
              names(height) <- paste('gate',seq(nh),sep='')
              names(latitude) <- paste('gate',seq(nh),sep='')
              names(longitude) <- paste('gate',seq(nh),sep='')
              dimnames(B) <- list(paste('gate',seq(nh),sep=''),c('x','y','z'))
              
              
              # save the results to file
              PP <- list(param=param,std=std,model=model,chisqr=chisqr,status=status,time_sec=time_sec,date=date,POSIXtime=POSIXtime,height=height,latitude=latitude,longitude=longitude,sites=sites,intersect=intersect,covar=covar,B=B,heightLimits.km=hlims/1000,contribSites=contribSites,mIon=c(30.5,16.0,1.0),MCMC=MCMC,timeLimits.s=iperLimits[k:(k+1)],functionCall=functionCall,apriori=apriori,resFile=resFile , resDir=odir,ViCoordinates='ENUmagnetic')
              if(nnn>nburnin){
                  save( PP , file=file.path(odir,resFile) )
              }

              cat(iperLimits[k+1],'\n')
#              print(PP$param[,1:12])
              
          }
          
      }
  }

