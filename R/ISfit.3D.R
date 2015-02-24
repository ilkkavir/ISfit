ISfit.3D <- function( ddirs='.' , odir='.' ,  heightLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , fitFun=leastSquare.lvmrq , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=10 , plotTest=FALSE , plotFit=FALSE , absCalib=FALSE , TiIsotropic=TRUE , TeIsotropic=TRUE , recursive=TRUE , aprioriFunction=ISapriori , scaleFun=acfscales , siteScales=NULL, calScale=1, MCMCsettings=list( niter=10000 , updatecov=100 , burninlength=5000 , outputlength=5000 ) , maxdev=2 , trueHessian=FALSE , nCores=1)
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
      #   plotTest        if TRUE, the direct theory and measurements are plotted at every single iteration step
      #   absCalib        TRUE if the remotes are absolutely calibrated, FALSE to allow for scaling of their calibration coefficients
      #   TiIsotropic     TRUE if ion thermal velocity distribution is modeled as isotropic, FALSE if bi-maxwellian
      #   recursive       logical, should the data directories be searched recursively
      #   scaleFun        function that returns acf scaling factors for each site
      #   siteScales      ACF scales for each site as returned by siteCalib. (Run first with siteScales=NULL, then run siteCalib
      #                   and use its output as siteScales in a second analysis run). This scaling affects only the relative site scales
      #                   actual electron density calibration is done wiht calScale
      #   calScale        additional scaling factor from ionosonde calibration applied to ALL ACF samples
      #   MCMCsettings    a list of input arguments for the modMCMC function
      #
      # OUTPUT:
      #   None, the results are written to files in odir.
      #

      if(nCores>1){
          if(plotTest) warning("Cannot use graphics with nCores > 1, setting plotTest=FALSE")
          if(plotFit) warning("Cannot use graphics with nCores > 1, setting plotFit=FALSE")
          plotTest <- FALSE
          plotFit <- FALSE
      }

      # create the output directory
      dir.create( odir , recursive=TRUE , showWarnings=FALSE)

      # copy the original function call, it will be stored in each data files
      functionCall <- match.call(expand.dots=TRUE)

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

      # walk through all integration periods
      for( k in seq( nIper ) ){
          runcmd <- expression(
              {
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
                      # if there is only one site this is trivial
                      # if absCalib==TRUE we can just pick any site
                      if( nd==1 | absCalib ){
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
                          # if there are none or several candidates we will need to ask the user
                          if( length(refsite) != 1){
                              if( length(refsite)==0){
                                  if(!siteselprev){
                                      cat('Could not find a monostatic site,\n')
                                      cat('available sites are:\n')
                                      cat('Site        TXfreq       TXlat       TXlon    TXheight        TXaz       TXele       RXlat       RXlon    RXheight        RXaz       RXele\n')
                                      for(s in seq(nd)){
                                          cat(sprintf("[%2i]",s),sprintf(" %10.2f",sites[s,2:12]),sprintf("\n"))
                                      }
                                      refsite <- as.numeric(readline("Please select number of the reference site: "))
                                      testnum <- as.numeric( readline("Use the same selection in all intgration periods? (1=yes, 0=no) :  ") )
                                      if(testnum) siteselprev <- refsite
                                  }else{
                                      cat('Could not find a monostatic site,\n')
                                      cat('available sites are:\n')
                                      cat('Site        TXfreq       TXlat       TXlon    TXheight        TXaz       TXele       RXlat       RXlon    RXheight        RXaz       RXele\n')
                                      for(s in seq(nd)){
                                          cat(sprintf("[%2i]",s),sprintf(" %10.2f",sites[s,2:12]),sprintf("\n"))
                                      }
                                      cat("Using previous selection (",siteselprev,")\n")
                                      refsite <- siteselprev
                                  }
                              }else{
                                  if(!siteselprev){
                                      cat('Found several monostatic sites,\n')
                                      cat('available sites are:\n')
                                      cat('Site        TXfreq       TXlat       TXlon    TXheight        TXaz       TXele       RXlat       RXlon    RXheight        RXaz       RXele\n')
                                      for(s in refsite){
                                          cat(sprintf("[%2i]",s),sprintf(" %10.2f",sites[s,2:12]),sprintf("\n"))
                                      }
                                      refsite <- as.numeric(readline("Please select number of the reference site: "))
                                      siteselprev <- refsite
                                  }else{
                                      cat('Found several monostatic sites,\n')
                                      cat('available sites are:\n')
                                      cat('Site        TXfreq       TXlat       TXlon    TXheight        TXaz       TXele       RXlat       RXlon    RXheight        RXaz       RXele\n')
                                      for(s in refsite){
                                          cat(sprintf("[%2i]",s),sprintf(" %10.2f",sites[s,2:12]),sprintf("\n"))
                                      }
                                      cat("Using previous selection (",siteselprev,")\n")
                                      refsite <- siteselprev
                                  }
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
                      t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='ut')
                      date <- c(t$year+1900,t$mon+1,t$mday,t$hour,t$min,t$sec)

                      # height gate limits
                      if( all( is.na( heightLimits.km ) ) ){
                          rlims <- sort( unique( ran[ which( sinds==refsite ) ] ) )
                          rlims <- c( rlims , max(rlims) + 1)
                          hlims <- rlims
                          for( h in seq(length(hlims))) hlims[h] <- range2llh(r=rlims[h],llhT=sites[refsite,3:5],llhR=sites[refsite,8:10],azelT=sites[refsite,6:7])['h']
                      }else{
                          hlims <- unique( heightLimits.km )*1000
                      }
                      nh <- length( hlims ) - 1

                      covar <- intersect <- list()
                      for(h in seq(nh)) covar[[h]] <- matrix(ncol=12+nd,nrow=12+nd)

                      model <- std <- param <- matrix(ncol=12+nd,nrow=nh)
                      latitude <- longitude <- height <- status <- chisqr <- rep(-1,nh)
                      B <- matrix(ncol=3,nrow=nh)
                      MCMC <- vector(length=nh,mode='list')

                      # convert all ranges to latitude, longitude, height
                      llh <- matrix(nrow=length(ran),ncol=3)
                      for( dind in seq(length(ran))){
                          if(!is.na(ran[dind])) llh[dind,] <- range2llh( r=ran[dind] , llhT=sites[sinds[dind],3:5] , llhR=sites[sinds[dind],8:10] , azelT=sites[sinds[dind],6:7])
                      }

                      # a list for site indices contributing at each height
                      contribSites <- apriori <- list()

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

                          acf.gate       <- acf[ gateinds ]
                          var.gate       <- var[ gateinds ]
                          lag.gate       <- lag[ gateinds ]
                          ran.gate       <- ran[ gateinds ]
                          llh.gate       <- llh[ gateinds , ]
                          sites.gate     <- sites[ sinds[gateinds] , ]

                          if(!is.matrix(llh.gate)) llh.gate <- matrix(llh.gate,nrow=1)
                          if(!is.matrix(sites.gate)) sites.gate <- matrix(sites.gate,nrow=1)

                          # remove NA values
                          nainds         <- is.na(acf.gate) | is.na(var.gate)

                          if(any(!nainds)){

                              acf.gate       <- acf.gate[ !nainds ]
                              var.gate       <- var.gate[ !nainds ]
                              lag.gate       <- lag.gate[ !nainds ]
                              llh.gate       <- llh.gate[ !nainds , ]
                              sites.gate     <- sites.gate[ !nainds , ]

                              if(!is.matrix(llh.gate)) llh.gate <- matrix(llh.gate,nrow=1)
                              if(!is.matrix(sites.gate)) sites.gate <- matrix(sites.gate,nrow=1)

                              # coordinates of the measurement volume
                              llhTarget      <- colMeans( llh.gate )
                              height[h]      <- llhTarget[3] / 1000
                              latitude[h]    <- llhTarget[1]
                              longitude[h]   <- llhTarget[2]


                             # beam intersections
                              intersect[[h]] <- list()
                              gainR <- aSite <- nlags.site <-  rep(NA,nd)
                              kSite <- list()
                              lag.site <- list()
                              acf.site <- list()
                              var.site <- list()
                              ind.site <- list()
                              for( s in seq(nd) ){

                                  if(any(is.na(sites[s,]))){
                                      gainR[s] <- NA
                                      aSite[s] <- NA
                                      kSite[[s]] <- c(NA,NA,NA)
                                  }else{
                                      # the beam widths and antenna types are stored in dscales
                                      intersect[[h]][[s]] <- beamIntersection( llhT=sites[s,3:5] , llhR=sites[s,8:10] , azelT=sites[s,6:7] , azelR=sites[s,11:12] , fwhmT=dscales[s,1] , fwhmR=dscales[s,3] , phArrT=dscales[s,2]>0 , phArrR=dscales[s,4]>0 , freq.Hz=sites[s,2] )

                                      # conversion from lat, lon, height to range in this gate
                                      rs1 <- height2range( llhT=sites[s,3:5] , azelT=sites[s,6:7] , llhR=sites[s,8:10] , h=hlims[h] )
                                      rs2 <- height2range( llhT=sites[s,3:5] , azelT=sites[s,6:7] , llhR=sites[s,8:10] , h=hlims[h+1] )

                                      # gain integral
                                      gainR[s] <- gategain( intersect[[h]][[s]] , c(rs1,rs2) , maxdev=maxdev)

                                      # scattering angle
                                      aSite[s] <- intersect[[h]][[s]][["phi"]]

                                      # scattering wave vector
                                      kSite[[s]] <- intersect[[h]][[s]][["k.ENU"]]
                                  }
                                  if(is.na(gainR[s])){
                                      lag.site[[s]] <- c()
                                      nlags.site[[s]] <- 0
                                      acf.site[[s]] <- c()
                                      var.site[[s]] <- c()
                                      ind.site[[s]] <- c()
                                  }else{

                                      # data points from this site
                                      data.site <- which(sites.gate[,1]==s)

                                      # lags measured at this site
                                      lag.site[[s]] <- unique(lag.gate[data.site])

                                      nlags.site[s] <- length(lag.site[[s]])

                                      # average data points from each lag value
                                      acf.site[[s]] <- rep(0+0i,nlags.site[s])
                                      var.site[[s]] <- rep(0,nlags.site[s])
                                      for( l in seq(nlags.site[s])){
                                          lagind <- which(lag.site[[s]][l]==lag.gate[data.site])
                                          for( lind in lagind ){
                                              acf.site[[s]][l] <- acf.site[[s]][l] + acf.gate[data.site][lind]/var.gate[data.site][lind]
                                              var.site[[s]][l] <- var.site[[s]][l] + 1/var.gate[data.site][lind]
                                          }
                                      }
                                      var.site[[s]] <- 1/var.site[[s]]
                                      acf.site[[s]] <- acf.site[[s]] * var.site[[s]]

                                      acf.site[[s]] <- acf.site[[s]] / gainR[s]
                                      var.site[[s]] <- var.site[[s]] / gainR[s]**2

                                      ind.site[[s]] <- rep(s,nlags.site[s])
                                  }

                              }

                              if(sum(nlags.site)>0){

                                  # magnetic field direction
                                  Btmp           <- igrf(date=date[1:3],lat=latitude[h],lon=longitude[h],height=height[h],isv=0,itype=1)
                                  # the model has y-axis to east and z-axis downwards, we have x towards east and z upwards
                                  B[h,]          <- c(Btmp$y,Btmp$x,-Btmp$z)

                                  # parameters from iri model
                                  ptmp           <- iriParams( time=date ,latitude=latitude[h],longitude=longitude[h],heights=height[h])

                                  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
                                  # This is approximately true for all ions, because ion density is much smaller than neutral density
                                  # the densities in outfmsis are in cm^-3
                                  ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,1])['NO+',] )


                                  # initial plasma parameter values, there should not be negative ones, as ion velocity is set to zero
                                  # Also heavier clusters are included in ion mass 30.5, H+ is used as is and the rest is given to O+
                                  parInit        <- pmax( c( ptmp['e-',1] , ptmp['Ti',1] , ptmp['Ti',1], ptmp['Te',1] , ptmp['Te',1] , ioncoll , 0 , 0 , 0 , sum(ptmp[c('O2+','NO+','cluster'),1])/ptmp['e-',1] , 0 , ptmp['H+',1]/ptmp['e-',1] , rep(1,nd) ) , 0 )
                                  parInit[c(10,12)] <- pmax( 0, parInit[c(10,12)])
                                  parInit[c(10,12)] <- pmin( 1, parInit[c(10,12)])
                                  if(sum(parInit[c(10,12)])>1) parInit[c(10,12)] <- parInit[c(10,12)] / sum(parInit[c(10,12)])
                                  parInit[11] <- 1 - sum( parInit[c(10,12)] )
                                  # switched to estimating the molecular ion abundance as 1 - O+ - H+. The above estimation has larger error at low altitudes where heavier ions actually exists
#                                  parInit        <- pmax( c( ptmp['e-',1] , ptmp['Ti',1] , ptmp['Ti',1], ptmp['Te',1] , ptmp['Te',1] , ioncoll , 0 , 0 , 0 , ifelse( (sum(ptmp[c('H+','O+'),1])/ptmp['e-',1]>=1) , 0 , 1-sum(ptmp[c('H+','O+'),1])/ptmp['e-',1]) , ifelse( (ptmp['O+',1]<0) , 0 , ptmp['O+',1]/ptmp['e-',1]) , ifelse( (ptmp['H+',1]<0) , 0 , ptmp['H+',1]/ptmp['e-',1] ) , rep(1,nd) ) , 0 )
                                  parInit[1]     <- max(parInit[1],1e9)

                                  # parameter scaling factors
                                  parScales      <- ISparamScales(parInit,3)

                                  # scale the initial parameter values
                                  initParam      <- scaleParams( parInit , parScales , inverse=F)

                                  # parameter value limits
                                  parLimits      <- ISparamLimits(3,nd)

                                  # scale the parameter limits
                                  limitParam     <- parLimits
                                  limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
                                  limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)

                                  # apriori information
                                  apriori[[h]]   <- aprioriFunction( initParam , nIon=3 , absCalib=absCalib , TiIsotropic=TiIsotropic , TeIsotropic=TeIsotropic , refSite=refsite , siteScales=sScales[,c((h+12),(h+12+nh))] )

                                  model[h,]      <- parInit

                                  fitpar   <- ISparamfit(
                                      acf             = unlist(acf.site),
                                      var             = unlist(var.site),
                                      lags            = unlist(lag.site),
                                      nData           = sum(nlags.site),
                                      fSite           = sites[,2],
                                      aSite           = aSite,
                                      kSite           = kSite,
                                      iSite           = unlist(ind.site),
                                      B               = B[h,],
                                      initParam       = initParam,
                                      invAprioriCovar = apriori[[h]]$invAprioriCovar,
                                      aprioriTheory   = apriori[[h]]$aprioriTheory,
                                      aprioriMeas     = apriori[[h]]$aprioriMeas,
                                      mIon            = c(30.5,16,1),
                                      nIon            = 3,
                                      paramLimits     = limitParam,
                                      directTheory    = ISdirectTheory,
                                      absLimit        = absLimit,
                                      diffLimit       = diffLimit,
                                      scaleFun        = scaleParams,
                                      scale           = parScales,
                                      plotTest        = plotTest,
                                      plotFit         = plotFit,
                                      maxLambda       = maxLambda,
                                      maxIter         = maxIter,
                                      fitFun          = fitFun,
                                      MCMCsettings    = MCMCsettings,
                                      trueHessian     = trueHessian
                                      )

                                  # scale back to physical units
                                  param[h,] <- scaleParams( fitpar$param , scale=parScales , inverse=TRUE )
                                  covar[[h]] <- scaleCovar( fitpar$covar , scale=parScales , inverse=TRUE)
                                  std[h,]   <- sqrt(diag(covar[[h]]))
                                  chisqr[h] <- fitpar[["chisqr"]]
                                  status[h] <- fitpar[["fitStatus"]]
                                  dimnames(covar[[h]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))

                                  contribSites[[h]] <- unique(unlist(ind.site))
                                  if(!is.null(fitpar$MCMC)){
                                      MCMC[[h]] <- fitpar$MCMC
                                      for( sr in seq(dim(MCMC[[h]][["pars"]][1]))) MCMC[[h]][["pars"]][k,] <-  scaleParams( MCMC[[h]][["pars"]][k,] , parScales , inverse=T )
                                      MCMC[[h]][["bestpar"]] <- scaleParams( MCMC[[h]][["bestpar"]] , parScales , inverse=TRUE )
                                      MCMC[[h]][["pars"]] <- t( apply( MCMC[[h]][["pars"]] , FUN=scaleParams , MARGIN=1  , scale=parScales , inverse=TRUE) )
                                      dimnames(MCMC[[h]][["pars"]]) <- list( NULL , c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')) )
                                      names(MCMC[[h]][["bestpar"]]) <- c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep=''))
                                  }
                              }
                          }
                      }

                      time_sec <- iperLimits[k+1]
                      POSIXtime <- as.POSIXlt(time_sec,origin='1970-01-01',tz='ut')
                      std[is.na(std)] <- Inf


                      dimnames(param) <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
                      dimnames(std)   <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
                      dimnames(model)   <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nd),sep='')))
                      names(height) <- paste('gate',seq(nh),sep='')
                      names(latitude) <- paste('gate',seq(nh),sep='')
                      names(longitude) <- paste('gate',seq(nh),sep='')
                      dimnames(B) <- list(paste('gate',seq(nh),sep=''),c('x','y','z'))


                      # save the results to file
                      PP <- list(param=param,std=std,model=model,chisqr=chisqr,status=status,time_sec=time_sec,date=date,POSIXtime=POSIXtime,height=height,latitude=latitude,longitude=longitude,sites=sites,intersect=intersect,covar=covar,B=B,heightLimits.km=hlims/1000,contribSites=contribSites,mIon=c(30.5,16.0,1.0),MCMC=MCMC,timeLimits.s=iperLimits[k:(k+1)],functionCall=functionCall,apriori=apriori)
                      resFile <- file.path( odir , paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "PP.Rdata" , sep=''))
                      save( PP , file=resFile )

                      cat(iperLimits[k+1],'\n')


                  }
              }
              )
          if(nCores>1){
              mcparallel(runcmd)
              if( (k%%nCores == 0) | k==nIper ) mccollect(wait=TRUE)
          }else{
              eval( runcmd )
          }

      }
  }
