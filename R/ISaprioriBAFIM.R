ISaprioriBAFIM <- function( PP , date , dateprev , latitude , longitude , height , nSite ,  nIon , absCalib=FALSE , TiIsotropic=FALSE , TeIsotropic=FALSE, refSite=1 , siteScales=NULL , hTeTi=100 , B=c(0,0,0) , ViPar0=FALSE , nCores=1, BAFIMpar=list(Ne=c(0,Inf,0.05,2.5e11),Ti=c(80,Inf,0.1,30),Te=c(100,Inf,0.1,30),Coll=c(0,0,.1,1),Vipar=c(80,Inf,0.05,2.5),Viperp=c(80,Inf,.05,10),Mp=c(150,500,.05,.01),Op=c(150,500,0.05,0.01),Hp=c(0,0,.05,.01),flipchem=c(150,350,.1)) , updateFile=TRUE , returnParams=FALSE , logNe , ... )
    {
        #
        #
        #
        #
        # NOTICE: the flipchem implementation is different from that used in the GUISDAP-BAFIM. Here we use a linear approximation. 
        #
        #
        #
        #
        # Prior model for plasma parameter fits by means of Bayesian Filtering
        #
        # INPUT:
        #  PP              fit output list from the previous integration period (an empty list in the first integration period)
        #  date            measurement time as c(year,month,day,hour,minute,seconds) (end of integration period)
        #  dateprev        measurement time as c(year,month,day,hour,minute,seconds) (start of integration period)
        #  latitude        geodetic latitudes of the measurement volumes (deg north)
        #  longitude       geodetic longitudes of the measurement volumes (deg east)
        #  height          heights in km
        #  nSite           number of receiver sites
        #  nIon            number of ion masses
        #  absCalib        Logical, all site scales are fixed to unity with small variance if absCalib==TRUE and
        #                  siteScales==NULL
        #  TiIsotropic     TRUE if isotropic Ti is assumed, FALSE for bimaxwellian ion velocity distribution
        #  TeIsotropic     TRUE if isotropic Te is assumed, FALSE for bimaxwellian electron velocity distribution
        #  refSite         reference site, whose scale is fixed to unity with small variance
        #  siteScales      a matrix of site scales and their variances or NULL
        #  hTeTi           Te=Ti below hTeTi [km]
        #  B               magnetic field (direction). The default is considered as missing value
        #  ViPar0          logical, force field-aligned ion velocity to zero
        #  nCores          number of cpu cores to use in forks
        #  BAFIMpar        a list of parameters that control the Bayesian filtering.
        #                      list(Ne=c(minh,maxh,length_scale,process_noise), Ti=c(...),Te=c(...),Coll=c(..),Vipar=c(...),Viperp=c(...),
        #                           Mp=c(...),Op=c(...),Hp=c(...),flipchem=c(hmin,hmax,flipchem_std))
        #                      altitudes in km, length scales and project noices as explained in the paper 
        #                      (that has not yet been submitted...)
        #  updateFile      Logical, should the upgraded PP list be written to the output file (default TRUE)
        #  returnParams    Logical, should the upgraded PP list be returned instead of the apriori model (Default FALSE)        
        #
        #  ...             arbitrary parameters to be passed forward to other functions, mainly for compatability reasons
        #
        # OUTPUT:
        #  if returnParams==FALSE (the default), a list with elements
        #    aprioriTheory      apriori theory matrix
        #    aprioriMeas        apriori "measurements"
        #    invAprioriCovar    inverse of apriori covariance matrix
        #  if returnParams==TRUE, an updated PP list with the range-smoothed parameter profiles included
        # 
        #  I. Virtanen 2012, 2013, 2023

        # initialize flipchem with the correct date if the model will be used
        fc <- NULL
        if((BAFIMpar$flipchem[2]-BAFIMpar$flipchem[1]) > 0){
            library(reticulate)
            Sys.setenv(RETICULATE_AUTOCREATE_PACKAGE_VENV="no")
            fcfile <- system.file('python','startflipchem.py',package='ISfit')
            source_python(fcfile)
            idate <- as.integer(date)
            fc <- startflipchem(idate[1],idate[2],idate[3],idate[4],idate[5],idate[6])
        }


        ## PP will be empty in the first time step
        if(length(PP)>0){

            # we will save the PP list again with modified arrays, make copies as necessary
            PP$paramFilter <- PP$param
            PP$stdFilter <- PP$std
            PP$covarFilter <- PP$covar

            # these will be updated with altitude-smoothed profiles:
            PP$paramRcorr <- PP$param
            PP$stdRcorr <- PP$std
            PP$covarRcorr <- PP$covar


        }

        ## number of heights
        nh <- length(height)

        ## IRI parameters. Call with just one lat, lon combination to speed up the analysis.
#        IRIlist <- mclapply(seq(nh) , FUN=iriParamsParFun , date=date,latitude=latitude,longitude=longitude,height=height,fitGate=rep(T,nh) , okData=rep(T,nh) , mc.cores=nCores)
        IRIpar <- iriParams( time=date,latitude=mean(latitude),longitude=mean(longitude),heights=height)

        ## Take a log of IRI Ne if we fit log10(Ne) instead of Ne. 
        if(logNe){
            IRIpar[1,] <- log10(pmax(IRIpar[1,],1e9))
        }

        ## Physically reasonable limits for the plasma parameters
        parLimits      <- ISparamLimits(3,nSite,logNe)

        ## A list for the prior values
        apriorilist <- list()

        ## Lists for prior values from IRI and produced with filtering in time. The appropriate one will be selected at each altiude in the end. 
        aprioriIRI <- aprioriBAFIM <- list()

        ## Form the IRI prior, loop over altitudes 
        for(h in seq(nh)){



            ################# IRI parameters ###################################

            
            ## IRI parameters in this gate
            ptmp <- IRIpar[,h]
            
            ## an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
            ## This is approximately true for all ions, because ion density is much smaller than neutral density
            ## Should be replaced with the more recent formulas by Schunk and Nagy.
            ioncoll        <- sum( ionNeutralCollisionFrequency( c( ifelse(logNe,10^ptmp[1],ptmp[1] ) , ptmp[2:length(ptmp)] ) )['NO+',] )


            ## Initial ion densities. Make sure that these are not negative.
            cH <- max(ptmp['H+'],0)
            cO <- max(ptmp['O+'],0)
            cM <- max(sum(ptmp[c('NO+','O2+','cluster')]),0)
            cTot <- cH + cO + cM
            ## make sure that we have reasonable values in the D region and topside. 
            if(cTot<1e7){
                if(h<150){
                    cH <- 0
                    cO <- 0
                    cM <- 1
                    cTot <- 1
                }else{
                    cM <- 0
                    cTot <- cO + cH
                }
            }

            ## A vector of initial plasma parameter values. 
            parInit <- pmax( c( ptmp['e-'] , ptmp['Ti'] , ptmp['Ti'], ptmp['Te'] , ptmp['Te'] , ioncoll , 0 , 0 , 0 , cM/cTot , cO/cTot , cH/cTot , rep(1,nSite) ) , 0 )
            
            ## Initial Ne must be at least 1e9 m^-3
            parInit[1]     <- max(parInit[1],ifelse(logNe,9,1e9))

            ## number of plasma parameters
            nPar <- length(parInit)

            ## Ion masses (M+, O+, H+)
            mIon <- c(30.5,16.0,1)

            
            ## parameter scaling factors
            parScales      <- ISparamScales(parInit,3,logNe)
            
            ## scale the initial parameter values
            aprioriParam      <- scaleParams( parInit , parScales , inverse=F)

            ## scale the parameter limits
            limitParam     <- parLimits
            limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
            limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)


            ## The apriori covariance matrix will be diagonal, so we begin with
            ## a vector of standard deviations, which is easier.
            aprioriStd                   <- vector(mode='numeric',length=nPar)

            ## Normalized process noise standard deviations. We use these as prior standard deviations in the first integration period
            dt <- abs( as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - as.double(ISOdate(dateprev[1],dateprev[2],dateprev[3],dateprev[4],dateprev[5],dateprev[6])))
            processStdScale <- scaleParams(c( BAFIMpar$Ne[4] , BAFIMpar$Ti[4] , BAFIMpar$Ti[4] ,BAFIMpar$Te[4] , BAFIMpar$Te[4] , BAFIMpar$Coll[4] , BAFIMpar$Viperp[4] , BAFIMpar$Viperp[4] , BAFIMpar$Vipar[4] , BAFIMpar$Mp[4], BAFIMpar$Op[4], BAFIMpar$Hp[4] ) ,  parScales[1:12],inverse=F)*sqrt(dt)

            
            ## The user input standard deviation in the fitted region (when starting the filter), small values elsewhere. 
            aprioriStd[1] <- ifelse( height[h]>=BAFIMpar$Ne[1] & height[h]<BAFIMpar$Ne[2] , processStdScale[1] , 1e-3 )
            aprioriStd[2] <- ifelse( height[h]>=BAFIMpar$Ti[1] & height[h]<BAFIMpar$Ti[2] , processStdScale[2] , 1e-3 )
            aprioriStd[3] <- ifelse( height[h]>=BAFIMpar$Ti[1] & height[h]<BAFIMpar$Ti[2] , processStdScale[3] , 1e-3 )
            aprioriStd[4] <- ifelse( height[h]>=BAFIMpar$Te[1] & height[h]<BAFIMpar$Te[2] , processStdScale[4] , 1e-3 )
            aprioriStd[5] <- ifelse( height[h]>=BAFIMpar$Te[1] & height[h]<BAFIMpar$Te[2] , processStdScale[5] , 1e-3 )
            aprioriStd[6] <- ifelse( height[h]>=BAFIMpar$Coll[1] & height[h]<BAFIMpar$Coll[2] , processStdScale[6] , 1e-3 )
            aprioriStd[7] <- ifelse( height[h]>=BAFIMpar$Viperp[1] & height[h]<BAFIMpar$Viperp[2] , processStdScale[7] , 1e-3 )
            aprioriStd[8] <- ifelse( height[h]>=BAFIMpar$Viperp[1] & height[h]<BAFIMpar$Viperp[2] , processStdScale[8] , 1e-3 )
            aprioriStd[9] <- ifelse( height[h]>=BAFIMpar$Vipar[1] & height[h]<BAFIMpar$Vipar[2] , processStdScale[9] , 1e-3 )
            aprioriStd[10] <- ifelse( height[h]>=BAFIMpar$Mp[1] & height[h]<BAFIMpar$Mp[2] , processStdScale[10] , 1e-3 )
            aprioriStd[11] <- ifelse( height[h]>=BAFIMpar$Op[1] & height[h]<BAFIMpar$Op[2] , processStdScale[11] , 1e-3 )
            aprioriStd[12] <- ifelse( height[h]>=BAFIMpar$Hp[1] & height[h]<BAFIMpar$Hp[2] , processStdScale[12] , 1e-3 )

            
            ## The final IRI prior
            aprioriIRI[[h]] <- list(aprioriParam=aprioriParam,limitParam=limitParam,parScales=parScales,aprioriCovar=diag(aprioriStd**2))

        }







        ############### Smooth the  fitted plasma parameter profiles in altitude ##########################

        ## Plasma scale heights from IRI
        kB <- 1.380649e-23
        amu <- 1.66053907e-27
        IRImol <- colSums(IRIpar[c('NO+','O2+','cluster'),])
        IRItot <- IRImol + IRIpar['O+',] + IRIpar['H+',]
        H <- kB * (IRIpar['Ti',] + IRIpar['Te',]) / 2 / ( amu * ( IRIpar['H+',]/IRItot + 16*IRIpar['O+',]/IRItot + 30.5*IRImol/IRItot ) * 9.82 * ( 6372/(6372+height) )**2 )
        
        
        
        if (length(PP)>0){
            
            ## time step duration
            dt <-abs( as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - PP$time_sec)
            
            if(dt==0){
                dt <- abs( as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - as.double(ISOdate(dateprev[1],dateprev[2],dateprev[3],dateprev[4],dateprev[5],dateprev[6])))
            }
            
            ## scaling factor for the correlation length scales
            hsAlt <- H/1000 * sqrt(dt)

            ## height gate widths
            dheights <- diff(PP$heightLimits.km)

            
            ## replace unrealistic values and failed fits with the previous predictions
            for(h in seq(nh)){
                
                okfit <- TRUE

                ## iri parameters at this altitude
                iriPars  <-  scaleParams(aprioriIRI[[h]]$aprioriParam,aprioriIRI[[h]]$parScales,inverse=T)

                if(PP$status[h] | PP$chisqr[h]>100 | any(PP$param[h,] < parLimits[1,]) | any(PP$param[h,] > parLimits[2,] )| any(is.na(PP$param[h,])) | any(PP$param[h,2:5] < .2*iriPars[2:5])){
                    okfit <- FALSE
                }


                ## if the previous fit failed
                ## neglect the correlations and just copy the prior mean and variance to avoid propagating the correlations
                ## below/above the lowest/highest measured altitude
                if(!okfit){

                    ## number of parameters
                    nParH <- length(PP$apriori[[h]]$aprioriParam)
                    ## the previous prediction
                    PP$param[h,] <- scaleParams(PP$apriori[[h]]$aprioriParam,PP$apriori[[h]]$parScales,inverse=T)
                    ## form the covariance matrix, first scaling factors for matrix normalization to avoid numerical instabilities
                    stdp <- sqrt(diag(PP$apriori[[h]]$invAprioriCovar[1:nParH,1:nParH]))
                    spre <- outer(stdp,stdp)
                    ## invert the precision matrix, scale with spre
                    Stmp <- solve(PP$apriori[[h]]$invAprioriCovar[1:nParH,1:nParH]/spre)/spre
                    PP$covar[[h]] <- scaleCovar(diag(diag(Stmp)),PP$apriori[[h]]$parScales,inverse=T)
                    ## standard deviations from diagonal of the covariance matrix
                    PP$std[h,] <- sqrt(diag(PP$covar[[h]]))

                }
                
            }

            
            ## need at least three gates for the smoothing
            if(nh>2){

                
                ## the smoothing must be done for unnormalized parameters, because the scales vary with altitude!
                
                
                ## Form a correlation prior in range (height) direction
                A <- matrix(0,nrow=(nh-1+nh-2),ncol=nh)
                SNe <- STipar <- STiperp <- STepar <- STeperp <- SColl <- SVix <- SViy <- SVipar <-SMp <-  SOp <- SHp <- A[,1]
                
                Aind <- 1
                
                ## The correlation powers solved from known variances, height steps, and correlation lengths
                corrP <- PP$std[,1:12]**2
                
                corrP[,1]  <- corrP[,1]  * dheights / (BAFIMpar$Ne[3]*hsAlt) # Ne
                corrP[,2]  <- corrP[,2]  * dheights / (BAFIMpar$Ti[3]*hsAlt) # Tipar
                corrP[,3]  <- corrP[,3]  * dheights / (BAFIMpar$Ti[3]*hsAlt) # Tiperp
                corrP[,4]  <- corrP[,4]  * dheights / (BAFIMpar$Te[3]*hsAlt) # Tepar
                corrP[,5]  <- corrP[,5]  * dheights / (BAFIMpar$Te[3]*hsAlt) # Teperp
                corrP[,6]  <- corrP[,6]  * dheights / (BAFIMpar$Coll[3]*hsAlt) # Collisions
                corrP[,7]  <- corrP[,7]  * dheights / (BAFIMpar$Viperp[3]*hsAlt) # Vix
                corrP[,8]  <- corrP[,8]  * dheights / (BAFIMpar$Viperp[3]*hsAlt) # Viy
                corrP[,9]  <- corrP[,9]  * dheights / (BAFIMpar$Vipar[3]*hsAlt) # ViB
                corrP[,10] <- corrP[,10] * dheights / (BAFIMpar$Mp[3]*hsAlt) # Molecular ions
                corrP[,11] <- corrP[,11] * dheights / (BAFIMpar$Op[3]*hsAlt) # O+
                corrP[,12] <- corrP[,12] * dheights / (BAFIMpar$Hp[3]*hsAlt) # H+

                
                
                ## The first order terms.
                ## M is always zero for the first and higher order terms
                ## The zeroth-order terms are added later
                
                for(hind in seq(1,nh-1)){

                    ## here '1' and '-1' are scaled according to the difference between the forward gate centre difference and dheights
                    A[Aind,c(0,1)+hind] <- c(1,-1)/(height[hind+1]-height[hind])*dheights[hind]

                    SNe[Aind]      <-  2 * corrP[hind,1]  * dheights[hind] / (BAFIMpar$Ne[3]*hsAlt[hind])
                    STipar[Aind]   <-  2 * corrP[hind,2]  * dheights[hind] / (BAFIMpar$Ti[3]*hsAlt[hind])
                    STiperp[Aind]  <-  2 * corrP[hind,3]  * dheights[hind] / (BAFIMpar$Ti[3]*hsAlt[hind])
                    STepar[Aind]   <-  2 * corrP[hind,4]  * dheights[hind] / (BAFIMpar$Te[3]*hsAlt[hind])
                    STeperp[Aind]  <-  2 * corrP[hind,5]  * dheights[hind] / (BAFIMpar$Te[3]*hsAlt[hind])
                    SColl[Aind]    <-  2 * corrP[hind,6]  * dheights[hind] / (BAFIMpar$Coll[3]*hsAlt[hind])
                    SVix[Aind]     <-  2 * corrP[hind,7]  * dheights[hind] / (BAFIMpar$Viperp[3]*hsAlt[hind])
                    SViy[Aind]     <-  2 * corrP[hind,8]  * dheights[hind] / (BAFIMpar$Viperp[3]*hsAlt[hind])
                    SVipar[Aind]   <-  2 * corrP[hind,9]  * dheights[hind] / (BAFIMpar$Vipar[3]*hsAlt[hind])
                    SMp[Aind]      <-  2 * corrP[hind,10] * dheights[hind] / (BAFIMpar$Mp[3]*hsAlt[hind])
                    SOp[Aind]      <-  2 * corrP[hind,11] * dheights[hind] / (BAFIMpar$Op[3]*hsAlt[hind])
                    SHp[Aind]      <-  2 * corrP[hind,12] * dheights[hind] / (BAFIMpar$Hp[3]*hsAlt[hind])
                    Aind           <- Aind + 1
                }
            
               #E The second order terms
                for(hind in seq(2,nh-1)){
                    ## '1', '-2' and  '1' are normalized according to the difference between gate centre separations and dheights
                    A[Aind,hind-1] <- 2/((height[hind+1]-height[hind])*(height[hind+1]-height[hind-1]))*dheights[hind]**2
                    A[Aind,hind] <- 2*(height[hind-1]-height[hind+1])/((height[hind+1]-height[hind])*(height[hind]-height[hind-1])*(height[hind+1]-height[hind-1]))*dheights[hind]**2
                    A[Aind,hind+1] <- 2/((height[hind]-height[hind-1])*(height[hind+1]-height[hind-1]))*dheights[hind]**2

                    SNe[Aind]     <- 8 * corrP[hind,1]  * ( dheights[hind] / (BAFIMpar$Ne[3]     * hsAlt[hind]) )**3
                    STipar[Aind]  <- 8 * corrP[hind,2]  * ( dheights[hind] / (BAFIMpar$Ti[3]     * hsAlt[hind]) )**3
                    STiperp[Aind] <- 8 * corrP[hind,3]  * ( dheights[hind] / (BAFIMpar$Ti[3]     * hsAlt[hind]) )**3
                    STepar[Aind]  <- 8 * corrP[hind,4]  * ( dheights[hind] / (BAFIMpar$Te[3]     * hsAlt[hind]) )**3
                    STeperp[Aind] <- 8 * corrP[hind,5]  * ( dheights[hind] / (BAFIMpar$Te[3]     * hsAlt[hind]) )**3
                    SColl[Aind]   <- 8 * corrP[hind,6]  * ( dheights[hind] / (BAFIMpar$Coll[3]   * hsAlt[hind]) )**3
                    SVix[Aind]    <- 8 * corrP[hind,7]  * ( dheights[hind] / (BAFIMpar$Viperp[3] * hsAlt[hind]) )**3
                    SViy[Aind]    <- 8 * corrP[hind,8]  * ( dheights[hind] / (BAFIMpar$Viperp[3] * hsAlt[hind]) )**3
                    SVipar[Aind]  <- 8 * corrP[hind,9]  * ( dheights[hind] / (BAFIMpar$Vipar[3]  * hsAlt[hind]) )**3
                    SMp[Aind]     <- 8 * corrP[hind,10] * ( dheights[hind] / (BAFIMpar$Mp[3]     * hsAlt[hind]) )**3
                    SOp[Aind]     <- 8 * corrP[hind,11] * ( dheights[hind] / (BAFIMpar$Op[3]     * hsAlt[hind]) )**3
                    SHp[Aind]     <- 8 * corrP[hind,12] * ( dheights[hind] / (BAFIMpar$Hp[3]     * hsAlt[hind]) )**3
                    
                    Aind <- Aind + 1
                }
                
                ## Combine all paramters from all gates in one large theory matrix
                nn <- dim(A)
                n1 <- nn[1]
                n2 <- nn[2]
                
                ## We have 12 parameters
                Acomb <- matrix(0,nrow=12*n1,ncol=12*n2)
                Scomb <- matrix(NaN,nrow=12*n1,ncol=1)
                for(ipar in seq(12)){
                    Acomb[ ((ipar-1)*n1+1) : (ipar*n1) , ((ipar-1)*n2+1) : (ipar*n2) ] <- A
                }

                ## form the Fisher information matrix
                Scomb <- c( SNe, STipar , STiperp , STepar , STeperp, SColl, SVix , SViy , SVipar , SMp , SOp , SHp )
                Qcomb <- t(Acomb)%*%diag(1/Scomb)%*%Acomb


            }else{
                #zero information if we did not smooth in range
                Qcomb <- matrix(0,ncol=nh*12,nrow=nh*12)
            }




            ## The zeroth order terms are measurements and their covariances from the previous step
            Cfit <- matrix(0,ncol=nh*12,nrow=nh*12)
            Mfit <- matrix(NaN,nrow=nh*12,ncol=1)
            for(ih in seq(nh)){
                fitCov <- PP$covar[[ih]]
                Cfit[ ((0:11)*nh + ih) , ((0:11)*nh + ih) ] <- fitCov[1:12,1:12]
                Mfit[ (0:11)*nh + ih ] <- PP$param[ih,1:12]
            }
            

            ## zeroth-order precision matrix that contains measurements from all heights
            Cdiagsqrt <- sqrt(diag(Cfit))
            Cscale <- outer(Cdiagsqrt,Cdiagsqrt)
            Qfit <- solve(Cfit/Cscale)/Cscale


            ## solve the the whole problem (zeroth, first, and second order terms).
            ## Normalize the variances to unit values to stabilise the matrix inversion
            Qsum <- Qfit + Qcomb
            Qdiagsqrt <- sqrt(diag(Qsum))
            Qscale <- outer(Qdiagsqrt,Qdiagsqrt)
            Cpost <- solve(Qsum/Qscale)/Qscale

            ## All smoothed profiles in one vector
            Xpost <- Cpost%*%Qfit%*%Mfit


            
            ## skip the smoothing if it obviously failed
            if ( any(is.na(Xpost)) | any(is.na(Cpost)) | any(Im(Xpost)!=0) | any(Im(Cpost)!=0) | any(diag(Cpost)<0)){
                Xpost <- Mfit
                Cpost <- Cfit
                print('Error in range smoothing, skipping..')
            }

        
            ## Pick the smoothed parameter profiles from Xpost
            NeCorr <- Xpost[1:nh];
            TiparCorr <- Xpost[(nh+1):(2*nh)];
            TiperpCorr <- Xpost[(2*nh+1):(3*nh)];
            TeparCorr <- Xpost[(3*nh+1):(4*nh)];
            TeperpCorr <- Xpost[(4*nh+1):(5*nh)];
            CollCorr <- Xpost[(5*nh+1):(6*nh)];
            VixCorr <- Xpost[(6*nh+1):(7*nh)];
            ViyCorr <- Xpost[(7*nh+1):(8*nh)];
            ViparCorr <- Xpost[(8*nh+1):(9*nh)];
            MpCorr <- Xpost[(9*nh+1):(10*nh)];
            OpCorr <- Xpost[(10*nh+1):(11*nh)];
            HpCorr <- Xpost[(11*nh+1):(12*nh)];
        
            #Standard deviations. NOTICE: we will pick also the full covariance matrices at each height later!
            NeErrCorr     <- sqrt(diag(Cpost[              1:nh,             1:nh]));
            TiparErrCorr  <- sqrt(diag(Cpost[(    nh+1): (2*nh),   (nh+1): (2*nh)]));
            TiperpErrCorr <- sqrt(diag(Cpost[(  2*nh+1): (3*nh), (2*nh+1): (3*nh)]));
            TeparErrCorr  <- sqrt(diag(Cpost[(  3*nh+1): (4*nh), (3*nh+1): (4*nh)]));
            TeperpErrCorr <- sqrt(diag(Cpost[(  4*nh+1): (5*nh), (4*nh+1): (5*nh)]));
            CollErrCorr   <- sqrt(diag(Cpost[(  5*nh+1): (6*nh), (5*nh+1): (6*nh)]));
            VixErrCorr    <- sqrt(diag(Cpost[(  6*nh+1): (7*nh), (6*nh+1): (7*nh)]));
            ViyErrCorr    <- sqrt(diag(Cpost[(  7*nh+1): (8*nh), (7*nh+1): (8*nh)]));
            ViparErrCorr  <- sqrt(diag(Cpost[(  8*nh+1): (9*nh), (8*nh+1): (9*nh)]));
            MpErrCorr     <- sqrt(diag(Cpost[( 9*nh+1) :(10*nh), (9*nh+1):(10*nh)]));
            OpErrCorr     <- sqrt(diag(Cpost[(10*nh+1) :(11*nh),(10*nh+1):(11*nh)]));
            HpErrCorr     <- sqrt(diag(Cpost[(11*nh+1) :(12*nh),(11*nh+1):(12*nh)]));



            ## some plots for debugging
            
            if(FALSE){
                
                                        #            plot(Xpost[(3*nh+1):(4*nh)],height,xlim=c(0,3000))
                                        #            lines(PP$param[,4],height)
                layout(matrix(seq(12),ncol=4))
                ## plot(log10(PP$param[,1]),height,xlim=c(10,12))
                ## lines(log10(PP$param[,1]+PP$std[,1]),height,col='blue')
                ## nesmooth <- NeCorr
                ## nesmooth[nesmooth<=1] <- 1
                ## lines(log10(nesmooth),height)
                ## lines(log10(nesmooth+NeErrCorr),height,col='red')
                
                
                
                plot((PP$param[,1]),height,xlim=c(0,ifelse(logNe,12,1e12)))
                lines((PP$param[,1]+PP$std[,1]),height,col='blue')
                nesmooth <- NeCorr
                nesmooth[nesmooth<=1] <- 1
                lines((nesmooth),height)
                lines((nesmooth+NeErrCorr),height,col='red')
                lines((nesmooth+sqrt(NeErrCorr**2+BAFIMpar$Ne[4]**2*dt)),height,col='green')
                
                plot(PP$param[,2],height,xlim=c(0,2000))
                lines(TiparCorr,height)
                lines(PP$param[,2]+PP$std[,2],height,col='blue')
                lines(TiparCorr+TiparErrCorr,height,col='red')
                lines((TiparCorr+sqrt(TiparErrCorr**2+BAFIMpar$Ti[4]**2*dt)),height,col='green')
                
                plot(PP$param[,3],height,xlim=c(0,2000))
                lines(TiperpCorr,height)
                lines(PP$param[,3]+PP$std[,3],height,col='blue')
                lines(TiperpCorr+TiperpErrCorr,height,col='red')
                lines((TiperpCorr+sqrt(TiperpErrCorr**2+BAFIMpar$Ti[4]**2*dt)),height,col='green')
                
                plot(PP$param[,4],height,xlim=c(0,2000))
                lines(TeparCorr,height)
                lines(PP$param[,4]+PP$std[,4],height,col='blue')
                lines(TeparCorr+TeparErrCorr,height,col='red')
                lines((TeparCorr+sqrt(TeparErrCorr**2+BAFIMpar$Te[4]**2*dt)),height,col='green')
                
                plot(PP$param[,5],height,xlim=c(0,2000))
                lines(TeperpCorr,height)
                lines(PP$param[,5]+PP$std[,5],height,col='blue')
                lines(TeperpCorr+TeperpErrCorr,height,col='red')
                lines((TeperpCorr+sqrt(TeperpErrCorr**2+BAFIMpar$Te[4]**2*dt)),height,col='green')
                
                
                plot(PP$param[,6],height,xlim=c(0,1e4))
                lines(CollCorr,height)
                lines(PP$param[,6]+PP$std[,6],height,col='blue')
                lines(CollCorr+CollErrCorr,height,col='red')
                lines((CollCorr+sqrt(CollErrCorr**2+BAFIMpar$Coll[4]**2*dt)),height,col='green')
                
                plot(PP$param[,7],height,xlim=c(-1,1)*100)
                lines(VixCorr,height)
                lines(PP$param[,7]+PP$std[,7],height,col='blue')
                lines(VixCorr+VixErrCorr,height,col='red')
                lines((VixCorr+sqrt(VixErrCorr**2+BAFIMpar$Viperp[4]**2*dt)),height,col='green')
                
                plot(PP$param[,8],height,xlim=c(-1,1)*100)
                lines(ViyCorr,height)
                lines(PP$param[,8]+PP$std[,8],height,col='blue')
                lines(ViyCorr+ViyErrCorr,height,col='red')
                lines((ViyCorr+sqrt(ViyErrCorr**2+BAFIMpar$Viperp[4]**2*dt)),height,col='green')
                
                plot(PP$param[,9],height,xlim=c(-1,1)*100)
                lines(ViparCorr,height)
                lines(PP$param[,9]+PP$std[,9],height,col='blue')
                lines(ViparCorr+ViparErrCorr,height,col='red')
                lines((ViparCorr+sqrt(ViparErrCorr**2+BAFIMpar$Vipar[4]**2*dt)),height,col='green')
                
                plot(PP$param[,10],height,xlim=c(0,1))
                lines(MpCorr,height)
                lines(PP$param[,10]+PP$std[,10],height,col='blue')
                lines(MpCorr+MpErrCorr,height,col='red')
                lines((MpCorr+sqrt(MpErrCorr**2+BAFIMpar$Mp[4]**2*dt)),height,col='green')
                
                plot(PP$param[,11],height,xlim=c(0,1))
                lines(OpCorr,height)
                lines(PP$param[,11]+PP$std[,11],height,col='blue')
                lines(OpCorr+OpErrCorr,height,col='red')
                lines((OpCorr+sqrt(OpErrCorr**2+BAFIMpar$Op[4]**2*dt)),height,col='green')
                
                plot(PP$param[,12],height,xlim=c(0,1))
                lines(HpCorr,height)
                lines(PP$param[,12]+PP$std[,12],height,col='blue')
                lines(HpCorr+HpErrCorr,height,col='red')
                lines((HpCorr+sqrt(HpErrCorr**2+BAFIMpar$Hp[4]**2*dt)),height,col='green')
                
                mtext(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6]),side=3,line=-2,outer=T)
            }
            ## end of the debugging plots

            
            ## collect the smoothed parameters and their error covariances at each height
            for (hind in seq(nh)){
                
                ## the range-smoothed parameters
                PP$paramRcorr[hind,1:12] <-c(NeCorr[hind],TiparCorr[hind],TiperpCorr[hind],TeparCorr[hind],TeperpCorr[hind],CollCorr[hind],VixCorr[hind],ViyCorr[hind],ViparCorr[hind],MpCorr[hind],OpCorr[hind],HpCorr[hind])

                ## standard deviations
                PP$stdRcorr[hind,1:12] <- c(NeErrCorr[hind],TiparErrCorr[hind],TiperpErrCorr[hind],TeparErrCorr[hind],TeperpErrCorr[hind],CollErrCorr[hind],VixErrCorr[hind],ViyErrCorr[hind],ViparErrCorr[hind],MpErrCorr[hind],OpErrCorr[hind],HpErrCorr[hind])

                ## error covariance matrices
                PP$covarRcorr[[hind]][1:12,1:12] <- Cpost[ ((0:11)*nh + hind) , ((0:11)*nh + hind) ]

                ## the prior model for the next time step
                aprioriBAFIM[[hind]] <- list()
                aprioriBAFIM[[hind]][['aprioriParam']] <- scaleParams( PP$paramRcorr[hind,1:12] , aprioriIRI[[hind]]$parScales[1:12] , inverse=FALSE )
                
                ## process noise standard deviation in normalized units
                processStd <- c( BAFIMpar$Ne[4] , BAFIMpar$Ti[4] , BAFIMpar$Ti[4] ,BAFIMpar$Te[4] , BAFIMpar$Te[4] , BAFIMpar$Coll[4] , BAFIMpar$Viperp[4] , BAFIMpar$Viperp[4] , BAFIMpar$Vipar[4] , BAFIMpar$Mp[4], BAFIMpar$Op[4], BAFIMpar$Hp[4] )*sqrt(dt)
                
                
                ## Error covariance matrix of the prediction (covariance of the smoothed parameters + process noise)
                aprioriBAFIM[[hind]][['aprioriCovar']] <- scaleCovar( PP$covarRcorr[[hind]][1:12,1:12] + diag(processStd[1:12])**2, aprioriIRI[[hind]]$parScales[1:12] , inverse=F )


            }

            
        }else{
            ## IRI parameters are used in the first iteration step
            aprioriBAFIM <- aprioriIRI
        }


        
        ## Pick the IRI / BAFIM priors according to the limits in BAFIMpar
 
        for(h in seq(nh)){

            ## first copy the IRI values
            aprioriParam <- aprioriIRI[[h]]$aprioriParam
            limitParam <- aprioriIRI[[h]]$limitParam
            parScales <- aprioriIRI[[h]]$parScales


            fitPar <- rep(FALSE,12)

            ## length of the parameter vector
            nPar <- length(aprioriParam)

            ## check if flipchem will be used in this gate
            flipchemfit <- FALSE
            if(height[h]>=BAFIMpar$flipchem[1] & height[h]<=BAFIMpar$flipchem[2]){
                flipchemfit <- TRUE
            }
            
            ## number of imaginary apriori "measurements"
            nApriori <- ifelse( flipchemfit , nPar + 7 ,  nPar + 6 )

            ## apriori theory matrix
            aprioriTheory <- matrix( 0 , nrow=nApriori , ncol=nPar )

            ## apriori measurement vector
            aprioriMeas <- aprioriParam

            ## first copy the IRI covariances
            aprioriCovar <- matrix( 0 , nrow=nApriori , ncol=nApriori )
            aprioriCovar[1:nPar,1:nPar] <- aprioriIRI[[h]]$aprioriCovar

            ## apriori parameter values
            aprioriTheory[1:nPar,1:nPar] <- diag(rep(1,nPar))

            ## Fill with the smoothed values where appropriate

            ## Ne
            if (height[h]>=BAFIMpar$Ne[1] & height[h]<BAFIMpar$Ne[2]){
                aprioriMeas[1] <- aprioriParam[1] <- aprioriBAFIM[[h]]$aprioriParam[1]
                fitPar[1] <- TRUE
            }

            ## Ti
            if (height[h]>=BAFIMpar$Ti[1] & height[h]<BAFIMpar$Ti[2]){
                aprioriMeas[2] <- aprioriParam[2] <- aprioriBAFIM[[h]]$aprioriParam[2]
                aprioriMeas[3] <- aprioriParam[3] <- aprioriBAFIM[[h]]$aprioriParam[3]
                fitPar[2:3] <- TRUE
            }else{
                ## remove IRI model values of Tiperp, these are controlled with the Tipar-Tiperp correlation
                aprioriMeas[3] <- 0
                aprioriTheory[3,] <- 0
            }

            ## Te
            if (height[h]>=BAFIMpar$Te[1] & height[h]<BAFIMpar$Te[2]){
                aprioriMeas[4] <- aprioriParam[4] <- aprioriBAFIM[[h]]$aprioriParam[4]
                aprioriMeas[5] <- aprioriParam[5] <- aprioriBAFIM[[h]]$aprioriParam[5]
                fitPar[4:5] <- TRUE
           }else{
                ## remove IRI model values of Teperp, these are controlled with the Tepar-Teperp correlation
                aprioriMeas[5] <- 0
                aprioriTheory[5,] <- 0
            }

            ## Collisions
            if (height[h]>=BAFIMpar$Coll[1] & height[h]<BAFIMpar$Coll[2]){
                aprioriMeas[6] <- aprioriParam[6] <- aprioriBAFIM[[h]]$aprioriParam[6]
                fitPar[6] <- TRUE
            }

            ## Vi perpendicular components
            if (height[h]>=BAFIMpar$Viperp[1] & height[h]<BAFIMpar$Viperp[2]){
                aprioriMeas[7] <- aprioriParam[7] <- aprioriBAFIM[[h]]$aprioriParam[7]
                aprioriMeas[8] <- aprioriParam[8] <- aprioriBAFIM[[h]]$aprioriParam[8]
                fitPar[7:8] <- TRUE
            }

            ## Vi parallel
            if (height[h]>=BAFIMpar$Vipar[1] & height[h]<BAFIMpar$Vipar[2]){
                aprioriMeas[9] <- aprioriParam[9] <- aprioriBAFIM[[h]]$aprioriParam[9]
                fitPar[9] <- TRUE
            }

            ## Molecular ions
            if (height[h]>=BAFIMpar$Mp[1] & height[h]<BAFIMpar$Mp[2]){
                aprioriMeas[10] <- aprioriParam[10] <- aprioriBAFIM[[h]]$aprioriParam[10]
                fitPar[10] <- TRUE
            }

            ## O+ ions
            if (height[h]>=BAFIMpar$Op[1] & height[h]<BAFIMpar$Op[2]){
                aprioriMeas[11] <- aprioriParam[11] <- aprioriBAFIM[[h]]$aprioriParam[11]
                fitPar[11] <- TRUE
            }
            
            ## H+ ions
            if (height[h]>=BAFIMpar$Hp[1] & height[h]<BAFIMpar$Hp[2]){
                aprioriMeas[12] <- aprioriParam[12] <- aprioriBAFIM[[h]]$aprioriParam[12]
                fitPar[12] <- TRUE
            }


            ## replace the IRI covariances with the prediction where appropriate
            aprioriCovar[1:12,1:12][fitPar,fitPar] <- aprioriBAFIM[[h]][["aprioriCovar"]][1:12,1:12][fitPar,fitPar]
            
            ## Force the prior values to be within the physically reasonable limits
            aprioriMeas[1:12] <- aprioriParam[1:12] <- pmax(aprioriMeas[1:12],limitParam[1,1:12])
            aprioriMeas[1:12] <- aprioriParam[1:12] <-  pmin(aprioriMeas[1:12],limitParam[2,1:12])

            
            ## The scaling factors should not be fitted if we have absolute calibration
            if(absCalib){
                diag(aprioriCovar)[(nIon+10):length(aprioriParam)] <- 1e-6
            }else{
                diag(aprioriCovar)[(nIon+10):length(aprioriParam)] <- 1
            }

            ## scaling factors from a calibration measurement
            if(!is.null(siteScales)){
                if(!is.matrix(siteScales)) siteScales <- matrix(siteScales,nrow=1)
                ssinds <- which(!is.na(rowSums(siteScales)))
                aprioriMeas[ssinds+nIon+9] <- siteScales[ssinds,1]  # user-given scaling factors
                if(absCalib){
                    diag(aprioriCovar)[ssinds+nIon+9] <- siteScales[ssinds,2]**2
                }
            }
            
            ## we must have one absolutely calibrated reference site
            diag(aprioriCovar)[nIon+9+refSite] <- 1e-6
            
            ## force certain parameter differences close to zero
            curRow                         <- nPar + 1



            ## the temperature ansitropies somewhat diffcult this way,
            ## it would perhaps be better to fit the field-aligned temperature
            ## and the difference Tperp - Tpar. This will require changes in a number
            ## of places but could be worth it...

            ## electron temperature anisotropy
            aprioriTheory[curRow,c(4,5)]   <- c(1,-1)
            aprioriMeas[curRow]            <- 0
            if(TeIsotropic){
                diag(aprioriCovar)[curRow] <- 1e-6
            }else{
                diag(aprioriCovar)[curRow] <- 1e6
            }
            curRow                         <- curRow + 1

            ## ion temperature anisotropy
            aprioriTheory[curRow,c(2,3)]   <- c(1,-1)
            aprioriMeas[curRow]            <- 0
            if(TiIsotropic){
                diag(aprioriCovar)[curRow] <- 1e-6
            }else{
                diag(aprioriCovar)[curRow] <- 1e6
            }
            curRow                         <- curRow + 1

            ## Sum of ion abundances must be one
            aprioriTheory[curRow,10:(nIon+9)] <- 1
            aprioriMeas[curRow] <- 1
            diag(aprioriCovar)[curRow] <- 1e-6
            curRow                         <- curRow + 1

            ## Te=Ti below hTeTi. Ne cannot be high when Ti>Te, either
            TeTiForce <- FALSE
            if(length(PP)>0){
                if(PP$param[h,1]>5e11 & PP$param[h,4]>PP$param[h,2]*1.05){
                    TeTiForce <- TRUE
                }
            }
            aprioriTheory[curRow,c(2,4)] <- c(1,-1)
            aprioriMeas[curRow] <- 0
            if(TeTiForce){
                diag(aprioriCovar)[curRow] <- ifelse(height[h]<hTeTi,1e-6,.01)
            }else{
                diag(aprioriCovar)[curRow] <- ifelse(height[h]<hTeTi,1e-6,1e6)
            }
            if(height[h]<hTeTi){
                aprioriTheory[4,] <- 0
                aprioriTheory[5,] <- 0
                aprioriMeas[c(4,5)] <- 0
            }
            curRow                         <- curRow + 1
            aprioriTheory[curRow,c(3,5)] <- c(1,-1)
            aprioriMeas[curRow] <- 0
            if(TeTiForce){
                diag(aprioriCovar)[curRow] <- ifelse(height[h]<hTeTi,1e-6,.01)
            }else{
                diag(aprioriCovar)[curRow] <- ifelse(height[h]<hTeTi,1e-6,1e6)
            }
            curRow                         <- curRow + 1
            
            ## optional ViPar=0
            aprioriTheory[curRow,c(7,8,9)] <- B[h,]/sum(sqrt(B[h,]^2))
            aprioriMeas[curRow] <- 0
            diag(aprioriCovar)[curRow] <- ifelse(ViPar0&all(B[h,]!=0),1e-6,1e4)
            curRow                         <- curRow + 1


            ## Optional flipchem input (linear approximation of the chemistry model)
            if (flipchemfit){
                fcApriori <- aprioriFlipchem( param=aprioriParam , flipchem=fc , flipchemStd=BAFIMpar$flipchem[3],  lat=latitude[h] , lon=longitude[h] , h=height[h] , scaleFun=scaleParams , scale=parScales , logNe , ... )
                aprioriTheory[curRow,] <- fcApriori$A
                aprioriMeas[curRow] <- fcApriori$m
                aprioriCovar[curRow,curRow] <- fcApriori$var
            }else{
                aprioriUpdateFunction <- NULL
            }

            ## Fisher information matrix for the iterative solver
            invAprioriCovar <- solve(aprioriCovar)

            ## The final prior model list for the iterative solver
            apriorilist[[h]] <- list(aprioriParam=aprioriParam,aprioriTheory=aprioriTheory,invAprioriCovar=invAprioriCovar,aprioriMeas=aprioriMeas,limitParam=limitParam,parScales=parScales,mIon=mIon,nIon=nIon,aprioriParamIRI=aprioriIRI[[h]]$aprioriParam,aprioriParamBAFIM=aprioriBAFIM[[h]]$aprioriParam,flipchem=fc,flipchemStd=BAFIMpar$flipchem[3],aprioriUpdateFunction=aprioriUpdateFunction)
        }

        ## If there was a previous fit
        if(length(PP)>0){
            
            ## Copy the range-smoothed data to the default ouputs (but the original ones are also there in the Filter-versions
            PP$param <- PP$paramRcorr
            PP$std <- PP$stdRcorr
            PP$covar <- PP$covarRcorr
            PP$BAFIMpar <- BAFIMpar
        
            # Overwrite the output file with the updated copy
            if(updateFile){
                save(PP,file=file.path(PP$resDir,PP$resFile))
            }
        }

        ## return either the plasma parameter list or the prior model, depending on the input argument
        if(returnParams){
            return(PP)
        }else{
            return(apriorilist)
        }

    }
