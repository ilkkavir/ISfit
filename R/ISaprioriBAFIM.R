ISaprioriBAFIM <- function( PP , date , latitude , longitude , height , nSite ,  nIon , absCalib=FALSE , TiIsotropic=FALSE , TeIsotropic=FALSE, refSite=1 , siteScales=NULL , hTeTi=100 , B=c(0,0,0) , ViPar0=FALSE , nCores=1, BAFIMpar=list(Ne=c(0,Inf,0.05,2.5e11),Ti=c(80,Inf,0.1,30),Te=c(100,Inf,0.1,30),Coll=c(0,0,.1,1),Vipar=c(80,Inf,0.05,2.5),Viperp=c(80,Inf,.05,10),Mp=c(150,500,.05,.01),Op=c(150,500,0.05,0.01),Hp=c(0,0,.05,.01)) , updateFile=TRUE , returnParams=FALSE , ... )
    {
        #
        #
        # Prior model for plasma parameter fits by means of Bayesian Filtering
        #
        # INPUT:
        #  PP              fit output list from the previous integration period (an empty list in the first integration period)
        #  date            measurement time as c(year,month,day,hour,minute,seconds)
        #  latitude        geodetic latitudes of the measurement volumes (deg north)
        #  longitude       geodetic longitudes of the measurement volumes (deg east)
        #  height               heights in km
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
            
        # THIS function has not been completed yet!

        nh <- length(height)

        # IRI parameters (would it be enough to call this with just one lat and lon?
#        IRIlist <- mclapply(seq(nh) , FUN=iriParamsParFun , date=date,latitude=latitude,longitude=longitude,height=height,fitGate=rep(T,nh) , okData=rep(T,nh) , mc.cores=nCores)
        IRIpar <- iriParams( time=date,latitude=mean(latitude),longitude=mean(longitude),heights=height) # IS this accurate enough for low-elevation measurements?


        # parameter value limits
        parLimits      <- ISparamLimits(3,nSite)


        
        apriorilist <- list()

        aprioriIRI <- aprioriBAFIM <- list()
        for(h in seq(nh)){



            ################# IRI parameters ###################################


            
            # parameters from iri model
            ptmp <- IRIpar[,h]
            
            # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
            # This is approximately true for all ions, because ion density is much smaller than neutral density
            ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp)['NO+',] )


            # initial plasma parameter values
            cH <- max(ptmp['H+'],0)
            cO <- max(ptmp['O+'],0)
            cM <- max(sum(ptmp[c('NO+','O2+','cluster')]),0)
            cTot <- cH + cO + cM
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
                          
            parInit <- pmax( c( ptmp['e-'] , ptmp['Ti'] , ptmp['Ti'], ptmp['Te'] , ptmp['Te'] , ioncoll , 0 , 0 , 0 , cM/cTot , cO/cTot , cH/cTot , rep(1,nSite) ) , 0 )
            
            
            parInit[1]     <- max(parInit[1],1e9)

            nPar <- length(parInit)
            
            mIon <- c(30.5,16.0,1)


            
            # parameter scaling factors
            parScales      <- ISparamScales(parInit,3)
            
            # scale the initial parameter values
            aprioriParam      <- scaleParams( parInit , parScales , inverse=F)

           # scale the parameter limits
            limitParam     <- parLimits
            limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
            limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)


            # the apriori covariance matrix will be diagonal, so we begin with
            # a vector of standard deviations, which is easier.
            aprioriStd                   <- vector(mode='numeric',length=nPar)

            # scaled process noise standard deviations. We use these as prior standard deviations in the first integration period
            processStdScale <- scaleParams(c( BAFIMpar$Ne[4] , BAFIMpar$Ti[4] , BAFIMpar$Ti[4] ,BAFIMpar$Te[4] , BAFIMpar$Te[4] , BAFIMpar$Coll[4] , BAFIMpar$Viperp[4] , BAFIMpar$Viperp[4] , BAFIMpar$Vipar[4] , BAFIMpar$Mp[4], BAFIMpar$Op[4], BAFIMpar$Hp[4] ) ,  parScales[1:12],inverse=F)

            
            # The user input standard deviation in the fitted region (when starting the filter), small values elsewhere. 
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

            
#            aprioriIRI[[h]] <- list(aprioriParam=aprioriParam,aprioriStd=aprioriStd,limitParam=limitParam,parScales=parScales,invAprioriCovar=diag(1/aprioriStd**2))
            aprioriIRI[[h]] <- list(aprioriParam=aprioriParam,aprioriStd=aprioriStd,limitParam=limitParam,parScales=parScales,aprioriCovar=diag(aprioriStd**2))

        }







        ############### smooth the plasma parameter profiles in altitude ##########################
        # Plasma scale heights
        kB <- 1.380649e-23
        amu <- 1.66053907e-27
        IRImol <- colSums(IRIpar[c('NO+','O2+','cluster'),])
        IRItot <- IRImol + IRIpar['O+',] + IRIpar['H+',]
        H <- kB * (IRIpar['Ti',] + IRIpar['Te',]) / 2 / ( amu * ( IRIpar['H+',]/IRItot + 16*IRIpar['O+',]/IRItot + 30.5*IRImol/IRItot ) * 9.82 * ( 6372/(6372+height) )**2 )


        
        if (length(PP)>0){

            # time step duration
            dt <-abs( as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - PP$time_sec)

            # scaling factor for the length scales
            hsAlt <- H/1000 * sqrt(dt)

            # height gate widths
            dheights <- diff(PP$heightLimits.km)
            
            # replace unrealistic values and failed fits with the previous predictions
            for(h in seq(nh)){
                okfit <- TRUE
                if(PP$status[h] | PP$chisqr[h]>100 | any(PP$param[h,] < parLimits[1,]) | any(PP$param[h,] > parLimits[2,] )| any(is.na(PP$param[h,])) ){
                    okfit <- FALSE
                }

                
                if(!okfit){
                    Qtmp <- PP$apriori[[h]]$invAprioriCovar
                    Atmp <- PP$apriori[[h]]$aprioriTheory
                    mtmp <- PP$apriori[[h]]$aprioriMeas
                    prec <- t(Atmp)%*%Qtmp%*%Atmp
                    stdp <- sqrt(diag(prec))
                    spre <- outer(stdp,stdp)
                    Stmp <- solve(prec/spre)/spre
                    xtmp <- c(Stmp%*%t(Atmp)%*%Qtmp%*%mtmp)
                    xtmp <- pmin(pmax(xtmp,PP$apriori[[h]]$limitParam[1,]),PP$apriori[[h]]$limitParam[2,])
                    PP$param[h,] <- scaleParams(xtmp,PP$apriori[[h]]$parScales,inverse=T)
                    PP$covar[[h]] <- scaleCovar(Stmp,PP$apriori[[h]]$parScales,inverse=T)
                    PP$std[h,] <- sqrt(diag(PP$covar[[h]]))
                }

            }


            # need at least three gates for the smoothing
            if(nh>2){


                # the smoothing must be done for unscaled parameters, because the scales vary with altitude!

                
                # Form a correlation prior in range (height) direction
                A <- matrix(0,nrow=(nh-1+nh-2),ncol=nh)
                SNe <- STipar <- STiperp <- STepar <- STeperp <- SColl <- SVix <- SViy <- SVipar <-SMp <-  SOp <- SHp <- A[,1]

                Aind <- 1
            
                # The correlation powers solved from known variances, height steps, and correlation lengths
                corrP <- PP$std[,1:12]**2
                
                corrP[,1] <- corrP[,1]*dheights/(BAFIMpar$Ne[3]*hsAlt) # Ne
                corrP[,2] <- corrP[,2]*dheights/(BAFIMpar$Ti[3]*hsAlt) # Tipar
                corrP[,3] <- corrP[,3]*dheights/(BAFIMpar$Ti[3]*hsAlt) # Tiperp
                corrP[,4] <- corrP[,4]*dheights/(BAFIMpar$Te[3]*hsAlt) # Tepar
                corrP[,5] <- corrP[,5]*dheights/(BAFIMpar$Te[3]*hsAlt) # Teperp
                corrP[,6] <- corrP[,6]*dheights/(BAFIMpar$Coll[3]*hsAlt) # Collisions
                corrP[,7] <- corrP[,7]*dheights/(BAFIMpar$Viperp[3]*hsAlt) # Vix
                corrP[,8] <- corrP[,8]*dheights/(BAFIMpar$Viperp[3]*hsAlt) # Viy
                corrP[,9] <- corrP[,9]*dheights/(BAFIMpar$Vipar[3]*hsAlt) # ViB
                corrP[,10] <- corrP[,10]*dheights/(BAFIMpar$Mp[3]*hsAlt) # Molecular ions
                corrP[,11] <- corrP[,11]*dheights/(BAFIMpar$Op[3]*hsAlt) # O+
                corrP[,12] <- corrP[,12]*dheights/(BAFIMpar$Hp[3]*hsAlt) # H+
                
            
                #The first order terms.
                # M is always zero for the first and higher order terms
                # The zeroth-order terms are added later

                for(hind in seq(1,nh-1)){
                    A[Aind,hind]   <- 1
                    A[Aind,hind+1] <- -1
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
            
               # The second order terms
               # NOTE: This is approximately OK also when the altitude resolution changes, because we assume that
               # the parameters are constant within a gate...
                for(hind in seq(2,nh-1)){
                    A[Aind,hind-1] <- 1
                    A[Aind,hind] <- -2
                    A[Aind,hind+1] <- 1

                    SNe[Aind] <- 8 * corrP[hind,1] * ( dheights[hind] / (BAFIMpar$Ne[3] * hsAlt[hind]) )**3
                    STipar[Aind] <- 8 * corrP[hind,2] * ( dheights[hind] / (BAFIMpar$Ti[3] * hsAlt[hind]) )**3
                    STiperp[Aind] <- 8 * corrP[hind,3] * ( dheights[hind] / (BAFIMpar$Ti[3] * hsAlt[hind]) )**3
                    STepar[Aind] <- 8 * corrP[hind,4] * ( dheights[hind] / (BAFIMpar$Te[3] * hsAlt[hind]) )**3
                    STeperp[Aind] <- 8 * corrP[hind,5] * ( dheights[hind] / (BAFIMpar$Te[3] * hsAlt[hind]) )**3
                    SColl[Aind] <- 8 * corrP[hind,6] * ( dheights[hind] / (BAFIMpar$Coll[3] * hsAlt[hind]) )**3
                    SVix[Aind] <- 8 * corrP[hind,7] * ( dheights[hind] / (BAFIMpar$Viperp[3] * hsAlt[hind]) )**3
                    SViy[Aind] <- 8 * corrP[hind,8] * ( dheights[hind] / (BAFIMpar$Viperp[3] * hsAlt[hind]) )**3
                    SVipar[Aind] <- 8 * corrP[hind,9] * ( dheights[hind] / (BAFIMpar$Vipar[3] * hsAlt[hind]) )**3
                    SMp[Aind] <- 8 * corrP[hind,10] * ( dheights[hind] / (BAFIMpar$Mp[3] * hsAlt[hind]) )**3
                    SOp[Aind] <- 8 * corrP[hind,11] * ( dheights[hind] / (BAFIMpar$Op[3] * hsAlt[hind]) )**3
                    SHp[Aind] <- 8 * corrP[hind,12] * ( dheights[hind] / (BAFIMpar$Hp[3] * hsAlt[hind]) )**3
                    
                    Aind <- Aind + 1
                }
                
                # combine all paramters in one large theory matrix
                nn <- dim(A)
                n1 <- nn[1]
                n2 <- nn[2]
                
                # we have 12 parameters
                Acomb <- matrix(0,nrow=12*n1,ncol=12*n2)
                Scomb <- matrix(NaN,nrow=12*n1,ncol=1)
                for(ipar in seq(12)){
                    Acomb[ ((ipar-1)*n1+1) : (ipar*n1) , ((ipar-1)*n2+1) : (ipar*n2) ] <- A
                }

                Scomb <- c( SNe, STipar , STiperp , STepar , STeperp, SColl, SVix , SViy , SVipar , SMp , SOp , SHp )
                Qcomb <- t(Acomb)%*%diag(1/Scomb)%*%Acomb


            }else{
                #zero information if we did not smooth in range
                Qcomb <- matrix(0,ncol=nh*12,nrow=nh*12)
            }




            # The zeroth order terms are measurements and their covariances from the previous step
            Cfit <- matrix(0,ncol=nh*12,nrow=nh*12)
            Mfit <- matrix(NaN,nrow=nh*12,ncol=1)
            for(ih in seq(nh)){
                fitCov <- PP$covar[[ih]]
                Cfit[ ((0:11)*nh + ih) , ((0:11)*nh + ih) ] <- fitCov[1:12,1:12]
                Mfit[ (0:11)*nh + ih ] <- PP$param[ih,1:12]
            }
            

            # zeroth-order precision matrix that contains measurements from all heights
            Cdiagsqrt <- sqrt(diag(Cfit))
            Cscale <- outer(Cdiagsqrt,Cdiagsqrt)
            Qfit <- solve(Cfit/Cscale)/Cscale


            # solve the the whole problem (zeroth, first, and second order terms).
            # Normalize the variances to unit values to stabilise the matrix inversion
            Qsum <- Qfit + Qcomb
            Qdiagsqrt <- sqrt(diag(Qsum))
            Qscale <- outer(Qdiagsqrt,Qdiagsqrt)
            Cpost <- solve(Qsum/Qscale)/Qscale

            Xpost <- Cpost%*%Qfit%*%Mfit


            
            # skip the smoothing if it obviously failed
            if ( any(is.na(Xpost)) | any(is.na(Cpost)) | any(Im(Xpost)!=0) | any(Im(Cpost)!=0) | any(diag(Cpost)<0)){
                Xpost <- Mfit
                Cpost <- Cfit
                print('Error in range smoothing, skipping..')
            }

            if(any(is.na(log10(PP$param[,1])))){
                print(PP$param[,1])
                }
            
        
            #Pick the parameter profiles
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
        
            #Standard deviations. NOTICE: we could pick the full covariance matrices at each height!
            NeErrCorr <- sqrt(diag(Cpost[1:nh,1:nh]));
            TiparErrCorr <- sqrt(diag(Cpost[(nh+1):(2*nh),(nh+1):(2*nh)]));
            TiperpErrCorr <- sqrt(diag(Cpost[(2*nh+1):(3*nh),(2*nh+1):(3*nh)]));
            TeparErrCorr <- sqrt(diag(Cpost[(3*nh+1):(4*nh),(3*nh+1):(4*nh)]));
            TeperpErrCorr <- sqrt(diag(Cpost[(4*nh+1):(5*nh),(4*nh+1):(5*nh)]));
            CollErrCorr <- sqrt(diag(Cpost[(5*nh+1):(6*nh),(5*nh+1):(6*nh)]));
            VixErrCorr <- sqrt(diag(Cpost[(6*nh+1):(7*nh),(6*nh+1):(7*nh)]));
            ViyErrCorr <- sqrt(diag(Cpost[(7*nh+1):(8*nh),(7*nh+1):(8*nh)]));
            ViparErrCorr <- sqrt(diag(Cpost[(8*nh+1):(9*nh),(8*nh+1):(9*nh)]));
            MpErrCorr <- sqrt(diag(Cpost[(9*nh+1):(10*nh),(9*nh+1):(10*nh)]));
            OpErrCorr <- sqrt(diag(Cpost[(10*nh+1):(11*nh),(10*nh+1):(11*nh)]));
            HpErrCorr <- sqrt(diag(Cpost[(11*nh+1):(12*nh),(11*nh+1):(12*nh)]));

#            plot(Xpost[(3*nh+1):(4*nh)],height,xlim=c(0,3000))
#            lines(PP$param[,4],height)
            layout(matrix(seq(12),ncol=4))
            ## plot(log10(PP$param[,1]),height,xlim=c(10,12))
            ## lines(log10(PP$param[,1]+PP$std[,1]),height,col='blue')
            ## nesmooth <- NeCorr
            ## nesmooth[nesmooth<=1] <- 1
            ## lines(log10(nesmooth),height)
            ## lines(log10(nesmooth+NeErrCorr),height,col='red')


            
            plot((PP$param[,1]),height,xlim=c(0,1e12))
            lines((PP$param[,1]+PP$std[,1]),height,col='blue')
            nesmooth <- NeCorr
            nesmooth[nesmooth<=1] <- 1
            lines((nesmooth),height)
            lines((nesmooth+NeErrCorr),height,col='red')

            
            plot(PP$param[,2],height,xlim=c(0,3000))
            lines(TiparCorr,height)
            lines(PP$param[,2]+PP$std[,2],height,col='blue')
            lines(TiparCorr+TiparErrCorr,height,col='red')
            
            plot(PP$param[,3],height,xlim=c(0,3000))
            lines(TiperpCorr,height)
            lines(PP$param[,3]+PP$std[,3],height,col='blue')
            lines(TiperpCorr+TiperpErrCorr,height,col='red')
            
            plot(PP$param[,4],height,xlim=c(0,3000))
            lines(TeparCorr,height)
            lines(PP$param[,4]+PP$std[,4],height,col='blue')
            lines(TeparCorr+TeparErrCorr,height,col='red')

            plot(PP$param[,5],height,xlim=c(0,3000))
            lines(TeperpCorr,height)
            lines(PP$param[,5]+PP$std[,5],height,col='blue')
            lines(TeperpCorr+TeperpErrCorr,height,col='red')
            
            
            plot(PP$param[,6],height,xlim=c(0,1e5))
            lines(CollCorr,height)
            lines(PP$param[,6]+PP$std[,6],height,col='blue')
            lines(CollCorr+CollErrCorr,height,col='red')
            
            plot(PP$param[,7],height,xlim=c(-1,1)*100)
            lines(VixCorr,height)
            lines(PP$param[,7]+PP$std[,7],height,col='blue')
            lines(VixCorr+VixErrCorr,height,col='red')

            plot(PP$param[,8],height,xlim=c(-1,1)*100)
            lines(ViyCorr,height)
            lines(PP$param[,8]+PP$std[,8],height,col='blue')
            lines(ViyCorr+ViyErrCorr,height,col='red')

            plot(PP$param[,9],height,xlim=c(-1,1)*100)
            lines(ViparCorr,height)
            lines(PP$param[,9]+PP$std[,9],height,col='blue')
            lines(ViparCorr+ViparErrCorr,height,col='red')

            plot(PP$param[,10],height,xlim=c(0,1))
            lines(MpCorr,height)
            lines(PP$param[,10]+PP$std[,10],height,col='blue')
            lines(MpCorr+MpErrCorr,height,col='red')

            plot(PP$param[,11],height,xlim=c(0,1))
            lines(OpCorr,height)
            lines(PP$param[,11]+PP$std[,11],height,col='blue')
            lines(OpCorr+OpErrCorr,height,col='red')

            plot(PP$param[,12],height,xlim=c(0,1))
            lines(HpCorr,height)
            lines(PP$param[,12]+PP$std[,12],height,col='blue')
            lines(HpCorr+HpErrCorr,height,col='red')
            
            
#            dt <- abs(as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - PP$time_sec)
            for (hind in seq(nh)){

                # the range-smoothed parameters
                PP$paramRcorr[hind,1:12] <-c(NeCorr[hind],TiparCorr[hind],TiperpCorr[hind],TeparCorr[hind],TeperpCorr[hind],CollCorr[hind],VixCorr[hind],ViyCorr[hind],ViparCorr[hind],MpCorr[hind],OpCorr[hind],HpCorr[hind])

                PP$stdRcorr[hind,1:12] <- c(NeErrCorr[hind],TiparErrCorr[hind],TiperpErrCorr[hind],TeparErrCorr[hind],TeperpErrCorr[hind],CollErrCorr[hind],VixErrCorr[hind],ViyErrCorr[hind],ViparErrCorr[hind],MpErrCorr[hind],OpErrCorr[hind],HpErrCorr[hind])

                PP$covarRcorr[[hind]][1:12,1:12] <- Cpost[ ((0:11)*nh + hind) , ((0:11)*nh + hind) ]

                # the prior model for the next time step
                aprioriBAFIM[[hind]] <- list()
                aprioriBAFIM[[hind]][['aprioriParam']] <- scaleParams( PP$paramRcorr[hind,1:12] , aprioriIRI[[hind]]$parScales[1:12] , inverse=FALSE )

                # process noise standard deviation in normalized units
                processStd <- c( BAFIMpar$Ne[4] , BAFIMpar$Ti[4] , BAFIMpar$Ti[4] ,BAFIMpar$Te[4] , BAFIMpar$Te[4] , BAFIMpar$Coll[4] , BAFIMpar$Viperp[4] , BAFIMpar$Viperp[4] , BAFIMpar$Vipar[4] , BAFIMpar$Mp[4], BAFIMpar$Op[4], BAFIMpar$Hp[4] )*sqrt(dt)

                
                # standard deviations of the smoothed values + the process noise
                aprioriBAFIM[[hind]][['aprioriStd']] <- scaleParams( PP$stdRcorr[hind,1:12] + processStd , aprioriIRI[[hind]]$parScales[1:12] , inverse=FALSE )

                # smoothed plasma parameter error covariance in this gate 
#                aprioriBAFIM[[hind]][['invAprioriCovar']] <- solve(scaleCovar( PP$covarRcorr[[hind]][1:12,1:12] + diag(processStd[1:12])**2, aprioriIRI[[hind]]$parScales[1:12] , inverse=F ))
                aprioriBAFIM[[hind]][['aprioriCovar']] <- scaleCovar( PP$covarRcorr[[hind]][1:12,1:12] + diag(processStd[1:12])**2, aprioriIRI[[hind]]$parScales[1:12] , inverse=F )

                
            }

            Cpred <- Cpost + diag(rep(processStd,each=nh)**2)
                
            stdCpred <- sqrt(diag(Cpred))
            sCpred <- outer(stdCpred,stdCpred)
            Qpred <- solve(Cpred/sCpred)/sCpred

            
            PP$BAFIM_G <- Cfit * t(Cpost * Qfit) * Qpred

        }else{
            # IRI parameters are used in the first iteration step
            aprioriBAFIM <- aprioriIRI
        }


        
        ## Pick the IRI / BAFIM priors according to the limits in BAFIMpar
 
        for(h in seq(nh)){
            

            aprioriParam <- aprioriIRI[[h]]$aprioriParam
            aprioriStd <- aprioriIRI[[h]]$aprioriStd
            limitParam <- aprioriIRI[[h]]$limitParam
            parScales <- aprioriIRI[[h]]$parScales

            fitPar <- rep(FALSE,12)

            # length of the parameter vector
            nPar                         <- length(aprioriParam)

            # number of imaginary apriori "measurements"
            nApriori                     <- ( nPar + 5 )

            # apriori theory matrix
            aprioriTheory                <- matrix( 0 , nrow=nApriori , ncol=nPar )

            # apriori measurement vector
            aprioriMeas                  <- aprioriParam#vector(mode='numeric',length=nApriori)

            # the apriori covariance matrix will be diagonal, so we begin with
            # a vector of standard deviations, which is easier.
            aprioriStd                   <- aprioriStd#vector(mode='numeric',length=nApriori)

            ## invAprioriCovar <- matrix( 0 , nrow=nApriori , ncol=nApriori )
            ## invAprioriCovar[1:nPar,1:nPar] <- aprioriIRI[[h]]$invAprioriCovar
            aprioriCovar <- matrix( 0 , nrow=nApriori , ncol=nApriori )
            aprioriCovar[1:nPar,1:nPar] <- aprioriIRI[[h]]$aprioriCovar

            # apriori parameter values
            aprioriTheory[1:nPar,1:nPar] <- diag(rep(1,nPar))

            # Fill with the smoothed values where appropriate            
            if (height[h]>=BAFIMpar$Ne[1] & height[h]<BAFIMpar$Ne[2]){
                aprioriMeas[1] <- aprioriParam[1] <- aprioriBAFIM[[h]]$aprioriParam[1]
                aprioriStd[1] <- aprioriBAFIM[[h]]$aprioriStd[1]

                fitPar[1] <- TRUE
                # does this really work like this??!??
                # we ar filling in also some kind of cross-information with the parameters we will not fit (and this off-diagonal terms are probably incorrectly weighted!!)
                # we may need to first form the actual covariance and then invert it!!! 
                ## invAprioriCovar[1:12,1] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,1]
                ## invAprioriCovar[1,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[1,1:12]
            }
            
            if (height[h]>=BAFIMpar$Ti[1] & height[h]<BAFIMpar$Ti[2]){
                aprioriMeas[2] <- aprioriParam[2] <- aprioriBAFIM[[h]]$aprioriParam[2]
                aprioriStd[2] <- aprioriBAFIM[[h]]$aprioriStd[2]
                aprioriMeas[3] <- aprioriParam[3] <- aprioriBAFIM[[h]]$aprioriParam[3]
                aprioriStd[3] <- aprioriBAFIM[[h]]$aprioriStd[3]
                fitPar[2:3] <- TRUE
                ## invAprioriCovar[1:12,2] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,2]
                ## invAprioriCovar[2,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[2,1:12]
                ## invAprioriCovar[1:12,3] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,3]
                ## invAprioriCovar[3,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[3,1:12]
            }else{
                # remove IRI model values of Tiperp, these are controlled with the Tipar-Tiperp correlation
                aprioriMeas[3] <- 0
                aprioriTheory[3,] <- 0
            }

            if (height[h]>=BAFIMpar$Te[1] & height[h]<BAFIMpar$Te[2]){
                aprioriMeas[4] <- aprioriParam[4] <- aprioriBAFIM[[h]]$aprioriParam[4]
                aprioriStd[4] <- aprioriBAFIM[[h]]$aprioriStd[4]
                aprioriMeas[5] <- aprioriParam[5] <- aprioriBAFIM[[h]]$aprioriParam[5]
                aprioriStd[5] <- aprioriBAFIM[[h]]$aprioriStd[5]
                fitPar[4:5] <- TRUE
                ## invAprioriCovar[1:12,4] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,4]
                ## invAprioriCovar[4,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[4,1:12]
                ## invAprioriCovar[1:12,5] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,5]
                ## invAprioriCovar[5,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[5,1:12]
           }else{
                # remove IRI model values of Teperp, these are controlled with the Tepar-Teperp correlation
                aprioriMeas[5] <- 0
                aprioriTheory[5,] <- 0
            }

            if (height[h]>=BAFIMpar$Coll[1] & height[h]<BAFIMpar$Coll[2]){
                aprioriMeas[6] <- aprioriParam[6] <- aprioriBAFIM[[h]]$aprioriParam[6]
                aprioriStd[6] <- aprioriBAFIM[[h]]$aprioriStd[6]
                fitPar[6] <- TRUE
                ## invAprioriCovar[1:12,6] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,6]
                ## invAprioriCovar[6,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[6,1:12]
            }

            if (height[h]>=BAFIMpar$Viperp[1] & height[h]<BAFIMpar$Viperp[2]){
                aprioriMeas[7] <- aprioriParam[7] <- aprioriBAFIM[[h]]$aprioriParam[7]
                aprioriStd[7] <- aprioriBAFIM[[h]]$aprioriStd[7]
                aprioriMeas[8] <- aprioriParam[8] <- aprioriBAFIM[[h]]$aprioriParam[8]
                aprioriStd[8] <- aprioriBAFIM[[h]]$aprioriStd[8]
                fitPar[7:8] <- TRUE
                ## invAprioriCovar[1:12,7] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,7]
                ## invAprioriCovar[7,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[7,1:12]
                ## invAprioriCovar[1:12,8] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,8]
                ## invAprioriCovar[8,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[8,1:12]
            }
            
            if (height[h]>=BAFIMpar$Vipar[1] & height[h]<BAFIMpar$Vipar[2]){
                aprioriMeas[9] <- aprioriParam[9] <- aprioriBAFIM[[h]]$aprioriParam[9]
                aprioriStd[9] <- aprioriBAFIM[[h]]$aprioriStd[9]
                fitPar[9] <- TRUE
                ## invAprioriCovar[1:12,9] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,9]
                ## invAprioriCovar[9,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[9,1:12]
            }
            
            if (height[h]>=BAFIMpar$Mp[1] & height[h]<BAFIMpar$Mp[2]){
                aprioriMeas[10] <- aprioriParam[10] <- aprioriBAFIM[[h]]$aprioriParam[10]
                aprioriStd[10] <- aprioriBAFIM[[h]]$aprioriStd[10]
                fitPar[10] <- TRUE
                ## invAprioriCovar[1:12,10] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,10]
                ## invAprioriCovar[10,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[10,1:12]
            }
            
            if (height[h]>=BAFIMpar$Op[1] & height[h]<BAFIMpar$Op[2]){
                aprioriMeas[11] <- aprioriParam[11] <- aprioriBAFIM[[h]]$aprioriParam[11]
                aprioriStd[11] <- aprioriBAFIM[[h]]$aprioriStd[11]
                fitPar[11] <- TRUE
                ## invAprioriCovar[1:12,11] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,11]
                ## invAprioriCovar[11,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[11,1:12]
            }
            
            if (height[h]>=BAFIMpar$Hp[1] & height[h]<BAFIMpar$Hp[2]){
                aprioriMeas[12] <- aprioriParam[12] <- aprioriBAFIM[[h]]$aprioriParam[12]
                aprioriStd[12] <- aprioriBAFIM[[h]]$aprioriStd[12]
                fitPar[12] <- TRUE
                ## invAprioriCovar[1:12,12] <- aprioriBAFIM[[h]]$invAprioriCovar[1:12,12]
                ## invAprioriCovar[12,1:12] <- aprioriBAFIM[[h]]$invAprioriCovar[12,1:12]
            }


            aprioriCovar[1:12,1:12][fitPar,fitPar] <- aprioriBAFIM[[h]][["aprioriCovar"]][1:12,1:12][fitPar,fitPar]
            
            aprioriMeas[1:12] <- aprioriParam[1:12] <- pmax(aprioriMeas[1:12],limitParam[1,1:12])
            aprioriMeas[1:12] <- aprioriParam[1:12] <-  pmin(aprioriMeas[1:12],limitParam[2,1:12])

            ## aprioriMeas[1:nPar]          <- aprioriParam
            
            ## aprioriStd[1]                <- 1e4                                        # electron density
            ## aprioriStd[2]                <- ifelse(height[h]<BAFIMpar$Ti[1],1e-3,2)                       # parallel ion temperature
            ## aprioriStd[3]                <- 2                                          # perpendicular ion temperature
            ## aprioriStd[4]                <- 2                                          # parallel electron temperature
            ## aprioriStd[5]                <- 2                                          # perpendicular electron temperature
            ## aprioriStd[6]                <- ifelse((height[h]>BAFIMpar$Coll[1])&(height[h]<BAFIMpar$Coll[2]),1,1e-3)   # ion-neutral collision frequency
            ## aprioriStd[7]                <- ifelse(height[h]<BAFIMpar$Viperp[1],.1,10)                        # ion velocity, x-component
            ## aprioriStd[8]                <- ifelse(height[h]<BAFIMpar$Viperp[1],.1,10)                        # ion velocity, y-component
            ## aprioriStd[9]                <- ifelse(height[h]<BAFIMpar$Vipar[1],.1,10)                        # ion velocity, z-component
            ## aprioriStd[10:(9+nIon)]      <- 1e-3                                       # ion abundances


            ## # remove model information about perpendicular temperatures
            ## aprioriTheory[3,] <- 0
            ## aprioriTheory[5,] <- 0
            ## aprioriMeas[c(3,5)] <- 0




            ####################### replace the IRI predictions with range-smoothed profiles from the previous fit where appropriate ##


            # replace the IRI value with results from the previous time step (or from the previous prior if the failed)
            # and the standard deviations with that from the previous fit + the process noise
            # 
            # NOTE: this is a development version, in the final version we will smooth the profiles first!
            ## if(length(PP)>0){
            ##     dt <- as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - PP$time_sec
            ##     processStdScale <- scaleParams(c( BAFIMpar$Ne[4] , BAFIMpar$Ti[4] , BAFIMpar$Ti[4] ,BAFIMpar$Te[4] , BAFIMpar$Te[4] , BAFIMpar$Coll[4] , BAFIMpar$Viperp[4] , BAFIMpar$Viperp[4] , BAFIMpar$Vipar[4] , BAFIMpar$Mp[4], BAFIMpar$Op[4], BAFIMpar$Hp[4] ),parScales[1:12],inverse=F)

            ##     if (!PP$status[h]){
            ##         parFit <- scaleParams(PP$param[h,],parScales,inverse=F)
            ##         stdFit <- scaleParams(PP$std[h,],parScales,inverse=F)
            ##     }else{
            ##         parFit <- PP$apriori[[h]]$aprioriParam
            ##         stdFit <- 1/sqrt(diag(PP$apriori[[h]]$invAprioriCovar))
            ##     }
            ##     if (height[h]>=BAFIMpar$Ne[1] & height[h]<BAFIMpar$Ne[2]){
            ##         aprioriMeas[1] <- aprioriParam[1] <- parFit[1]
            ##         aprioriStd[1] <- stdFit[1] + processStdScale[1]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Ti[1] & height[h]<BAFIMpar$Ti[2]){
            ##         aprioriMeas[2] <- aprioriParam[2] <- parFit[2]
            ##         aprioriStd[2] <- stdFit[2] + processStdScale[2]*sqrt(dt)
            ##         aprioriMeas[3] <- aprioriParam[3] <- parFit[3]
            ##         aprioriStd[3] <- stdFit[3] + processStdScale[3]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Te[1] & height[h]<BAFIMpar$Te[2]){
            ##         aprioriMeas[4] <- aprioriParam[4] <- parFit[4]
            ##         aprioriStd[4] <- stdFit[4] + processStdScale[4]*sqrt(dt)
            ##         aprioriMeas[5] <- aprioriParam[5] <- parFit[5]
            ##         aprioriStd[5] <- stdFit[5] + processStdScale[5]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Coll[1] & height[h]<BAFIMpar$Coll[2]){
            ##         aprioriMeas[6] <- aprioriParam[6] <- parFit[6]
            ##         aprioriStd[6] <- stdFit[6] + processStdScale[6]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Viperp[1] & height[h]<BAFIMpar$Viperp[2]){
            ##         aprioriMeas[7] <- aprioriParam[7] <- parFit[7]
            ##         aprioriStd[7] <- stdFit[7] + processStdScale[7]*sqrt(dt)
            ##         aprioriMeas[8] <- aprioriParam[8] <- parFit[8]
            ##         aprioriStd[8] <- stdFit[8] + processStdScale[8]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Vipar[1] & height[h]<BAFIMpar$Vipar[2]){
            ##         aprioriMeas[9] <- aprioriParam[9] <- parFit[9]
            ##         aprioriStd[9] <- stdFit[9] + processStdScale[9]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Mp[1] & height[h]<BAFIMpar$Mp[2]){
            ##         aprioriMeas[10] <- aprioriParam[10] <- parFit[10]
            ##         aprioriStd[10] <- stdFit[10] + processStdScale[10]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Op[1] & height[h]<BAFIMpar$Op[2]){
            ##         aprioriMeas[11] <- aprioriParam[11] <- parFit[11]
            ##         aprioriStd[11] <- stdFit[11] + processStdScale[11]*sqrt(dt)
            ##     }
            ##     if (height[h]>=BAFIMpar$Hp[1] & height[h]<BAFIMpar$Hp[2]){
            ##         aprioriMeas[12] <- aprioriParam[12] <- parFit[12]
            ##         aprioriStd[12] <- stdFit[12] + processStdScale[12]*sqrt(dt)
            ##     }
            ## }





            
            
            if(absCalib){
                aprioriStd[(nIon+10):length(aprioriParam)] <- 1e-3 # fix all sites to the same ACF scale
                diag(aprioriCovar)[(nIon+10):length(aprioriParam)] <- 1e-6
                ## diag(invAprioriCovar)[(nIon+10):length(aprioriParam)] <- 1e6
            }else{
                aprioriStd[(nIon+10):length(aprioriParam)] <- 1   # allow scaling for other sites
                diag(aprioriCovar)[(nIon+10):length(aprioriParam)] <- 1
                ## diag(invAprioriCovar)[(nIon+10):length(aprioriParam)] <- 1
            }
            if(!is.null(siteScales)){
                if(!is.matrix(siteScales)) siteScales <- matrix(siteScales,nrow=1)
                ssinds <- which(!is.na(rowSums(siteScales)))
                aprioriMeas[ssinds+nIon+9] <- siteScales[ssinds,1]  # user-given scaling factors
                if(absCalib){
                    aprioriStd[ssinds+nIon+9] <- siteScales[ssinds,2]
                    diag(aprioriCovar)[ssinds+nIon+9] <- siteScales[ssinds,2]**2
                    ## diag(invAprioriCovar)[ssinds+nIon+9] <- 1/siteScales[ssinds,2]**2
                }
            }
            
            aprioriStd[nIon+9+refSite]     <- 1e-3                 # do not allow scaling at the reference site
            diag(aprioriCovar)[nIon+9+refSite] <- 1e-6
            ## diag(invAprioriCovar)[nIon+9+refSite] <- 1e6
            
            # force certain parameter differences close to zero
            curRow                         <- nPar + 1



            # the temperature ansitropies somewhat diffcult this way,
            # it would perhaps be better to fit the field-aligned temperature
            # and the difference Tperp - Tpar. This will require changes in a number
            # of places but could be worth it...

            # electron temperature anisotropy
            aprioriTheory[curRow,c(4,5)]   <- c(1,-1)
            aprioriMeas[curRow]            <- 0
            if(TeIsotropic){
                aprioriStd[curRow]             <- 1e-3
                diag(aprioriCovar)[curRow] <- 1e-6
                ## diag(invAprioriCovar)[curRow] <- 1e6
            }else{
                aprioriStd[curRow]             <- 1e3
                diag(aprioriCovar)[curRow] <- 1e6
                ## diag(invAprioriCovar)[curRow] <- 1e-6
            }
            curRow                         <- curRow + 1

            # ion temperature anisotropy
            aprioriTheory[curRow,c(2,3)]   <- c(1,-1)
            aprioriMeas[curRow]            <- 0
            if(TiIsotropic){
                aprioriStd[curRow]             <- 1e-3
                diag(aprioriCovar)[curRow] <- 1e-6
                ## diag(invAprioriCovar)[curRow] <- 1e6
            }else{
                aprioriStd[curRow]             <- 1e3
                diag(aprioriCovar)[curRow] <- 1e6
                ## diag(invAprioriCovar)[curRow] <- 1e-6
            }
            curRow                         <- curRow + 1

            # Sum of ion abundances must be one
            aprioriTheory[curRow,10:(nIon+9)] <- 1
            aprioriMeas[curRow] <- 1
            aprioriStd[curRow] <- 1e-3
            diag(aprioriCovar)[curRow] <- 1e-6
            ## diag(invAprioriCovar)[curRow] <- 1e6
            curRow                         <- curRow + 1

            # Te=Ti below hTeTi
            aprioriTheory[curRow,c(2,4)] <- c(1,-1)
            aprioriMeas[curRow] <- 0
            aprioriStd[curRow] <- ifelse(height[h]<hTeTi,1e-3,1e3)
            diag(aprioriCovar)[curRow] <- ifelse(height[h]<hTeTi,1e-6,1e6)
            ## diag(invAprioriCovar)[curRow] <- ifelse(height[h]<hTeTi,1e6,1e-6)
            if(height[h]<hTeTi){
                aprioriTheory[4,] <- 0
                aprioriTheory[5,] <- 0
                aprioriMeas[c(4,5)] <- 0
            }
            curRow                         <- curRow + 1
            
            # optional ViPar=0
            aprioriTheory[curRow,c(7,8,9)] <- B[h,]/sum(sqrt(B[h,]^2))
            aprioriMeas[curRow] <- 0
            aprioriStd[curRow] <- ifelse(ViPar0&all(B[h,]!=0),1e-3,100)
            diag(aprioriCovar)[curRow] <- ifelse(ViPar0&all(B[h,]!=0),1e-6,1e4)
            ##diag(invAprioriCovar)[curRow] <- ifelse(ViPar0&all(B[h,]!=0),1e6,1e-4)

            invAprioriCovar <- solve(aprioriCovar)

#            apriorilist[[h]] <- list(aprioriParam=aprioriParam,aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas,limitParam=limitParam,parScales=parScales,mIon=mIon,nIon,aprioriParamIRI=aprioriIRI[[h]]$aprioriParam,aprioriParamBAFIM=aprioriBAFIM[[h]]$aprioriParam)
            apriorilist[[h]] <- list(aprioriParam=aprioriParam,aprioriTheory=aprioriTheory,invAprioriCovar=invAprioriCovar,aprioriMeas=aprioriMeas,limitParam=limitParam,parScales=parScales,mIon=mIon,nIon,aprioriParamIRI=aprioriIRI[[h]]$aprioriParam,aprioriParamBAFIM=aprioriBAFIM[[h]]$aprioriParam)
        }

        
        if(length(PP)>0){
            
            # copy the range-smoothed data to the default ouputs (but the original ones are also there in the Filter-versions
            PP$param <- PP$paramRcorr
            PP$std <- PP$stdRcorr
            PP$covar <- PP$covarRcorr
        
            # overwrite the output file with the updated copy
            if(updateFile){
                save(PP,file=file.path(PP$resDir,PP$resFile))
            }
        }

        if(returnParams){
            return(PP)
        }else{
            return(apriorilist)
        }

    }
