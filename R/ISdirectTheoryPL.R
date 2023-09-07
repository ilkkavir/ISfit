ISdirectTheoryPL <- function( param , scaleFun , nData , mIon , nSite , iSite , fSite , aSite , kSite ,  xSite ,  B , fAmb , plnSite=plns , pliSite=pliSite, plfSite=plfSite , plaSite=plaSite , plkSite=plkSite , ... ){
#
# Direct theory function for plasma paramter fits with combined ion and plasma line data
#
# INPUT:
#  param        a vector of plasma parameters
#                 c( Ne , Ti , Te/Ti , nu_en , v , pm1 , pm2 , ... )
#   scaleFun    A function that is used for scaling the possibly normalised plasma parameters to physical units
#   nData       Number of data points that will be generated
#   mIon        Ion masses
#   nSite       number of ion line data site
#   iSite       ion line site indices
#   fSite       Transmitter frequency
#   aSite       Scattering angle (difference between incident and scattered wave propagation directions) in degrees
#   kSite       scattering wave vectors
#   xSite       Frequency axis
#   B           magnetic field direction vector
#   fAmb        A matrix of frequency ambiguity functions for all data points. Each column contains one frequency ambiguity
#               function, and they are given in the order in which the data points should be generated. All kinds of scaling
#               factors must be included in the ambiguity functions, this function does not e.g. scale according to frequency
#               step size!
#   plnSite     number of plasma line sites
#   pliSite     plasma line site indices
#   plfSite     transmitter frequencies of plasma line sites
#   plaSite     scattering angles for plasma line sites
#   plkSite     scattering wave vectors for plasma line sites
#
# OUTPUT:
#  dirtheData   A complex vector of direct theory values. length(dirtheData)=nData + 2*plnSite
#
# I. Virtanen 2012
#

  # scale parameters
  sparam         <- scaleFun( param , ... , inverse=T)

  # Conversion to parameter list
  parlist <- ISparamVec2List( sparam , mIon )

  # spectra at each site
  nf <- sapply( xSite , length )
  sSite <- matrix( 0 , nrow=max(nf) , ncol=nSite )
  for(k in seq(nSite)){
      if(nf[k]>0) sSite[1:nf[k],k] <- ISspectrum.3D( ele=parlist[["ele"]] , ion=parlist[["ion"]] , Bdir=B , kdir=kSite[[k]] , fradar=fSite[[k]] , scattAngle=aSite[[k]] , freq=xSite[[k]] )
  }

  # ACFs
  ACF <-  .Call( "crossprods" ,
                as.integer( nData ),
                fAmb,
                as.integer( max(nf) ),
                as.integer( iSite - 1),
                sSite,
                parlist$cSite
                )

# the plasma line frequencies

    if(plnSite > 0){
        for(k in seq(plnSite)){
           # solve f_R +/- from equation (1) of Nicolls et al., 2006 and substitute here!
        }
    }

    
}
