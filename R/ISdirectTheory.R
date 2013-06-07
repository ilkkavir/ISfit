ISdirectTheory <- function( param , scaleFun , nData , mIon , nSite , iSite , fSite , aSite , kSite ,  xSite ,  B , fAmb , ... ){
#
# Direct theory function for one-dimensional plasma paramter fits. Contains the option for several sites, because
# that allows multi-frequency studies with co-located radars (e.g. EISCAT UHF and VHF systems)
#
# INPUT:
#  param        a vector of plasma parameters
#                 c( Ne , Ti , Te/Ti , nu_en , v , pm1 , pm2 , ... )
#   scaleFun    A function that is used for scaling the possibly normalised plasma parameters to physical units
#   nData       Number of data points that will be generated
#   mIon        Ion masses
#   fSite       Transmitter frequency
#   aSite       Scattering angle (difference between incident and scattered wave propagation directions) in degrees
#   xSite       Frequency axis
#   fAmb        A matrix of frequency ambiguity functions for all data points. Each column contains one frequency ambiguity
#               function, and they are given in the order in which the data points should be generated. All kinds of scaling
#               factors must be included in the ambiguity functions, this function does not e.g. scale according to frequency
#               step size!
#
# OUTPUT:
#  dirtheData   A complex vector of direct theory values. length(dirtheData)=nData
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
      sSite[1:nf[k],k] <- ISspectrum.3D( ele=parlist[["ele"]] , ion=parlist[["ion"]] , Bdir=B , kdir=kSite[[k]] , fradar=fSite[[k]] , scattAngle=aSite[[k]] , freq=xSite[[k]] )
  }

  # ACFs
  return( .Call( "crossprods" ,
                as.integer( nData ),
                fAmb,
                as.integer( max(nf) ),
                as.integer( iSite - 1),
                sSite,
                parlist$cSite
                )
         )
  
} # ISdirectTheory
