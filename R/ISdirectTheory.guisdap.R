ISdirectTheory.guisdap <- function( param , scaleFun , nData , mIon , fSite , aSite ,  xSite ,  fAmb , ... ){
#
# Direct theory function for one-dimensional plasma paramter fits. Contains the option for several sites, because
# that allows multi-frequency studies with co-located radars (e.g. EISCAT UHF and VHF systems)
#
# INPUT:
#  param        a vector of plasma parameters
#                 c( Ne , Ti , Te/Ti , nu_en , v , pm2 , ... )

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


  # spectra at each site
  sSite          <- ISspectrum.guisdap(p=sparam,pm0=mIon,fradar=fSite,scattAngle=aSite,freq=xSite)

  # ACFs
  dirtheData     <- sSite%*%fAmb

  return(dirtheData)
  
} # ISdirectTheory
