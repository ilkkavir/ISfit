ISdirectTheory.general <- function( param , scaleFun , nData , nIon , nSite , iSite , fSite , aSite , kSite , xSite , B , fAmb , ... ){
#
# Direct theory function for plasma paramter fits.
#
# INPUT:
#  param        a vector of plasma parameters and ACF scales:
#                 c( Ne , Te_par , Te_perp , nu_en , ve_x , ve_y , ve_z ,
#               m1 , N1 , T1_par , T1_perp , nu_1n , v1_x , v1_y , v1_z ,
#               m2 , N2 , T2_par , ...
#               ...
#               mnIon , ...        , vnIon_z ,
#               s_1 , s_2 , ... , s_Nsite)
#   scaleFun    A function that is used for scaling the possibly normalised plasma parameters to physical units
#   nData       Number of data points that will be generated
#   nIon        Number of ion masses
#   nSite       Number of receiver sites
#   iSite       Site indices, an integer vector giving the site number for each data point. length(iSite)=nData
#   fSite       Transmitter frequency at each site. Including this enables combined UHF&VHF analysis. length(fSite)=nSite
#   aSite       Scattering angle (difference between incident and scattered wave propagation directions) in degrees. length(aSite)=nSite
#   kSite       A list of scattering wave vectors at each site. In the same cartesian system that is used for the velocities.
#               length(kSite)=nSite , length(kSite[[i]])=3
#   xSite       A list of frequency axes for spectrum calculation at each site. length(xSite)=nSite
#   B           Magnetic field direction vector, common for all sites because the same volume is being observed. length(B)=3
#   fAmb        A matrix of frequency ambiguity functions for all data points. Each column contains one frequency ambiguity
#               function, and they are given in the order in which the data points should be generated. If different frequency
#               axis lengths are used at different sites, the shorter ones are zero-padded. All kinds of scaling factors must
#               be included in the ambiguity functions, this function does not e.g. scale according to frequency step size!
#
# OUTPUT:
#  dirtheData   A complex vector of direct theory values. length(dirtheData)=nData
#
# I. Virtanen 2012
#

  # scale parameters
  sparam         <- scaleFun( param , ... , inverse=T)

  # conversion to parameter lists
  parlist        <- ISparamVec2List.general(sparam,nIon)

  # spectra at each site
  nf             <- sapply( xSite , length )
  sSite          <- matrix(0,nrow=max(nf),ncol=nSite)
  for(k in seq(nSite)){
    sSite[1:nf[k],k]     <- ISspectrum.3D( ele=parlist$ele , ion=parlist$ion , Bdir=B , kdir=kSite[[k]] , fradar=fSite[[k]] , scattAngle=aSite[[k]] , freq=xSite[[k]] )
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

} # ISdirectTheory.general
