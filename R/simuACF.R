simuACF <- function(ele=c(1e11,300,300,0,0,0,0),ion=list(c(30.5,.7e11,300,300,0,0,0,0),c(16,.3e11,300,300,0,0,0,0),c(1,0,300,300,0,0,0,0)),Bdir=c(0,0,1),kdir=c(0,0,1),fradar=993e6,scattAngle=180,freq=seq(-1000,1000)*100,lags=seq(50)*10e-6){
#
# A simulated ACF for given parameters.
#
# INPUT:
#  ele        c(Ne,Te_par,Te_perp,nu_en,ve_x,ve_y,ve_z)
#  ion        list( c(m1,N1,T1_par,T1_perp,nu_1n,v1_x,v1_y,v1_z) , c(m2,N2,T2_par,T2_perp,nu_2n,v2_x,v2_y,v2_z) , ... )
#  Bdir       c(x,y,z) Magnetic field direction
#  kdir       c(x,y,z) Scattering wave vector direction 
#  fradar     radar frequency in Hz
#  scattAngle scattering angle (the angle between incident and scattered wave vectors, 180 for backscattering)
#  freq       frequency axis
#  lags       time-lags
#
#  ion masses in amu
#  ion densities can be given as absolute values, but they are treated as relative abundances, so that sum(Ni) = Ne
#  only singly-charged ions
#  radar frequency in Hz
#  scattering angle in degrees
#  frequency axis in Hz
#  time-lags in s
#
# OUTPUT:
#  acf        a complex autocorrelation function
#
#
#  I. Virtanen 2012
#

  # IS spectrum
  s <- ISspectrum.3D(ele=ele,ion=ion,kdir=kdir,fradar=fradar,scattAngle=scattAngle,freq=freq,Bdir=Bdir)
  # The ACF via frequency ambiguity functions
  nl <- length(lags)
  acf <- vector(mode='complex',length=nl)

  for(k in seq(nl)) acf[k] <- sum(s * frequencyAmbiguity( lags[k] , freq ) )

  return(acf)
  
} # simuACF
