writeVascii <- function( PPlist , fname='Vi.dat' ){
    #
    # Write ion velocity components to ascii files
    #
    #
    #
    #
    # IV 2016
    #

    fcon <- file( fname , 'w' )

    dims <- dim(PPlist[["param"]])
    nh <- dims[1]
    nt <- dims[3]

    cat('#
# Ion velocity estimates from multistatic incoherent scatter analysis.
#
# The columns are:
#
#    1. unix time at end of integration period (seconds since 1970-01-01)
#    2. Height (km)
#  3-8. Ion velocity components (m/s) toward
#      3. Geographic East
#      4. Geographic North
#      5. Zenith (positive upwards)
#      6. Geomagnetic Eeast
#      7. Geomagnetic North (perpendicular to B, i.e. not horizontal)
#      8. Along Magetic field (positive downward)
#  9-14. Covariances of the velocity estimates in geographic coordinates (m^2/s^2)
#      9. Variance East
#     10. Variance North
#     11. Variance Up
#     12. Covariance East-North
#     13. Covariance East-Up
#     14. Covariance North-Up
# 15-17. Same as 9-11, but for the geomagnetic coordinates. 
#
#
# IV 2016
# 
#    Timestamp   Height   East  North     Up  Beast Bnorth   Bpar    Var(East)   Var(North)      Var(Up)     Cov(E,N)    Cov(E-Up)    Cov(N-Up)   Var(Beast)  Var(Bnorth)    Var(Bpar)\n',file=fcon,append=FALSE)

    for( t in seq(nt) ){
        for( h in seq(nh) ){
            cat(sprintf("%14.3f %8.3f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f\n",as.numeric(PPlist[["POSIXtime"]][[t]]),PPlist[["height"]][h,t],PPlist[["param"]][h,'Vix',t],PPlist[["param"]][h,'Viy',t],PPlist[["param"]][h,'Viz',t],PPlist[["param"]][h,'ViBx',t],PPlist[["param"]][h,'ViBy',t],PPlist[["param"]][h,'ViB',t],PPlist[["covar"]][h,'Vix','Vix',t],PPlist[["covar"]][h,'Viy','Viy',t],PPlist[["covar"]][h,'Viz','Viz',t],PPlist[["covar"]][h,'Vix','Viy',t],PPlist[["covar"]][h,'Vix','Viz',t],PPlist[["covar"]][h,'Viy','Viz',t],PPlist[["std"]][h,'ViBx',t]^2,PPlist[["std"]][h,'ViBy',t]^2,PPlist[["std"]][h,'ViB',t]^2),file=fcon,append=TRUE)
        }
    }
    
    close(fcon)
}
