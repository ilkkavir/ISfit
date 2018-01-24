writeVNeasciiVKSKAIRA <- function( PPlist , fname='Vi.dat' ){
    #
    # Write ion velocity components and electron densities to ascii files
    #
    # This version is for 4-static VHF + KIR + SOD + KAIRA measurements,
    # pick the vertical velocity from
    # projections along the VHF beam, and the horizontal component towards KAIRA
    #
    # The VHF must be at first row of PPlist$sites, KIR at the second row, and SOD at the third row
    #
    #
    # IV 2016, 2017
    #

    fcon <- file( fname , 'w' )

    dims <- dim(PPlist[["param"]])
    nh <- dims[1]
    nt <- dims[3]

    cat('#
# Ion velocity estimates from multistatic incoherent scatter analysis.
# This is a special version for bistatic observations with VHF and KAIRA
#
# The columns are:
#
#    1. unix time at end of integration period (seconds since 1970-01-01)
#    2. Height (km)
#  3-4. Ion velocity components (m/s)
#      3. Zenith (positive upwards)
#      4. Horizontal component toward KAIRA
#  5-6. Variances of the velocity estimates  (m^2/s^2)
#      5. Variance zenith
#      6. Variance horizontal toward KAIRA
#  7-8. Electron density and its variance
#      7. Electron density (m^-3)
#      8. Variance of electron density (m^-6)
#
#
# IV 2016, 2017
#
#    Timestamp   Height     Up hKAIRA      Var(Up)  Var(hKAIRA)      Ne Var(Ne)\n',file=fcon,append=FALSE)


    # list of dimension names for horizontal components toward kaira.
    nkvhor <- c('ViR4hor','ViR5hor','ViR6hor','ViR7hor','ViR8hor','ViR9hor','ViR10hor','ViR11hor','ViR12hor')

    for( t in seq(nt) ){
        for( h in seq(nh) ){
            # try to find the correct horizontal component
            vkhor <- NA
            vkhorstd <- NA
            for (bn in nkvhor){
                if(!is.na(PPlist[["param"]][h,bn,t])){
                    vkhor <- PPlist[["param"]][h,bn,t]
                    vkhorstd <- PPlist[["std"]][h,bn,t]
                    break
                }
            }
            # ViR1 = -ViUp for vertical VHF
            cat(sprintf("%14.3f %8.3f %6.0f %6.0f %12.0f %12.0f %3.2e %3.2e\n",as.numeric(PPlist[["POSIXtime"]][[t]]),PPlist[["height"]][h,t],-PPlist[["param"]][h,'ViR1',t],vkhor,PPlist[["std"]][h,'ViR1',t]**2,vkhorstd**2,PPlist[["param"]][h,'Ne',t],PPlist[["std"]][h,'Ne',t]**2),file=fcon,append=TRUE)
        }
    }

    close(fcon)
}
