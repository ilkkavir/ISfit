writeEascii <- function( Elist , fname='Efield.dat' ){
    #
    # Write electric field estimates and their covariances in ascii files
    #
    #
    #
    #
    # IV 2016
    #

    fcon <- file( fname , 'w' )
    nl <- dim(Elist[["E"]])[1]

    cat('#
# Electric field estimates from multistatic incoherent scatter analysis.
#
# The columns are:
#
# 1. unix time at end of integration period (seconds since 1970-01-01)
# 2. Electric field magnetic East component (V/m)
# 3. Electric field magnetic North component (V/m).
#    Perpendicular to the geomagnetic field,
#    not horizontal in geographic coordinates!
# 4. Variance of the East component (V/m)^2
# 5. Variance of the North component (V/m)^2
# 6. Error covariance of the two components (V/m)^2
#
# IV 2016
# 
#Timestamp     East    North    Var(East)    Var(North)  Cov(Eeast,North)\n',file=fcon,append=FALSE)

    for( l in seq(nl) ){
        cat(sprintf("%14.3f %7.4f %7.4f %12.9f %12.9f %12.9f\n",Elist[["time"]][l],Elist[["E"]][l,1],Elist[["E"]][l,2],Elist[["Ecov"]][l,1,1],Elist[["Ecov"]][l,2,2],Elist[["Ecov"]][l,1,2]),file=fcon,append=TRUE)
    }
    
    close(fcon)

}
