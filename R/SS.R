##
## Sum-of-squares
##
##
##
##
##

SS <- function( p , measData , measVar , directTheory , aprioriTheory , aprioriMeas , invAprioriCovar , ... )
    {

        dirtheData <- directTheory( p , ... )

        return(c( chiSquare( measData , measVar , dirtheData , p , aprioriTheory , aprioriMeas , invAprioriCovar)))

    }
