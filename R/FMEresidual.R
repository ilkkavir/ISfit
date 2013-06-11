FMEresidual <- function( p , directTheory , measData , measVar , aprioriTheory , aprioriMeas , invAprioriCovar , ... )
    {
        return(directTheory( p , ... )/sqrt(measVar))
        
    }
