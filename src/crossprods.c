#include "ISfit.h"

SEXP crossprods( SEXP nData , SEXP fAmb , SEXP nFreq , SEXP iSite , SEXP sSite , SEXP cSite )
{
  int nd = *INTEGER( nData );
  Rcomplex * ambig = COMPLEX( fAmb );
  int nf = *INTEGER( nFreq );
  int    * is = INTEGER( iSite );
  double * ss = REAL( sSite );
  double * cs = REAL( cSite );

  SEXP dirtheData;
  Rcomplex * dData;
  R_len_t  k,l,n,m;

  // allocate the direct theory vector
  PROTECT( dirtheData = allocVector( CPLXSXP , nd ) );

  // a pointer to the allocated vector
  dData = COMPLEX( dirtheData );

  // walk through data points
  for ( k = 0 ; k < nd ; ++k ){

    // initialise to zero
    dData[k].r = .0;
    dData[k].i = .0;

    // walk through frequency points
    for( l = 0 ; l < nf ; ++l ){

      // frequency ambiguity * spectrum * scaling for this site
      dData[k].r += ambig[ k*nf + l ].r * ss[ is[k]*nf + l ] * cs[is[k]];
      dData[k].i += ambig[ k*nf + l ].i * ss[ is[k]*nf + l ] * cs[is[k]];
    }

  }

  UNPROTECT(1);

  return(dirtheData);


}
