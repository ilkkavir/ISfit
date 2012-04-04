#include <stdlib.h>

void crossprods(double *dirtheDataR, double *dirtheDataI, int *nData, double *fAmbR, double *fAmbI, int *nFreq, int *iSite, double *sSite, double *cSite)
{

  register int k,l,n,m;

  for( k=0 ; k < *nData ; ++k){
    l = *nFreq * iSite[k];
    m = *nFreq * k;
    dirtheDataR[k] = 0.0;
    dirtheDataI[k] = 0.0;
    for( n=0 ; n < *nFreq ; ++n){
      dirtheDataR[k] = dirtheDataR[k] + fAmbR[ m + n] * sSite[ l + n ];
      dirtheDataI[k] = dirtheDataI[k] + fAmbI[ m + n] * sSite[ l + n ];
    }
    dirtheDataR[k] = dirtheDataR[k] * cSite[ iSite[ k ] ];
    dirtheDataI[k] = dirtheDataI[k] * cSite[ iSite[ k ] ];
  }

}
