#include "ISfit.h"

static const R_CallMethodDef callMethods[2] = {
  { "crossprods" , (DL_FUNC) & crossprods , 6 } ,
  { NULL , NULL , 0}
};

void R_init_ISfit( DllInfo *info )
{
  R_registerRoutines( info , NULL , callMethods , NULL , NULL);
}
