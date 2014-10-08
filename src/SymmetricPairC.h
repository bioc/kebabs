#ifndef __SymmetricPair_C_H__

#define __SymmetricPair_C_H__


#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP symmetricPairKernelC(SEXP kmSIR, SEXP selXR, SEXP selYR,
                                     SEXP sizeXR, SEXP sizeYR, SEXP kernelTypeR,
                                     SEXP symmetricR);

#endif

#endif
