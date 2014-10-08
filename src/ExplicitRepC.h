#ifndef __ExplicitRep_C_H__

#define __ExplicitRep_C_H__

# ifdef __R_Interfaces_Only__

#include <Rcpp.h>

RcppExport SEXP genExplRepC(SEXP xR, SEXP isXStringSetR, SEXP selXR, SEXP annCharsetR, SEXP annXR,
                            SEXP maxSeqLengthR, SEXP kernelTypeR, SEXP kR, SEXP mR, SEXP bioCharsetR,
                            SEXP featuresR, SEXP motifsR, SEXP motifLengthsR, SEXP maxMotifLengthR,
                            SEXP maxPatternLengthR, SEXP nodeLimitR, SEXP presenceR, SEXP reverseComplementR,
                            SEXP normalizedR, SEXP lowercaseR, SEXP unmappedR, SEXP useRowNamesR,
                            SEXP useColNamesR, SEXP zeroFeaturesR, SEXP sparseR);
extern "C" {

void freeHeapCallocsC(SEXP kernelTypeR);

}

#endif

#endif
