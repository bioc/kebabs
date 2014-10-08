#ifndef __FeatureWeightsPosDep_C_H__

#define __FeatureWeightsPosDep_C_H__


#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP getFeatureWeightsPosDepC(SEXP svR, SEXP selSVR, SEXP offsetSVR, SEXP isXStringSetR, SEXP maxSeqLengthR,
                                         SEXP svmIndexR, SEXP weightLimitR, SEXP kernelTypeR, SEXP kR, SEXP mR,
                                         SEXP bioCharsetR, SEXP motifsR, SEXP motifLengthsR, SEXP maxMotifLengthR,
                                         SEXP maxPatternLengthR, SEXP nodeLimitR, SEXP coefsR, SEXP reverseComplementR,
                                         SEXP posSpecR, SEXP minPosR, SEXP maxPosR, SEXP normalizedR, SEXP lowercaseR,
                                         SEXP unmappedR);

extern "C" {

void freeHeapFeatureWeightsC();

}

#endif

#endif

