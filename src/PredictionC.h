#ifndef __Prediction_C_H__

#define __Prediction_C_H__


#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP getPosDepPredOrProfC(SEXP featureWeightsR, SEXP weightLimitR, SEXP bR, SEXP xR, SEXP isXStringSetR,
                                     SEXP selXR, SEXP sizeXR, SEXP offsetXR, SEXP maxSeqLengthR, SEXP bioCharsetR,
                                     SEXP kernelTypeR, SEXP kR, SEXP mR, SEXP motifsR, SEXP motifLengthsR,
                                     SEXP maxMotifLengthR, SEXP maxPatternLengthR, SEXP nodeLimitR, SEXP posSpecificR,
                                     SEXP distWeightR, SEXP lowercaseR, SEXP unmappedR, SEXP reverseComplementR,
                                     SEXP normalizedR, SEXP getPredProfileR, SEXP firstWeightedPosR, SEXP minPosR,
                                     SEXP maxPosR);

#endif

#endif

