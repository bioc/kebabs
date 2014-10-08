#ifndef __PredictionProfile_C_H__

#define __PredictionProfile_C_H__


#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP generatePredictionProfilesC(SEXP xR, SEXP bioVectorR, SEXP selXR, SEXP numSamplesR,
                        SEXP annCharsetR, SEXP annXR, SEXP maxSeqLengthR, SEXP unmappedR,
                        SEXP reverseComplementR, SEXP kernelTypeR, SEXP kR, SEXP mR, SEXP bioCharsetR,
                        SEXP featureWeightsR, SEXP svmIndexR, SEXP motifsR, SEXP motifLengthsR,
                        SEXP maxMotifLengthR, SEXP maxPatternLengthR, SEXP fwMotifsR,
                        SEXP fwMotifLengthsR, SEXP fwMaxMotifLengthR, SEXP fwMaxPatternLengthR,
                        SEXP nodeLimitR, SEXP lowercaseR, SEXP normalizedR, SEXP presenceR);

#endif

#endif

