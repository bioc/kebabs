#ifndef __Mismatch_C_H__

#define __Mismatch_C_H__


#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP mismatchKernelMatrixC(SEXP xR, SEXP yR, SEXP selXR, SEXP selYR, SEXP sizeXR, SEXP sizeYR,
                                      SEXP isXStringSetR, SEXP symmetricR, SEXP bioCharsetR, SEXP lowercaseR,
                                      SEXP unmappedR, SEXP maxSeqLengthR, SEXP kR, SEXP mR, SEXP normalizedR,
                                      SEXP presenceR);

#else

extern "C"
{
    #include "ByteStringVector.h"
}

void freeHeapMismatch();

void genPredProfileMismatch(Rcpp::NumericMatrix pprof, ByteStringVector x, Rcpp::IntegerVector selX,
                            int numSamples, ByteStringVector annCharset, ByteStringVector annX,
                            int maxSeqLength, bool unmapped, int kernelType, int k, int m, int bioCharset,
                            Rcpp::NumericMatrix featureWeights, int svmIndex, bool lowercase, bool normalized,
                            bool presence);

RcppExport SEXP genExplRepMismatch(ByteStringVector x, int sizeX, Rcpp::IntegerVector selX, int maxSeqLength,
                                   int bioCharset, ByteStringVector features, int k, int m, bool presence,
                                   bool normalized, bool unmapped, bool lowercase, bool useRowNames,
                                   bool useColNames, bool zeroFeatures, bool sparse);

void mismatchFreeHeapC();

#endif

#endif
