#ifndef __GappyPair_C_H__

#define __GappyPair_C_H__

#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP gappyPairKernelMatrixC(SEXP xR, SEXP yR, SEXP selXR, SEXP selYR, SEXP sizeXR,
                    SEXP sizeYR, SEXP isXStringSetR, SEXP symmetricR, SEXP offsetXR,
                    SEXP offsetYR, SEXP annCharsetR, SEXP annXR, SEXP annYR, SEXP bioCharsetR,
                    SEXP ignoreLowerR, SEXP unmappedR, SEXP maxSeqLengthR, SEXP kR, SEXP mR,
                    SEXP posSpecR, SEXP distWeightR, SEXP normalizedR, SEXP presenceR,
                    SEXP reverseComplementR);

#else

extern "C"
{
    #include "ByteStringVector.h"
    #include "khash.h"
}

void freeHeapGappyPair();

void genPredProfileGappyPair(Rcpp::NumericMatrix pprof, ByteStringVector x, Rcpp::IntegerVector selX, int numSamples,
                             ByteStringVector annCharset, ByteStringVector annX, int maxSeqLength, bool unmapped,
                             bool reverseComplement, int kernelType, int k, int m, int bioCharset,
                             Rcpp::NumericMatrix featureWeights, int svmIndex, bool lowercase, bool normalized,
                             bool presence);

void getFeaturesOfSVGappyPair(SEXP **pdFeatWeights, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, ByteStringVector x,
                              int sizeX, Rcpp::IntegerVector selX, Rcpp::IntegerVector offsetX, int maxSeqLength,
                              Rcpp::NumericVector coefs, bool reverseComplement, bool posSpecific, double weightLimit,
                              int k, int m, int minPos, int maxPos, uint64_t dimFeatureSpace, struct alphaInfo *alphaInf,
                              bool normalized, int featIndexSize, uint64_t *numKeys, void **keys);

void genFeatVectorsGappyPair(ByteStringVector x, int sizeX, Rcpp::IntegerVector selX, Rcpp::IntegerVector offsetX,
                             ByteStringVector annX, ByteStringVector annCharset, int maxSeqLength,
                             int k, int m, struct alphaInfo *alphaInf, uint64_t dimFeatureSpace,
                             bool presence, bool normalized, bool unmapped, bool reverseComplement, bool posSpecific,
                             int sortType, int numPositions, uint64_t **startIndex, void **featVectorIndex,
                             int32_t **featVectorValue, uint32_t **kernelValue, int *indexSize);

uint64_t * featureNamesToIndexGappyPair(SEXP featureNames, int numFeatures, ByteStringVector annCharset,
                                        Rcpp::IntegerVector annotationIndexMap, int k, int m, bool reverseComplement,
                                        struct alphaInfo *alphaInf);

RcppExport SEXP genExplRepGappyPair(ByteStringVector x, int sizeX, Rcpp::IntegerVector selX,
                                    ByteStringVector annCharset, ByteStringVector annX,
                                    int maxSeqLength, int bioCharset, ByteStringVector features,
                                    int k, int m, bool presence, bool reverseComplement,
                                    bool normalized, bool unmapped, bool lowercase, bool useRowNames,
                                    bool useColNames, bool zeroFeatures, bool sparse);

#endif

#endif
