#ifndef __Motif_C_H__

#define __Motif_C_H__


#include <Rcpp.h>

# ifdef __R_Interfaces_Only__

RcppExport SEXP validateMotifsC(SEXP motifsR, SEXP motifLengthsR);

RcppExport SEXP motifKernelMatrixC(SEXP xR, SEXP yR, SEXP selXR, SEXP selYR, SEXP sizeXR, SEXP sizeYR,
                                   SEXP isXStringSetR, SEXP symmetricR, SEXP offsetXR, SEXP offsetYR,
                                   SEXP annCharsetR, SEXP annXR, SEXP annYR, SEXP motifsR, SEXP motifLengthsR,
                                   SEXP nodeLimitR, SEXP maxMotifLengthR, SEXP maxPatternLengthR, SEXP bioCharsetR,
                                   SEXP ignoreLowerR, SEXP unmappedR, SEXP maxSeqLengthR, SEXP posSpecR,
                                   SEXP distWeightR, SEXP normalizedR, SEXP presenceR);

#else

extern "C"
{
    #include "ByteStringVector.h"
}

void freeHeapMotif();

void genPredProfileMotif(Rcpp::NumericMatrix profiles, ByteStringVector x, Rcpp::IntegerVector selX, int numSamples,
                         ByteStringVector annCharset, ByteStringVector annX, int maxSeqLength, bool unmapped,
                         int kernelType, int k, int m, int bioCharset, Rcpp::NumericMatrix featureWeights,
                         int svmIndex, ByteStringVector motifs, Rcpp::IntegerVector *motifLengths, int maxMotifLength,
                         int maxPatternLength, ByteStringVector fwMotifs, Rcpp::IntegerVector *fwMotifLengths,
                         int fwMaxMotifLength, int fwMaxPatternLength, int nodeLimit, bool lowercase,
                         bool normalized, bool presence);

uint64_t * featureNamesToIndexMotif(SEXP featureNames, int numFeatures, void **pMotifTree, int *freeNode,
                                    ByteStringVector motifs, Rcpp::IntegerVector *motifLengths, int maxMotifLength,
                                    int maxPatternLength, int nodeLimit, struct alphaInfo *alphaInf);

void getFeaturesOfSVMotif(SEXP **pdFeatWeights, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, ByteStringVector x,
                          int sizeX, Rcpp::IntegerVector selX, Rcpp::IntegerVector offsetX, int maxSeqLength,
                          Rcpp::NumericVector coefs, bool posSpecific, double weightLimit, ByteStringVector motifs,
                          Rcpp::IntegerVector motifLengths, int maxMotifLength, int maxPatternLength, int nodeLimit,
                          int minPos, int maxPos, uint64_t dimFeatureSpace, struct alphaInfo *alphaInf, bool normalized,
                          int featIndexSize, uint64_t *numKeys, void **keys);

void genFeatVectorsMotif(ByteStringVector x, int sizeX, Rcpp::IntegerVector selX, Rcpp::IntegerVector offsetX,
                         int maxSeqLength, void **pMotifTree, int *freeNode, ByteStringVector motifs,
                         Rcpp::IntegerVector motifLengths, int maxMotifLength, int maxPatternLength,
                         int nodeLimit, struct alphaInfo *alphaInf, bool presence, bool normalized,
                         bool posSpecific, int sortType, uint64_t **startIndex, uint32_t **featVectorIndex,
                         int32_t **featVectorValue, uint32_t **kernelValue);

RcppExport void findUnweightedPositions(ByteStringVector motifs, Rcpp::IntegerVector *unweightedPosStart,
                                        uint32_t **unweightedPos);

RcppExport SEXP genExplRepMotif(ByteStringVector x, int sizeX, Rcpp::IntegerVector selX, ByteStringVector annCharset,
                                ByteStringVector annX, int maxSeqLength, int bioCharset, ByteStringVector motifs,
                                Rcpp::IntegerVector motifLengths, int maxMotifLength, int maxPatternLength,
                                int nodeLimit, bool presence, bool normalized, bool unmapped, bool lowercase,
                                bool useRowNames, bool useColNames, bool zeroFeatures, bool sparse);

#endif

#endif
