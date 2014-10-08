//
// C Routines for Position Dependent Feature Weights
//
// Source : FeatureWeightsPosDepC.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014   J o h a n n e s  P a l m e
//

#include <Rcpp.h>
#include "Kebabs.h"
#include "KernelUtils.h"
#include "SpectrumC.h"
#include "MismatchC.h"
#include "GappyPairC.h"
#include "MotifC.h"

extern "C"
{
    #include "ByteStringVector.h"
    #include "khash.h"
}

using namespace Rcpp;

// KHASH_MAP_INIT_INT64(pdfw, double) in KernelUtils.h
// KHASH_MAP_INIT_INT64(pdfi, double) in KernelUtils.h
static khash_t(pdfw) *pPDFeatWeightsHMap;
static khash_t(pdfi) *pPDFeatIndexHMap;

template<typename T>
bool getWeightsPerPosition(T maxUnSignedIndex, SEXP *pdFeatWeights, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap,
                           bool posSpecific, int k, uint64_t dimFeatureSpace, uint64_t numKeys, T *keys)
{
    int i, numCols, currElem, colIndex, currCol;
    uint64_t featIndex;
    khiter_t iter;
    SEXP slot_i, slot_p, slot_x;

    slot_i = GET_SLOT(*pdFeatWeights, Rf_install(MATRIX_I_SLOT));
    slot_p = GET_SLOT(*pdFeatWeights, Rf_install(MATRIX_P_SLOT));
    slot_x = GET_SLOT(*pdFeatWeights, Rf_install(MATRIX_X_SLOT));
    numCols = INTEGER(GET_SLOT(*pdFeatWeights, Rf_install(MATRIX_DIM_SLOT)))[1];

    // same representation for feature weights both for position specific and position dependent kernels
    // for distance weighted kernels the actual feature weights per position are computed before prediction
    // profile generation when the relevant position range is known for prediction samples
    currCol = -1;
    currElem = 0;

    for (i = 0; i < (int) numKeys; i++)
    {
        iter = kh_get(pdfw, pdfwmap, keys[i]);

        if (iter == kh_end(pdfwmap))
        {
            // $$$ TODO remove print
            Rprintf("key %llu not found in hashmap during determination of feature weights\n", keys[i]);
            return(FALSE);
        }

        colIndex = keys[i] / dimFeatureSpace;
        featIndex = keys[i] % dimFeatureSpace;

        while (currCol < colIndex)
            INTEGER(slot_p)[++currCol] = currElem;

        REAL(slot_x)[currElem] = kh_value(pdfwmap, iter);

        iter = kh_get(pdfi, pdfimap, featIndex);

        if (iter == kh_end(pdfimap))
        {
            // $$$ TODO remove print
            Rprintf("pattern %llu not found in hashmap during determination of feature weights\n", featIndex);
            return(FALSE);
        }

        INTEGER(slot_i)[currElem++] = kh_value(pdfimap, iter);
    }

    while (currCol <= numCols)
        INTEGER(slot_p)[++currCol] = currElem;

    return(TRUE);
}

bool getPDFeatureWeights(SEXP *pdFeatWeights, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, bool posSpecific,
                         int k, uint64_t dimFeatureSpace, int featIndexSize, uint64_t numKeys, void *keys)
{
    uint8_t maxUIndex8 = MAXUINT8;
    uint16_t maxUIndex16 = MAXUINT16;
    uint32_t maxUIndex32 = MAXUINT32;
    uint64_t maxUIndex64 = MAXUINT64;


    switch (featIndexSize)
    {
        case 1:
        {
            return(getWeightsPerPosition(maxUIndex8, pdFeatWeights, pdfwmap, pdfimap, posSpecific, k, dimFeatureSpace,
                                         numKeys, (uint8_t *) keys));
        }

        case 2:
        {
            return(getWeightsPerPosition(maxUIndex16, pdFeatWeights, pdfwmap, pdfimap, posSpecific, k, dimFeatureSpace,
                                         numKeys, (uint16_t *) keys));
        }

        case 3:
        case 4:
        {
            return(getWeightsPerPosition(maxUIndex32, pdFeatWeights, pdfwmap, pdfimap, posSpecific, k, dimFeatureSpace,
                                         numKeys, (uint32_t *) keys));
        }

        default:
        {
            return(getWeightsPerPosition(maxUIndex64, pdFeatWeights, pdfwmap, pdfimap, posSpecific, k, dimFeatureSpace,
                                         numKeys, (uint64_t *) keys));
        }
    }
}

// function returns numeric matrix
RcppExport SEXP getFeatureWeightsPosDepC(SEXP svR, SEXP selSVR, SEXP offsetSVR, SEXP isXStringSetR, SEXP maxSeqLengthR,
                                         SEXP svmIndexR, SEXP weightLimitR, SEXP kernelTypeR, SEXP kR, SEXP mR,
                                         SEXP bioCharsetR, SEXP motifsR, SEXP motifLengthsR, SEXP maxMotifLengthR,
                                         SEXP maxPatternLengthR, SEXP nodeLimitR, SEXP coefsR, SEXP reverseComplementR,
                                         SEXP posSpecR, SEXP minPosR, SEXP maxPosR, SEXP normalizedR, SEXP lowercaseR,
                                         SEXP unmappedR)
{
    ByteStringVector sv;
    ByteStringVector motifs;
    int sizeSV, indexSize;
    uint64_t dimFeatureSpace, temp, numKeys;
    void *keys = NULL;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    khash_t(pdfw) *pdfwmap;
    khash_t(pdfi) *pdfimap;
    SEXP *pPDFeatWeights, pdFeatWeights;

    int k = as<int>(kR);
    int m = as<int>(mR);
    int kernelType = as<int>(kernelTypeR);
    int bioCharset = as<int>(bioCharsetR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    int minPos = as<int>(minPosR);
    int maxPos = as<int>(maxPosR);
    int maxMotifLength = as<int>(maxMotifLengthR);
    double weightLimit = as<double>(weightLimitR);
    bool isXStringSet = as<bool>(isXStringSetR);
    bool lowercase = as<bool>(lowercaseR);
    bool normalized = as<bool>(normalizedR);
    bool unmapped = as<bool>(unmappedR);
    bool posSpecific = as<bool>(posSpecR);
    bool reverseComplement = as<bool>(reverseComplementR);
    IntegerVector selSV(selSVR);
    IntegerVector offsetSV(offsetSVR);
    NumericVector coefs(coefsR);

    sizeSV = selSV.size();
    motifs.length = 0;
    pPDFeatWeights = &pdFeatWeights;

    if (kernelType == MOTIF)
        motifs = charVector2ByteStringVec(motifsR);

    if (isXStringSet)
        sv = XStringSet2ByteStringVec(svR);
    else
        sv = charVector2ByteStringVec(svR);

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    // no annotation for position dependent kernels
    dimFeatureSpace = getDimFeatureSpace(kernelType, k, m, alphaInf.numAlphabetChars, 0,
                                         motifs.length, maxMotifLength);

    indexSize = 1;
    temp = dimFeatureSpace * (maxPos - minPos + 1) - 1;

    while ((temp >>= 8) > 0)
        indexSize++;

    pdfwmap = kh_init(pdfw);
    pPDFeatWeightsHMap = pdfwmap;
    pdfimap = kh_init(pdfi);
    pPDFeatIndexHMap = pdfimap;

    // get sparse feature vector
    switch (kernelType)
    {
        case SPECTRUM:
        {
            getFeaturesOfSVSpectrum(&pPDFeatWeights, pdfwmap, pdfimap, sv, sizeSV, selSV, offsetSV,
                                    maxSeqLength, coefs, reverseComplement, posSpecific, weightLimit,
                                    k, minPos, maxPos, dimFeatureSpace, &alphaInf, normalized,
                                    indexSize, &numKeys, &keys);
            break;
        }

        // case MIXED_SPECTRUM:
        //     break;

        // case MISMATCH:
        // mismatch does not support position dependency

        // case WEIGHTED_DEGREE:
        //     break;

        case GAPPY_PAIR:
        {
            getFeaturesOfSVGappyPair(&pPDFeatWeights, pdfwmap, pdfimap, sv, sizeSV, selSV, offsetSV,
                                     maxSeqLength, coefs, reverseComplement, posSpecific, weightLimit,
                                     k, m, minPos, maxPos, dimFeatureSpace, &alphaInf, normalized,
                                     indexSize, &numKeys, &keys);
            break;
        }

        case MOTIF:
        {
            motifs = charVector2ByteStringVec(motifsR);
            IntegerVector motifLengths(motifLengthsR);
            int maxMotifLength = as<int>(maxMotifLengthR);
            int maxPatternLength = as<int>(maxPatternLengthR);
            int nodeLimit = as<int>(nodeLimitR);

            getFeaturesOfSVMotif(&pPDFeatWeights, pdfwmap, pdfimap, sv, sizeSV, selSV, offsetSV,
                                 maxSeqLength, coefs, posSpecific, weightLimit, motifs, motifLengths,
                                 maxMotifLength, maxPatternLength, nodeLimit, minPos, maxPos,
                                 dimFeatureSpace, &alphaInf, normalized, indexSize, &numKeys, &keys);
            break;
        }
    }

    if (pPDFeatWeights == NULL)
        return(R_NilValue);

    if (getPDFeatureWeights(pPDFeatWeights, pdfwmap, pdfimap, posSpecific, k, dimFeatureSpace, indexSize,
                            numKeys, keys))
        return(*pPDFeatWeights);
    else
        return(R_NilValue);
}

extern "C" {

void freeHeapFeatureWeightsC()
{
    if (pPDFeatWeightsHMap != NULL)
    {
        kh_destroy(pdfw, pPDFeatWeightsHMap);
        pPDFeatWeightsHMap = NULL;
    }

    if (pPDFeatIndexHMap != NULL)
    {
        kh_destroy(pdfi, pPDFeatIndexHMap);
        pPDFeatIndexHMap = NULL;
    }
}

}
