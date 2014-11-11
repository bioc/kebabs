//
// C Routines for Prediction
//
// Source : PredictionC.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014   J o h a n n e s  P a l m e
//

#include "Kebabs.h"
#include "KernelUtils.h"
#include "khash.h"
#include "FeatureVectorC.h"
#include "SpectrumC.h"
#include "MismatchC.h"
#include "GappyPairC.h"
#include "MotifC.h"
#include "SparseMatrixHash.h"

#if __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

KHASH_MAP_INIT_INT(fwa32, uint32_t)
KHASH_MAP_INIT_INT64(fwa64, uint32_t)
static khash_t(fwa32) *pAccessHMap32;
static khash_t(fwa64) *pAccessHMap64;

#if __clang__
#pragma clang diagnostic pop
#endif

using namespace Rcpp;

uint64_t * featureNamesToIndex(SEXP featureNames, int numFeatures, ByteStringVector annCharset, IntegerVector annotationIndexMap,
                               int kernelType, int k, int m, void **motifTree, int *freeNode, ByteStringVector motifs,
                               IntegerVector *motifLengths, int maxMotifLength, int maxPatternLength, int nodeLimit,
                               bool reverseComplement, struct alphaInfo *alphaInf)
{
    switch (kernelType)
    {
        case SPECTRUM:

            return(featureNamesToIndexSpectrum(featureNames, numFeatures, annCharset, annotationIndexMap,
                                               k, reverseComplement, alphaInf));
            break;

        case MISMATCH:

//            return(featureNamesToIndexMismatch(featureNames, annCharset, k, m, motifs,
//                                        alphaInf));
            break;

        case GAPPY_PAIR:

            return(featureNamesToIndexGappyPair(featureNames, numFeatures, annCharset, annotationIndexMap,
                                                k, m, reverseComplement, alphaInf));
            break;

        case MOTIF:

            return(featureNamesToIndexMotif(featureNames, numFeatures, motifTree, freeNode, motifs, motifLengths,
                                            maxMotifLength, maxPatternLength, nodeLimit, alphaInf));
            break;
    }

    return(NULL);
}

void genFeatureVectorsPosDep(ByteStringVector x, int sizeX, IntegerVector selX, IntegerVector offsetX,
                             ByteStringVector annX, ByteStringVector annCharset, int maxSeqLength,
                             int kernelType, int k, int m, void **pMotifTree, int *freeNode,
                             ByteStringVector motifs, IntegerVector *motifLengths, int maxMotifLength,
                             int maxPatternLength, int nodeLimit, struct alphaInfo *alphaInf,
                             uint64_t dimFeatureSpace, bool presence, bool normalized, bool unmapped,
                             bool reverseComplement, bool posSpecific, int sortType, int numPositions,
                             uint64_t **startIndex, void **featVectorIndex, int32_t **featVectorValue,
                             uint32_t **kernelValue, int *indexSize)
{
    switch (kernelType)
    {
        case SPECTRUM:

            genFeatVectorsPosDepSpectrum(x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength, k, alphaInf,
                                         dimFeatureSpace, presence, normalized, unmapped, reverseComplement,
                                         posSpecific, KBS_SORT_BY_POSITION, numPositions, startIndex,
                                         featVectorIndex, featVectorValue, kernelValue, indexSize);
            break;

        // mismatch does not support position dependent kernels

        case GAPPY_PAIR:

            genFeatVectorsGappyPair(x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength, k, m, alphaInf,
                                    dimFeatureSpace, presence, normalized, unmapped, reverseComplement,
                                    posSpecific, KBS_SORT_BY_POSITION, numPositions, startIndex,
                                    featVectorIndex, featVectorValue, kernelValue, indexSize);
            break;

        case MOTIF:

            genFeatVectorsMotif(x, sizeX, selX, offsetX, maxSeqLength, pMotifTree, freeNode, motifs,
                                *motifLengths, maxMotifLength, maxPatternLength, nodeLimit, alphaInf,
                                presence, normalized, posSpecific, KBS_SORT_BY_POSITION, startIndex,
                                (uint32_t **) featVectorIndex, featVectorValue, kernelValue);
            break;
    }

}

RcppExport SEXP getPosDepPredOrProfC(SEXP featureWeightsR, SEXP weightLimitR, SEXP bR, SEXP xR, SEXP isXStringSetR,
                                     SEXP selXR, SEXP sizeXR, SEXP offsetXR, SEXP maxSeqLengthR, SEXP bioCharsetR,
                                     SEXP kernelTypeR, SEXP kR, SEXP mR, SEXP motifsR, SEXP motifLengthsR,
                                     SEXP maxMotifLengthR, SEXP maxPatternLengthR, SEXP nodeLimitR, SEXP posSpecificR,
                                     SEXP distWeightR, SEXP lowercaseR, SEXP unmappedR, SEXP reverseComplementR,
                                     SEXP normalizedR, SEXP getPredProfileR, SEXP firstWeightedPosR, SEXP minPosR,
                                     SEXP maxPosR)
{
    const void *vmax;
    int i, j, l, *slot_i, *slot_p, maxMotifLength, maxPatternLength, numRows, numCols, freeNode;
    int distWeightLength, position, numPositions, offset, indexSize, result, rowIndex,upper1, upper2;
    int lengthPred, nrowProf, ncolProf, nrowFeatWeights, ncolFeatWeights, numDots, numUnweightedPos;
    uint32_t *kernelValue, *unweightedPositions, **unweightedPos, startIndexUnweighted;
    int32_t *featVectorValue;
    uint64_t dimFeatureSpace, *rowIndices, *startIndex, key, featureIndex;
    double *slot_x, part, normFactor, sumAlphas, weight;
    bool bit64 = FALSE, implicitPosition = FALSE, presence = FALSE;
    void *featVectorIndex, *pMotifTree;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    khiter_t iter;
    khash_t(fwa32) *hmap32;
    khash_t(fwa64) *hmap64;
    SEXP rownames;

    ByteStringVector x;
    ByteStringVector motifs;
    IntegerVector *pMotifLengths;
    ByteStringVector dummy;
    dummy.length = 0;
    distWeightLength = 0;
    pMotifTree = NULL;
    hmap32 = NULL;
    hmap64 = NULL;

    int k = as<int>(kR);
    int m = as<int>(mR);
    int sizeX = as<int>(sizeXR);
    int minPos = as<int>(minPosR);
    int maxPos = as<int>(maxPosR);
    int firstWeightedPos = as<int>(firstWeightedPosR);
    int kernelType = as<int>(kernelTypeR);
    int bioCharset = as<int>(bioCharsetR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    int nodeLimit = as<int>(nodeLimitR);
    double b = as<double>(bR);
    double weightLimit = as<double>(weightLimitR);
    bool reverseComplement = as<bool>(reverseComplementR);
    bool getPredProfile = as<bool>(getPredProfileR);
    bool isXStringSet = as<bool>(isXStringSetR);
    bool posSpecific = as<bool>(posSpecificR);
    bool normalized = as<bool>(normalizedR);
    bool lowercase = as<bool>(lowercaseR);
    bool unmapped = as<bool>(unmappedR);

    numPositions = maxPos - minPos + 1;
    unweightedPos = &unweightedPositions;

    if (!getPredProfile)
    {
        lengthPred = sizeX;
        nrowProf = 0;
        ncolProf = 0;
    }
    else
    {
        lengthPred = 0;
        nrowProf = sizeX;
        ncolProf = numPositions;
    }

    NumericVector pred(lengthPred);
    NumericMatrix pprof(nrowProf,ncolProf);

    vmax = vmaxget();

    if (kernelType == SPECTRUM)
        implicitPosition = TRUE;
    else
    {
        implicitPosition = FALSE;

        if (kernelType == MOTIF)
            indexSize = 4;
    }

    if (!isXStringSet)
        x = charVector2ByteStringVec(xR);
    else
        x = XStringSet2ByteStringVec(xR);

    if (kernelType == MOTIF)
    {
        IntegerVector motifLengths(motifLengthsR);
        pMotifLengths = &motifLengths;
        maxMotifLength = as<int>(maxMotifLengthR);
        maxPatternLength = as<int>(maxPatternLengthR);
        motifs = charVector2ByteStringVec(motifsR);
    }
    else
    {
        pMotifLengths = NULL;
        maxMotifLength = 0;
        maxPatternLength = 0;
        motifs.length = 0;
    }

    IntegerVector selX(selXR);
    IntegerVector offsetX(offsetXR);
    NumericVector distWeight(distWeightR);
    IntegerVector unweightedPosStart(motifs.length + 1);
    IntegerVector noAnnot(0);
    IntegerVector sel(1);

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    // annotation not relevant for pos dependent kernels
    dimFeatureSpace = getDimFeatureSpace(kernelType, k, m, alphaInf.numAlphabetChars, 0,
                                         motifs.length, maxMotifLength);

    slot_i = INTEGER(GET_SLOT(featureWeightsR, Rf_install(MATRIX_I_SLOT)));
    slot_p = INTEGER(GET_SLOT(featureWeightsR, Rf_install(MATRIX_P_SLOT)));
    slot_x = REAL(GET_SLOT(featureWeightsR, Rf_install(MATRIX_X_SLOT)));
    numRows = INTEGER(GET_SLOT(featureWeightsR, Rf_install(MATRIX_DIM_SLOT)))[0];
    numCols = INTEGER(GET_SLOT(featureWeightsR, Rf_install(MATRIX_DIM_SLOT)))[1];
    rownames = VECTOR_ELT(GET_SLOT(featureWeightsR, Rf_install(MATRIX_DIMNAMES_SLOT)), 0);

    rowIndices = featureNamesToIndex(rownames, numRows, dummy, noAnnot, kernelType, k, m, &pMotifTree, &freeNode,
                                     motifs, pMotifLengths, maxMotifLength, maxPatternLength, nodeLimit,
                                     reverseComplement, &alphaInf);

    if (!posSpecific)
    {
        nrowFeatWeights = numRows;
        ncolFeatWeights = maxPos - minPos + 1;
        distWeightLength = distWeight.size();
    }
    else
    {
        nrowFeatWeights = 0;
        ncolFeatWeights = 0;
    }

    // alloc feature weights matrix for distance weighted kernels
    NumericMatrix featureWeights(nrowFeatWeights, ncolFeatWeights);

    if (posSpecific)
    {
        bit64 = (dimFeatureSpace * (maxPos - minPos + 1)) > MAXUINT32;

        // add access hashmap to sparse feature weights matrix
        if (bit64)
        {
            hmap64 = (khash_t(fwa64) *) generateAccessHashmap(numRows, numCols, rowIndices,
                                                              dimFeatureSpace, firstWeightedPos, slot_i,
                                                              slot_p, slot_x, bit64);
            pAccessHMap64 = hmap64;
        }
        else
        {
            hmap32 = (khash_t(fwa32) *) generateAccessHashmap(numRows, numCols, rowIndices,
                                                              dimFeatureSpace, firstWeightedPos, slot_i,
                                                              slot_p, slot_x, bit64);
            pAccessHMap32 = hmap32;
        }
    }
    else
    {
        // generate feature index to row index hashmap
        hmap64 = kh_init(fwa64);
        pAccessHMap64 = hmap64;

        for (i = 0; i < numRows; i++)
        {
            iter = kh_put(fwa64, hmap64, rowIndices[i], &result);

            if (result != -1)
                kh_value(hmap64, iter) = i;
        }

        // generate full weight matrix - walk through dgCMatrix columnwise
        for (i = 0; i < numCols; i++)
        {
            for (j = slot_p[i]; j < slot_p[i + 1]; j++)
            {
                // copy weight vector forward and reverse around position
                // clipped to relevant index range
                rowIndex = slot_i[j];
                sumAlphas = slot_x[j];

                if (i + 1 >= distWeightLength)
                    upper1 = distWeightLength;
                else
                    upper1 = i + 1;

                if (i + distWeightLength <= ncolFeatWeights)
                    upper2 = distWeightLength;
                else
                    upper2 = ncolFeatWeights - i;

                for (l = 0; l < upper1; l++)
                    featureWeights(rowIndex, i - l) += sumAlphas * distWeight[l];

                for (l = 1; l < upper2; l++)
                    featureWeights(rowIndex, i + l) += sumAlphas * distWeight[l];
            }
        }
    }

    // generate unweighted positions for motif kernel for prediction profile
    if (kernelType == MOTIF && getPredProfile)
        findUnweightedPositions(motifs, &unweightedPosStart, unweightedPos);

    normFactor = 1;

    for (i = 0; i < sizeX; i++)
    {
        sel[0] = selX[i]; // pass samples individually

        if (offsetX.length() > 0)
            offset = offsetX[sel[0]];
        else
            offset = 0;

        genFeatureVectorsPosDep(x, 1, sel, offsetX, dummy, dummy, maxSeqLength, kernelType, k, m, &pMotifTree,
                                &freeNode, motifs, pMotifLengths, maxMotifLength, maxPatternLength, nodeLimit,
                                &alphaInf, dimFeatureSpace, presence, normalized, unmapped, reverseComplement,
                                posSpecific, KBS_SORT_BY_POSITION, numPositions, &startIndex, &featVectorIndex,
                                &featVectorValue, &kernelValue, &indexSize);

        if (normalized)
            normFactor = 1.0 / sqrt(kernelValue[0]);

        if (posSpecific)
        {
            if (!getPredProfile)
            {
                // compute predictions
                pred[i] = b;

                for (j = 0; j < (int) startIndex[1]; j++)
                {
                    if (implicitPosition)
                    {
                        switch (indexSize) {
                            case 1:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint8_t *) featVectorIndex)[j];
                                break;
                            case 2:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint16_t *) featVectorIndex)[j];
                                break;
                            case 3:
                            case 4:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint32_t *) featVectorIndex)[j];
                                break;
                            default:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint64_t *) featVectorIndex)[j];
                                break;
                        }
                    }
                    else
                    {
                        switch (indexSize) {
                            case 1:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint8_t *) featVectorIndex)[j];
                                break;
                            case 2:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint16_t *) featVectorIndex)[j];
                                break;
                            case 3:
                            case 4:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint32_t *) featVectorIndex)[j];
                                break;
                            default:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint64_t *) featVectorIndex)[j];
                                break;
                        }
                    }

                    if (bit64)
                    {
                        iter = kh_get(fwa64, hmap64, key);

                        if (iter != kh_end(hmap64))
                            pred[i] += normFactor * slot_x[kh_value(hmap64, iter)];
                    }
                    else
                    {
                        iter = kh_get(fwa32, hmap32, (uint32_t) key);

                        if (iter != kh_end(hmap32))
                            pred[i] += normFactor * slot_x[kh_value(hmap32, iter)];
                    }
                }
            }
            else
            {
                // generate prediction profile
                if (implicitPosition)
                {
                    for (j = 0; j < (int) startIndex[1]; j++)
                    {
                        switch (indexSize) {
                            case 1:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint8_t *) featVectorIndex)[j];
                                break;
                            case 2:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint16_t *) featVectorIndex)[j];
                                break;
                            case 3:
                            case 4:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint32_t *) featVectorIndex)[j];
                                break;
                            default:
                                key = (offset + j + 1) * dimFeatureSpace + ((uint64_t *) featVectorIndex)[j];
                                break;
                        }

                        if (bit64)
                        {
                            iter = kh_get(fwa64, hmap64, key);

                            if (iter != kh_end(hmap64))
                            {
                                part = normFactor * slot_x[kh_value(hmap64, iter)] / k;
                                position = offset + j + 1 - minPos;

                                for (l = position; l < position + k; l++)
                                    pprof(i, l) += part;
                            }
                        }
                        else
                        {
                            iter = kh_get(fwa32, hmap32, (uint32_t) key);

                            if (iter != kh_end(hmap32))
                            {
                                part = normFactor * slot_x[kh_value(hmap32, iter)] / k;
                                position = offset + j + 1 - minPos;

                                for (l = position; l < position + k; l++)
                                    pprof(i, l) += part;
                            }
                        }
                    }
                }
                else
                {
                    for (j = 0; j < (int) startIndex[1]; j++)
                    {
                        switch (indexSize) {
                            case 1:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint8_t *) featVectorIndex)[j];
                                featureIndex = ((uint8_t *) featVectorIndex)[j];
                                break;
                            case 2:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint16_t *) featVectorIndex)[j];
                                featureIndex = ((uint16_t *) featVectorIndex)[j];
                                break;
                            case 3:
                            case 4:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint32_t *) featVectorIndex)[j];
                                featureIndex = ((uint32_t *) featVectorIndex)[j];
                                break;
                            default:
                                key = featVectorValue[j] * dimFeatureSpace + ((uint64_t *) featVectorIndex)[j];
                                featureIndex = ((uint64_t *) featVectorIndex)[j];
                                break;
                        }

                        if (bit64)
                        {
                            iter = kh_get(fwa64, hmap64, key);

                            if (iter != kh_end(hmap64))
                            {
                                position = featVectorValue[j] - minPos;

                                if (kernelType == GAPPY_PAIR)
                                {
                                    part = normFactor * slot_x[kh_value(hmap64, iter)] / 2 * k;
                                    numDots = key % (m + 1);

                                    for (l = position; l < position + k; l++)
                                        pprof(i, l) += part;

                                    for (l = position + k + numDots; l < position + 2 * k + numDots; l++)
                                        pprof(i, l) += part;
                                }
                                else if (kernelType == MOTIF)
                                {
                                    numUnweightedPos = unweightedPosStart[featureIndex + 1] -
                                                       unweightedPosStart[featureIndex];

                                    if (numUnweightedPos > 0)
                                        startIndexUnweighted = unweightedPosStart[featureIndex];
                                    else
                                        startIndexUnweighted = MAXUINT32;

                                    part = normFactor * slot_x[kh_value(hmap64, iter)] /
                                           ((*pMotifLengths)[featureIndex] - numUnweightedPos);

                                    for (l = position; l < position + (*pMotifLengths)[featureIndex]; l++)
                                    {
                                        while ((startIndexUnweighted <
                                                (uint32_t) (unweightedPosStart)[featureIndex + 1]) &&
                                               (l - position) > (int) unweightedPositions[startIndexUnweighted])
                                            startIndexUnweighted++;

                                        if (startIndexUnweighted == (uint32_t) unweightedPosStart[featureIndex + 1])
                                            startIndexUnweighted = MAXUINT32;

                                        if (startIndexUnweighted == MAXUINT32 ||
                                            ((int) unweightedPositions[startIndexUnweighted] != (l - position)))
                                            pprof(i, l) += part;
                                    }
                                }
                            }
                        }
                        else
                        {
                            iter = kh_get(fwa32, hmap32, (uint32_t) key);

                            if (iter != kh_end(hmap32))
                            {
                                position = featVectorValue[j] - minPos;

                                if (kernelType == GAPPY_PAIR)
                                {
                                    part = normFactor * slot_x[kh_value(hmap32, iter)] / 2 * k;
                                    numDots = key % (m + 1);

                                    for (l = position; l < position + k; l++)
                                        pprof(i, l) += part;

                                    for (l = position + k + numDots; l < position + 2 * k + numDots; l++)
                                        pprof(i, l) += part;
                                }
                                else if (kernelType == MOTIF)
                                {
                                    numUnweightedPos = unweightedPosStart[featureIndex + 1] -
                                                       unweightedPosStart[featureIndex];

                                    if (numUnweightedPos > 0)
                                        startIndexUnweighted = unweightedPosStart[featureIndex];
                                    else
                                        startIndexUnweighted = MAXUINT32;

                                    part = normFactor * slot_x[kh_value(hmap32, iter)] /
                                           ((*pMotifLengths)[featureIndex] - numUnweightedPos);

                                    for (l = position; l < position + (*pMotifLengths)[featureIndex]; l++)
                                    {
                                        while ((startIndexUnweighted <
                                                (uint32_t) (unweightedPosStart)[featureIndex + 1]) &&
                                               (l - position) > (int) unweightedPositions[startIndexUnweighted])
                                            startIndexUnweighted++;

                                        if (startIndexUnweighted == (uint32_t) unweightedPosStart[featureIndex + 1])
                                            startIndexUnweighted = MAXUINT32;

                                        if (startIndexUnweighted == MAXUINT32 ||
                                            ((int) unweightedPositions[startIndexUnweighted] != (l - position)))
                                            pprof(i, l) += part;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // handling for distance weighted kernels
            if (!getPredProfile)
            {
                // compute predictions
                pred[i] = b;

                for (j = 0; j < (int) startIndex[1]; j++)
                {
                    switch (indexSize) {
                        case 1:
                            featureIndex = ((uint8_t *) featVectorIndex)[j];
                            break;
                        case 2:
                            featureIndex = ((uint16_t *) featVectorIndex)[j];
                            break;
                        case 3:
                        case 4:
                            featureIndex = ((uint32_t *) featVectorIndex)[j];
                            break;
                        default:
                            featureIndex = ((uint64_t *) featVectorIndex)[j];
                            break;
                    }

                    iter = kh_get(fwa64, hmap64, featureIndex);

                    if (iter != kh_end(hmap64))
                    {
                        if (implicitPosition)
                            position = offset + j - minPos + 1;
                        else
                            position = featVectorValue[j] - minPos;

                        weight = featureWeights(kh_value(hmap64, iter), position);

                        // late weight pruning for distance weighted kernels
                        if (fabs(weight) > weightLimit)
                            pred[i] += normFactor * weight;
                    }
                }
            }
            else
            {
                // generate prediction profile
                for (j = 0; j < (int) startIndex[1]; j++)
                {
                    switch (indexSize) {
                        case 1:
                            featureIndex = ((uint8_t *) featVectorIndex)[j];
                            break;
                        case 2:
                            featureIndex = ((uint16_t *) featVectorIndex)[j];
                            break;
                        case 3:
                        case 4:
                            featureIndex = ((uint32_t *) featVectorIndex)[j];
                            break;
                        default:
                            featureIndex = ((uint64_t *) featVectorIndex)[j];
                            break;
                    }

                    iter = kh_get(fwa64, hmap64, featureIndex);

                    if (iter != kh_end(hmap64))
                    {
                        if (implicitPosition)
                        {
                            position = offset + j + 1 - minPos;
                            weight = featureWeights(kh_value(hmap64, iter), position);

                            // late weight pruning for distance weighted kernels
                            if (fabs(weight) > weightLimit)
                            {
                                part = normFactor * featureWeights(kh_value(hmap64, iter), position) / k;

                                for (l = position; l < position + k; l++)
                                    pprof(i, l) += part;
                            }
                        }
                        else
                        {
                            position = featVectorValue[j] - minPos;
                            weight = featureWeights(kh_value(hmap64, iter), position);

                            // late weight pruning for distance weighted kernels
                            if (fabs(weight) > weightLimit)
                            {
                                if (kernelType == GAPPY_PAIR)
                                {
                                    part = normFactor * weight / 2 * k;
                                    numDots = featureIndex % (m + 1);

                                    for (l = position; l < position + k; l++)
                                        pprof(i, l) += part;

                                    for (l = position + k + numDots; l < position + 2 * k + numDots; l++)
                                        pprof(i, l) += part;
                                }
                                else if (kernelType == MOTIF)
                                {
                                    numUnweightedPos = unweightedPosStart[featureIndex + 1] -
                                    unweightedPosStart[featureIndex];

                                    if (numUnweightedPos > 0)
                                        startIndexUnweighted = unweightedPosStart[featureIndex];
                                    else
                                        startIndexUnweighted = MAXUINT32;

                                    part = normFactor * weight / ((*pMotifLengths)[featureIndex] - numUnweightedPos);

                                    for (l = position; l < position + (*pMotifLengths)[featureIndex]; l++)
                                    {
                                        while ((startIndexUnweighted <
                                                (uint32_t) (unweightedPosStart)[featureIndex + 1]) &&
                                               (l - position) > (int) (*unweightedPos)[startIndexUnweighted])
                                            startIndexUnweighted++;

                                        if (startIndexUnweighted == (uint32_t) unweightedPosStart[featureIndex + 1])
                                            startIndexUnweighted = MAXUINT32;

                                        if (startIndexUnweighted == MAXUINT32 ||
                                            ((int) (*unweightedPos)[startIndexUnweighted] != (l - position)))
                                            pprof(i, l) += part;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (posSpecific)
    {
        if (bit64)
        {
            kh_destroy(fwa64, hmap64);
            pAccessHMap64 = NULL;
        }
        else
        {
            kh_destroy(fwa32, hmap32);
            pAccessHMap32 = NULL;
        }
    }
    else
    {
        kh_destroy(fwa64, hmap64);
        pAccessHMap64 = NULL;
    }

    // dealloc of C heap done by on.exit hook

    vmaxset(vmax);

    if (getPredProfile)
        return(pprof);
    else
        return(pred);
}
