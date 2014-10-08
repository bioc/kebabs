//
// C Routines for Prediction Profile Generation
//
// Source : PredictionProfileC.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014   J o h a n n e s  P a l m e
//

#include "Kebabs.h"
#include "KernelUtils.h"
#include "khash.h"
#include "SpectrumC.h"
#include "MismatchC.h"
#include "GappyPairC.h"
#include "MotifC.h"

using namespace Rcpp;

RcppExport SEXP generatePredictionProfilesC(SEXP xR, SEXP bioVectorR, SEXP selXR, SEXP numSamplesR,
                        SEXP annCharsetR, SEXP annXR, SEXP maxSeqLengthR, SEXP unmappedR,
                        SEXP reverseComplementR, SEXP kernelTypeR, SEXP kR, SEXP mR, SEXP bioCharsetR,
                        SEXP featureWeightsR, SEXP svmIndexR, SEXP motifsR, SEXP motifLengthsR,
                        SEXP maxMotifLengthR, SEXP maxPatternLengthR, SEXP fwMotifsR,
                        SEXP fwMotifLengthsR, SEXP fwMaxMotifLengthR, SEXP fwMaxPatternLengthR,
                        SEXP nodeLimitR, SEXP lowercaseR, SEXP normalizedR, SEXP presenceR)
{
    const void *vmax;
    int maxMotifLength, maxPatternLength, fwMaxMotifLength, fwMaxPatternLength;
    ByteStringVector x;
    ByteStringVector annX;
    ByteStringVector annCharset;
    ByteStringVector motifs;
    ByteStringVector fwMotifs;
    IntegerVector *pMotifLengths;
    IntegerVector *pfwMotifLengths;

    int k = as<int>(kR);
    int m = as<int>(mR);
    int numSamples = as<int>(numSamplesR);
    int kernelType = as<int>(kernelTypeR);
    int bioCharset = as<int>(bioCharsetR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    int nodeLimit = as<int>(nodeLimitR);
    int svmIndex = as<int>(svmIndexR);
    bool bioVector = as<bool>(bioVectorR);
    bool normalized = as<bool>(normalizedR);
    bool presence = as<bool>(presenceR);
    bool lowercase = as<bool>(lowercaseR);
    bool unmapped = as<bool>(unmappedR);
    bool reverseComplement = as<bool>(reverseComplementR);

    NumericMatrix pprof(numSamples, maxSeqLength);

    vmax = vmaxget();

    annCharset.length = 0;
    annX.length = 0;

    if (bioVector)
        x = charVector2ByteStringVec(xR);
    else
        x = XStringSet2ByteStringVec(xR);

    if (!Rf_isNull(annXR))
    {
        annCharset = charVector2ByteStringVec(annCharsetR);
        annX = charVector2ByteStringVec(annXR);
    }

    if (kernelType == MOTIF)
    {
        IntegerVector motifLengths(motifLengthsR);
        pMotifLengths = &motifLengths;
        maxMotifLength = as<int>(maxMotifLengthR);
        maxPatternLength = as<int>(maxPatternLengthR);
        motifs = charVector2ByteStringVec(motifsR);
        IntegerVector fwMotifLengths(fwMotifLengthsR);
        pfwMotifLengths = &fwMotifLengths;
        fwMaxMotifLength = as<int>(fwMaxMotifLengthR);
        fwMaxPatternLength = as<int>(fwMaxPatternLengthR);
        fwMotifs = charVector2ByteStringVec(fwMotifsR);
    }
    else
    {
        pMotifLengths = NULL;
        maxMotifLength = 0;
        maxPatternLength = 0;
        motifs.length = 0;
        pfwMotifLengths = NULL;
        fwMaxMotifLength = 0;
        fwMaxPatternLength = 0;
        fwMotifs.length = 0;
    }

    IntegerVector selX(selXR);
    NumericMatrix featureWeights(featureWeightsR);

    switch (kernelType)
    {
        case SPECTRUM:

            genPredProfileSpectrum(pprof, x, selX, numSamples, annCharset, annX, maxSeqLength, unmapped,
                                   reverseComplement, kernelType, k, bioCharset, featureWeights, svmIndex,
                                   lowercase, normalized, presence);
            break;

        case MISMATCH:

            genPredProfileMismatch(pprof, x, selX, numSamples, annCharset, annX, maxSeqLength, unmapped,
                                   kernelType, k, m, bioCharset, featureWeights, svmIndex, lowercase,
                                   normalized, presence);
            break;

        case GAPPY_PAIR:

            genPredProfileGappyPair(pprof, x, selX, numSamples, annCharset, annX, maxSeqLength, unmapped,
                                    reverseComplement, kernelType, k, m, bioCharset, featureWeights, svmIndex,
                                    lowercase, normalized, presence);
            break;

        case MOTIF:

            genPredProfileMotif(pprof, x, selX, numSamples, annCharset, annX, maxSeqLength, unmapped,
                                kernelType, k, m, bioCharset, featureWeights, svmIndex, motifs,
                                pMotifLengths, maxMotifLength, maxPatternLength, fwMotifs, pfwMotifLengths,
                                fwMaxMotifLength, fwMaxPatternLength, nodeLimit, lowercase, normalized,
                                presence);
            break;

    }

    // dealloc of C heap done by on.exit hook

    vmaxset(vmax);

    return(pprof);
}
