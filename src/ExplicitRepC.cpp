//
// C Routines for Explicit Representation
//
// Source : ExplicitRepC.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014-2015   J o h a n n e s  P a l m e
//


#include <Rcpp.h>
#include "Kebabs.h"
#include "SpectrumC.h"
#include "MismatchC.h"
#include "GappyPairC.h"
#include "MotifC.h"

extern "C"
{
    #include "ByteStringVector.h"
}

using namespace Rcpp;


RcppExport SEXP genExplRepC(SEXP xR, SEXP isXStringSetR, SEXP selXR, SEXP annCharsetR, SEXP annXR,
                            SEXP maxSeqLengthR, SEXP kernelTypeR, SEXP kR, SEXP mR, SEXP bioCharsetR,
                            SEXP featuresR, SEXP motifsR, SEXP motifLengthsR, SEXP maxMotifLengthR,
                            SEXP maxPatternLengthR, SEXP nodeLimitR, SEXP presenceR, SEXP reverseComplementR,
                            SEXP normalizedR, SEXP lowercaseR, SEXP unmappedR, SEXP useRowNamesR,
                            SEXP useColNamesR, SEXP zeroFeaturesR, SEXP sparseR)
{
    ByteStringVector x;
    ByteStringVector annX;
    ByteStringVector annCharset;
    ByteStringVector motifs;
    ByteStringVector features;
    int sizeX;
    int k = as<int>(kR);
    int m = as<int>(mR);
    int kernelType = as<int>(kernelTypeR);
    int bioCharset = as<int>(bioCharsetR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    bool isXStringSet = as<bool>(isXStringSetR);
    bool lowercase = as<bool>(lowercaseR);
    bool normalized = as<bool>(normalizedR);
    bool unmapped = as<bool>(unmappedR);
    bool presence = as<bool>(presenceR);
    bool reverseComplement = as<bool>(reverseComplementR);
    bool useRowNames = as<bool>(useRowNamesR);
    bool useColNames = as<bool>(useColNamesR);
    bool zeroFeatures = as<bool>(zeroFeaturesR);
    bool sparse = as<bool>(sparseR);

    annCharset.length = 0;
    annX.length = 0;
    features.length = 0;
    motifs.length = 0;


    if (isXStringSet)
        x = XStringSet2ByteStringVec(xR);
    else
        x = charVector2ByteStringVec(xR);

    IntegerVector selX(selXR);
    sizeX = selX.size();


    if (!Rf_isNull(featuresR))
        features = charVector2ByteStringVec(featuresR);
    else
        features.length = -1;

    if (!Rf_isNull(annXR))
    {
        annCharset = charVector2ByteStringVec(annCharsetR);
        annX = charVector2ByteStringVec(annXR);
    }

    switch (kernelType)
    {
        case SPECTRUM:
        {
            return(genExplRepSpectrum(x, sizeX, selX, annCharset, annX, maxSeqLength,
                                      bioCharset, features, k, presence, reverseComplement,
                                      normalized, unmapped, lowercase, useRowNames, useColNames,
                                      zeroFeatures, sparse));
            break;
        }
//        case MIXED_SPECTRUM:
//            break;

        case MISMATCH:
            return(genExplRepMismatch(x, sizeX, selX, maxSeqLength, bioCharset, features,
                                      k, m, presence, normalized, unmapped, lowercase,
                                      useRowNames, useColNames, zeroFeatures, sparse));
            break;

//        case WEIGHTED_DEGREE:
//            break;

        case MOTIF:
        {
            motifs = charVector2ByteStringVec(motifsR);
            IntegerVector motifLengths(motifLengthsR);
            int maxMotifLength = as<int>(maxMotifLengthR);
            int maxPatternLength = as<int>(maxPatternLengthR);
            int nodeLimit = as<int>(nodeLimitR);

            return(genExplRepMotif(x, sizeX, selX, annCharset, annX, maxSeqLength, bioCharset,
                                   motifs, motifLengths, maxMotifLength, maxPatternLength,
                                   nodeLimit, presence, normalized, unmapped, lowercase,
                                   useRowNames, useColNames, zeroFeatures, sparse));
            break;
        }

        case GAPPY_PAIR:
        {
            return(genExplRepGappyPair(x, sizeX, selX, annCharset, annX, maxSeqLength,
                                       bioCharset, features, k, m, presence, reverseComplement,
                                       normalized, unmapped, lowercase, useRowNames, useColNames,
                                       zeroFeatures, sparse));
            break;
        }
    }

    return(NULL);
}

extern "C" {

void freeHeapCallocsC(SEXP kernelTypeR)
{
    int kernelType = as<int>(kernelTypeR);

    switch (kernelType)
    {
        case SPECTRUM:
            freeHeapSpectrum();
            break;
        case MISMATCH:
            freeHeapMismatch();
            break;
        case GAPPY_PAIR:
            freeHeapGappyPair();
            break;
        case MOTIF:
            freeHeapMotif();
            break;
    }
}

}
