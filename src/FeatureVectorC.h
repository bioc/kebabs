#ifndef __FeatureVector_C_H__

#define __FeatureVector_C_H__


#include "SpectrumC.h"

extern "C"
{
    #include "ByteStringVector.h"
}

// generate the sparse feature vector for one or multiple samples for a given kernel
// the representation is compressed, the value array contains the position for
// position dependent kernels
template<typename T>
void generateFeatureVectorsT(T maxUnsignedIndex, ByteStringVector x, int sizeX, Rcpp::IntegerVector selX,
                             Rcpp::IntegerVector offsetX, ByteStringVector annX, ByteStringVector annCharset,
                             int maxSeqLength, int kernelType, int k, int m, struct alphaInfo *alphaInf,
                             uint64_t dimFeatureSpace, bool presence, bool normalized, bool unmapped,
                             bool posSpecific, int sortType, T **featVectorIndex,
                             int32_t **featVectorValue, uint32_t **startIndex)
{
    switch (kernelType)
    {
        case SPECTRUM:
        {
            genFeatureVectorsSpectrum(maxUnsignedIndex, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
                                      k, alphaInf, presence, normalized, unmapped, posSpecific, sortType,
                                      startIndex, featVectorIndex, featVectorValue);
            break;
        }
            //        case MIXED_SPECTRUM:
            //            break;

        case MISMATCH:
        {
            //            genFeatureVectorsMismatch(maxUnsignedIndex, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
            //                                      k, m, presence, normalized, unmapped, elemIndex, featVectorIndex,
            //                                      featVectorValue, startIndex));
            break;
        }
            //        case WEIGHTED_DEGREE:
            //            break;

        case MOTIF:
        {
            //            motifs = charVector2ByteStringVec(motifsR);
            //            IntegerVector motifLengths(motifLengthsR);
            //            int maxMotifLength = as<int>(maxMotifLengthR);
            //            int maxPatternLength = as<int>(maxPatternLengthR);
            //            int nodeLimit = as<int>(nodeLimitR);

            //            genFeatureVectorsMotif(maxUnsignedIndex, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
            //                                   motifs, motifLengths, maxMotifLength, maxPatternLength, nodeLimit,
            //                                   presence, normalized, unmapped, elemIndex, featVectorIndex,
            //                                   featVectorValue, startIndex));
            break;
        }

        case GAPPY_PAIR:
        {
            //            genFeatureVectorsGappyPair(maxUnsignedIndex, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
            //                                       k, m, presence, normalized, unmapped, elemIndex, featVectorIndex,
            //                                       featVectorValue, startIndex));
            break;
        }
    }

    return;
}

#endif
