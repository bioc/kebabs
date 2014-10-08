#ifndef __KernelUtils_CPP_H__

#define __KernelUtils_CPP_H__


#include <Rcpp.h>
#include "Utils.h"

extern "C"
{
    #include "ByteStringVector.h"
    #include "Kebabs.h"
    #include "khash.h"
    #include "string.h"
}

struct allIndMaps {
    int *dna;
    int *dnaLower;
    int *rna;
    int *rnaLower;
    int *aa;
    int *aaLower;
    int *all;
    int *allLower;
    int *reverse;
    int *unmapped;
    int *reverseUnmapped;
};

struct alphaInfo {
    int                  seqType;
    bool              lowercase;
    bool              unmapped;
    char              *currAlphabet;
    int                  numAlphabetChars;
    int                  maxAlphaIndex;
    const int          *indexMap;
    int                  *reverseIndexMap;
    int               *seqIndexMap;
};

// hash maps for position specific feature weight computation
KHASH_MAP_INIT_INT64(pdfw, double)
KHASH_MAP_INIT_INT64(pdfi, uint32_t)

void getAlphabetInfo(int seqType, bool lowercase, bool unmapped,
                     struct alphaInfo *info, struct allIndMaps *allIndexMaps);

void initAnnotationMaps(ByteStringVector annCharset, Rcpp::IntegerVector *annIndexMap,
                        Rcpp::IntegerVector *revAnnMap);

void initMatrixWithNA(Rcpp::NumericMatrix matrix, int sizeX, int sizeY);

uint64_t getDimFeatureSpace(int kernelType, int k, int m, int numAlphabetChars, int numAnnChars,
                            int numMotifs, int maxMotifLength);

RcppExport SEXP generateEmptyExplicitRep(int sizeX, bool sparse);

// parallel mergesort of 2 arrays
template<typename T>
void merge(T *array, int lower, int middle, int upper, T *buffer)
{
    int i, j, k;
    // move both array sections to buffer
    for (i = middle + 1; i > lower; i--)
        buffer[i-1] = array[i-1];

    // $$$ TODO use memcpy
//    memcpy((void *) &(buffer[lower]), (void *) &(array[lower]), (upper - lower + 1) * sizeof(T));
//    i = lower;

    for (j = middle; j < upper; j++)
        buffer[upper + middle - j] = array[j + 1];

    // sort both array sections back to original array
    for (k = lower; k <= upper; k++)
    {
        if (buffer[j] < buffer[i])
            array[k] = buffer[j--];
        else
            array[k] = buffer[i++];
    }
}

template<typename T>
void mergesort(T *array,int lower, int upper, T *buffer)
{
    if (upper > lower)
    {
        int middle = (upper + lower) / 2;
        mergesort(array, lower, middle, buffer);
        mergesort(array, middle + 1, upper, buffer);
        merge(array, lower, middle, upper, buffer);
    }
}

template<typename T>
void sortArray(T maxUnSignedIndex, T *array, int noCols, int maxNoRows)
{
    int i, end, limit;
    const void *vmax;

    vmax = vmaxget();

    T *buffer = (T *) R_alloc(maxNoRows, sizeof(T));

    for (i = 0; i < noCols; i++)
    {
        end = i * maxNoRows;
        limit = end + maxNoRows;

        while(array[end] != maxUnSignedIndex && end < limit)
            end++;

        end--;

        mergesort(&array[i * maxNoRows], 0, end - i * maxNoRows, buffer);
    }

    vmaxset(vmax);
}

// parallel mergesort of 2 arrays
template<typename T1, typename T2>
void merge2(T1 *array1, T2 *array2, int lower, int middle, int upper, T1 *buffer1, T2 *buffer2)
{
    int i, j, k;
    // move both array sections to buffer
    for (i = middle + 1; i > lower; i--)
    {
        buffer1[i - 1] = array1[i - 1];
        buffer2[i - 1] = array2[i - 1];
    }
    for (j = middle; j < upper; j++)
    {
        buffer1[upper + middle - j] = array1[j + 1];
        buffer2[upper + middle - j] = array2[j + 1];
    }
    // sort both array sections back to original array
    for (k = lower; k <= upper; k++)
    {
        if (buffer1[j] < buffer1[i])
        {
            array1[k] = buffer1[j];
            array2[k] = buffer2[j--];
        }
        else
        {
            array1[k] = buffer1[i];
            array2[k] = buffer2[i++];
        }
    }
}

template<typename T1, typename T2>
void mergesort2(T1 *array1, T2 *array2, int lower, int upper, T1 *buffer1, T2 *buffer2)
{
    if (upper > lower)
    {
        int middle = (upper + lower) / 2;
        mergesort2(array1, array2, lower, middle, buffer1, buffer2);
        mergesort2(array1, array2, middle + 1, upper, buffer1, buffer2);
        merge2(array1, array2, lower, middle, upper, buffer1, buffer2);
    }
}

template<typename T1, typename T2>
void sort2Arrays(T1 maxUnSignedIndex, T1 *array1, T2 *array2, int noCols, int maxNoRows, uint64_t *startIndex)
{
    int i, j;
    const void *vmax;

    vmax = vmaxget();

    T1 *buffer1 = (T1 *) R_alloc(maxNoRows, sizeof(T1));
    T2 *buffer2 = (T2 *) R_alloc(maxNoRows, sizeof(T2));

    for (i = 0; i < noCols; i++)
    {
        if (startIndex == NULL)
        {
            j = i * maxNoRows;

            while(array1[j] != maxUnSignedIndex && j < (i + 1)*maxNoRows)
                j++;

            j--;

            mergesort2(&array1[i * maxNoRows], &array2[i * maxNoRows], 0, j - i * maxNoRows, buffer1, buffer2);
        }
        else
        {
            mergesort2(&array1[startIndex[i]], &array2[startIndex[i]], 0, startIndex[i + 1] - startIndex[i] - 1,
                       buffer1, buffer2);
        }

    }

    vmaxset(vmax);
}

template<typename T>
void computeKernelMatrix(T maxUnSignedIndex, T *featVectorIndex, int32_t *featVectorValue,
                         Rcpp::NumericMatrix km, double *normValues, int fDim, int sizeX,
                         int sizeY, bool normalized, bool symmetric)
{
    uint32_t endx, endy, ix, iy;
    int i, j;
    double kv;

    // calculate kernel matrix
    if (symmetric)
    {
        // calculate kernel matrix
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            if (normalized)
            {
                if (normValues[i] == 0)
                    km(i, i) = 0;
                else
                    km(i, i) = 1;
            }
            else
                km(i, i) = normValues[i];

            for (j = i + 1; j < sizeY; j++)
            {
                ix = i * fDim;
                iy = j * fDim;
                endx = ix + fDim;
                endy = iy + fDim;
                kv = 0;

                while (ix < endx && iy < endy)
                {
                    if (featVectorIndex[ix] == maxUnSignedIndex ||
                        featVectorIndex[iy] == maxUnSignedIndex)
                    {
                        break;
                    }

                    if (featVectorIndex[ix] < featVectorIndex[iy])
                        ix++;
                    else if (featVectorIndex[ix] > featVectorIndex[iy])
                        iy++;
                    else
                    {
                        kv += featVectorValue[ix] * featVectorValue[iy];
                        ix++;
                        iy++;
                    }
                }

                if (normalized)
                {
                    if (kv == 0)
                        km(i, j) = 0;
                    else
                        km(i, j) = kv / normValues[i] / normValues[j];
                }
                else
                    km(i, j) = kv;

                km(j, i) = km(i, j);
            }
        }
    }
    else
    {
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            for (j = 0; j < sizeY; j++)
            {
                ix = i * fDim;
                iy = (sizeX + j) * fDim;
                endx = ix + fDim;
                endy = iy + fDim;
                kv = 0;

                while (ix < endx && iy < endy)
                {
                    if (featVectorIndex[ix] == maxUnSignedIndex ||
                        featVectorIndex[iy] == maxUnSignedIndex)
                        break;

                    if (featVectorIndex[ix] < featVectorIndex[iy])
                        ix++;
                    else if (featVectorIndex[ix] > featVectorIndex[iy])
                        iy++;
                    else
                    {
                        kv += featVectorValue[ix] * featVectorValue[iy];
                        ix++;
                        iy++;
                    }
                }

                if (normalized)
                {
                    if (kv == 0)
                        km(i, j) = 0;
                    else
                        km(i, j) = kv / normValues[i] / normValues[sizeX + j];
                }
                else
                    km(i, j) = kv;
            }
        }
    }
}

template<typename T>
void computeKernelMatrixPos(T maxUnSignedIndex, T *featVectorIndex, int32_t *featVectorValue,
                            Rcpp::NumericMatrix km, double *normValues, uint32_t fDim, int maxNumPatterns,
                            int sizeX, int sizeY, bool normalized, bool symmetric, bool computePosition,
                            Rcpp::NumericVector distWeight)
{
    uint32_t endx, endy, ix, iy, prevIndex;
    int i, j, j1, j2, posX, posY, distWeightLength, startIndex, numSamples, yOffset, numPatterns1, numPatterns2;
    double kv, distance;

    if (symmetric)
    {
        numSamples = sizeX;
        yOffset = 0;
    }
    else
    {
        numSamples = sizeX + sizeY;
        yOffset = sizeX;
        startIndex = 0;
    }


    if (distWeight.length() == 0)
    {
        // calculate kernel matrix
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            if (symmetric)
            {
                if (normalized)
                    km(i, i) = 1;
                else
                    km(i, i) = normValues[i];

                startIndex = i + 1;
            }

            if (computePosition)
                posX = 1 - featVectorValue[i];

            for (j = startIndex; j < sizeY; j++)
            {
                ix = i * fDim;
                iy = (yOffset + j) * fDim;
                endx = ix + fDim;
                endy = iy + fDim;
                prevIndex = 0;
                kv = 0;

                if (computePosition)
                {
                    posY = 1 - featVectorValue[yOffset + j];

                    if (posX > posY)
                        iy += posX - posY;
                    else
                        ix += posY - posX;

                    while (ix < endx && iy < endy)
                    {
                        if (featVectorIndex[ix] == maxUnSignedIndex ||
                            featVectorIndex[iy] == maxUnSignedIndex)
                            break;

                        if (featVectorIndex[ix] == featVectorIndex[iy])
                            kv++;

                        ix++;
                        iy++;
                    }
                }
                else
                {
                    while (ix < endx && iy < endy)
                    {
                        if (featVectorIndex[ix] == maxUnSignedIndex ||
                            featVectorIndex[iy] == maxUnSignedIndex)
                            break;

                        if (featVectorValue[ix] < featVectorValue[iy])
                            ix++;
                        else if (featVectorValue[ix] > featVectorValue[iy])
                            iy++;
                        else
                        {
                            if (maxNumPatterns == 1)
                            {
                                if (featVectorIndex[ix] == featVectorIndex[iy])
                                    kv++;

                                ix++;
                                iy++;
                            }
                            else
                            {
                                numPatterns1 = 0;
                                numPatterns2 = 0;

                                for (j1 = 0; j1 < maxNumPatterns; j1++)
                                {
                                    if (featVectorIndex[ix + j1] == maxUnSignedIndex)
                                        break;

                                    if (featVectorValue[ix + j1] != featVectorValue[ix])
                                        break;

                                    numPatterns1++;

                                    for (j2 = 0; j2 < maxNumPatterns; j2++)
                                    {
                                        if (featVectorIndex[ix + j2] == maxUnSignedIndex)
                                            break;

                                        if (featVectorValue[iy + j2] != featVectorValue[ix])
                                            break;

                                        if (j1 == 0)
                                            numPatterns2++;

                                        if (featVectorIndex[ix + j1] == featVectorIndex[iy + j2])
                                            kv++;
                                    }
                                }

                                ix += numPatterns1;
                                iy += numPatterns2;
                            }
                        }
                    }
                }

                if (normalized)
                {
                    if (kv == 0)
                        km(i, j) = 0;
                    else
                        km(i, j) = kv / normValues[i] / normValues[j + yOffset];
                }
                else
                    km(i, j) = kv;

                if (symmetric)
                    km(j, i) = km(i, j);
            }
        }
    }
    else // distance weighting
    {
        // sort feature vectors
        sort2Arrays(maxUnSignedIndex, featVectorIndex, featVectorValue, numSamples, fDim, NULL);

        distWeightLength = distWeight.size();

        // get kv of samples
        for (i = 0; i < numSamples; i++)
        {
            R_CheckUserInterrupt();

            ix = i * fDim;
            iy = i * fDim;
            prevIndex = 0;
            endx = ix + fDim;
            kv = 0;

            while (ix < endx && iy < endx)
            {
                if (featVectorIndex[ix] == maxUnSignedIndex ||
                    featVectorIndex[iy] == maxUnSignedIndex)
                    break;

                if (featVectorIndex[ix] < featVectorIndex[iy])
                    ix++;
                else if (featVectorIndex[ix] > featVectorIndex[iy])
                    iy++;
                else
                {
                    prevIndex = iy;

                    while (featVectorIndex[ix] == featVectorIndex[iy])
                    {
                        distance = abs(featVectorValue[iy++] - featVectorValue[ix]);

                        if (distance < distWeightLength)
                            kv = kv + distWeight[distance];
                    }

                    iy = prevIndex;
                    ix++;
                }
            }

            if (normalized)
            {
                normValues[i] = sqrt(kv);

                if (symmetric)
                {
                    if (kv == 0)
                        km(i,i) = 0;
                    else
                        km(i, i) = 1;
                }
            }
            else
            {
                if (symmetric)
                    km(i, i) = kv;
            }
        }

        // calculate kernel matrix
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            if (symmetric)
                startIndex = i + 1;

            for (j = startIndex; j < sizeY; j++)
            {
                ix = i * fDim;
                endx = ix + fDim;
                iy = (yOffset + j) * fDim;
                endy = (yOffset + j + 1) * fDim;
                prevIndex = 0;
                kv = 0;

                while (ix < endx && iy < endy)
                {
                    if (featVectorIndex[ix] == maxUnSignedIndex ||
                        featVectorIndex[iy] == maxUnSignedIndex)
                        break;

                    if (featVectorIndex[ix] < featVectorIndex[iy])
                        ix++;
                    else if (featVectorIndex[ix] > featVectorIndex[iy])
                        iy++;
                    else
                    {
                        prevIndex = iy;

                        while (featVectorIndex[ix] == featVectorIndex[iy])
                        {
                            distance = abs(featVectorValue[iy++] - featVectorValue[ix]);

                            if (distance < distWeightLength)
                                kv = kv + distWeight[distance];
                        }

                        iy = prevIndex;
                        ix++;
                    }
                }

                if (normalized)
                {
                    if (kv == 0)
                        km(i, j) = 0;
                    else
                        km(i, j) = kv / normValues[i] / normValues[j + yOffset];
                }
                else
                    km(i, j) = kv;

                if (symmetric)
                    km(j, i) = km(i, j);
            }
        }
    }
}

#endif
