//
// C Routines for Mismatch Kernel
//
// Source : Mismatch.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014   J o h a n n e s  P a l m e
//

#include "Kebabs.h"
#include "KernelUtils.h"

extern "C"
{
    #include "BitArray.h"
    #include "khash.h"
    #include "ksort.h"
}

#define LOWER_LONG            0x00000000FFFFFFFF
#define USER_INTERRUPT_LIMIT  100000
#define INDEX_MAP_INIT_FF     0xFF

using namespace Rcpp;

struct hmData
{
    double featWeight;
    uint32_t unweightedPosIndex;
};

#if __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

// hash map for mapping 64 bit feature indices to feature weights
KHASH_MAP_INIT_INT64(fw, struct hmData)
// hash map for mapping 64 bit feature indices
KHASH_MAP_INIT_INT64(fim, uint32_t)
// hash set for 64 bit feature indices
KHASH_MAP_INIT_INT64(ind, uint8_t)
// hash map for feature counting in case of normalization
KHASH_MAP_INIT_INT64(fc, uint32_t)
// sorting of 64 bit feature indices
KSORT_INIT(mism, uint64_t, ks_lt_generic)

struct ifMutateFeature {
    int sampleIndex;
    int k;
    int m;
    int numNonzeroFeatures;
    bool normalized;
    bool presence;
    bool populatedHashmap;
    bool storeFeatures;
    bool countNonzeroFeatures;
    bool calcKernelValue;
    uint32_t *featIndexMap;
    uint32_t *featureCounts;
    khash_t(fc) *fchmap;
    khash_t(fim) *hmap;
    uint64_t *powAlpha;
    uint64_t *featuresInSample;
    double kernelValue;
    NumericMatrix *pErd;
};

static uint32_t *pfeatIndexMap;
static uint32_t *pIndexMap;
static khash_t(fim) *pFeatureHMap;
static khash_t(ind) *pIndexHMap;
static khash_t(fw) *hmap;
static khash_t(fc) *pFeatureCountsHMap;
static uint32_t *pFeatureCounts;
static double *pNormValues;

#if __clang__
#pragma clang diagnostic pop
#endif


void buildSubtree(const char * x, int pos, int index, int curr, int k, int m,
                  int level, int noMism, double *sum, struct prefTree *pTree,
                  int maxNoOfNodes, int *freeNode, bool unmapped, bool presence,
                  bool *printWarning, struct indexBlock *nullBlock,
                  struct alphaInfo *alphaInf)
{
    int i, misMatchNew, indexNew, currNew;

    if (noMism == m)
    {
        if (pTree->node[curr].ib.idx[index] > 0)
        {
            curr = pTree->node[curr].ib.idx[index];
            if (level == k - 1)
            {
                if (pTree->node[curr].leaf == TRUE)
                {
                    if (!presence)
                    {
                        *sum = *sum - (double)pow(pTree->node[curr].value,2)
                                    + (double)pow(pTree->node[curr].value + 1,2);
                        pTree->node[curr].value++;
                    }
                }
                else
                {
                    if (*printWarning)
                    {
                        Rprintf("Invalid leaf reached:\n");
                        Rprintf("    curr: %d, index: %d, pos: %d\n",
                                curr, index, pos);
                    }
                }
                return;
            }
        }
        else
        {
            pTree->node[curr].ib.idx[index] = *freeNode;
            curr = *freeNode;
            if (curr == maxNoOfNodes)
            {
                if (*printWarning)
                {
                    Rprintf("Maximum number of nodes exceeded\n");
                    *printWarning = FALSE;
                }
                return;
            }

            *freeNode = *freeNode + 1;

            if (level == k-1)
            {
                pTree->node[curr].leaf = TRUE;
                pTree->node[curr].value = 1;
                *sum = *sum + 1;
                return;
            }
            else
            {
                pTree->node[curr].ib = *nullBlock;
                pTree->node[curr].leaf = FALSE;
            }
        }

        pos++;
        indexNew = alphaInf->seqIndexMap[(int) x[pos]];

        if (indexNew > -1)
        {
            buildSubtree(x, pos, indexNew, curr, k, m, level+1, noMism, sum, pTree, maxNoOfNodes,
                         freeNode, unmapped, presence, printWarning, nullBlock, alphaInf);
        }

        return;
    }
    else
    {
        for (i = 0; i < alphaInf->numAlphabetChars; i++)
        {
            if (pTree->node[curr].ib.idx[i] > 0)
            {
                currNew = pTree->node[curr].ib.idx[i];
                if (level == k - 1)
                {
                    if (pTree->node[currNew].leaf == TRUE)
                    {
                        if (!presence)
                        {
                            *sum = *sum - (double)pow(pTree->node[currNew].value,2)
                                        + (double)pow(pTree->node[currNew].value + 1,2);
                            pTree->node[currNew].value++;
                        }
                    }
                    else
                    {
                        if (*printWarning)
                        {
                            Rprintf("Invalid leaf reached:\n");
                            Rprintf("    curr: %d, index: %d, pos: %d\n",
                                    currNew, index, pos);
                        }
                    }
                    continue;
                }
            }
            else
            {
                pTree->node[curr].ib.idx[i] = *freeNode;
                currNew = *freeNode;
                if (currNew == maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }
                    return;
                }

                *freeNode = *freeNode + 1;
                if (level == k - 1)
                {
                    pTree->node[currNew].leaf = TRUE;
                    pTree->node[currNew].value = 1;
                    *sum = *sum + 1;
                    continue;
                }
                else
                {
                    pTree->node[currNew].ib = *nullBlock;
                    pTree->node[currNew].leaf = FALSE;
                }
            }

            indexNew = alphaInf->seqIndexMap[(int) x[pos + 1]];

            if (indexNew > -1)
            {
                misMatchNew = noMism;

                if (index != i)
                    misMatchNew++;

                buildSubtree(x, pos+1, indexNew, currNew, k, m, level+1, misMatchNew, sum, pTree, maxNoOfNodes,
                             freeNode, unmapped, presence, printWarning, nullBlock, alphaInf);
            }
        }
    }

    return;
}

double createMismatchTree(const char *s, int slen, int k, int m, struct prefTree *pTree,
                          int maxNoOfNodes, int *freeNode, bool unmapped, bool presence,
                          bool *printWarning, struct indexBlock *nullBlock,
                          struct alphaInfo *alphaInf)
{
    double kernelValue = 0;
    int curr, index;
    int len = slen - k;

    // init first node
    pTree->node[0].ib = *nullBlock;
    pTree->node[0].leaf = FALSE;

    // build mismatch tree for kmers of sequence
    for (int i=0; i<=len; i++)
    {
        curr = 0;

        index = alphaInf->seqIndexMap[(int) s[i]];

        if (index > -1)
        {
            buildSubtree(s, i, index, curr, k, m, 0, 0, &kernelValue, pTree, maxNoOfNodes, freeNode,
                         unmapped, presence, printWarning, nullBlock, alphaInf);
        }
    }

    return(kernelValue);
}

void traverseSubtree(const char * x, int length, int index, int pos, int curr, int k, int m,
                     int level, int noMism, double *kernelValue, struct prefTree *pTree,
                     int maxNoOfNodes, int *freeNode, bool unmapped, bool presence,
                     bool *printWarning, struct alphaInfo *alphaInf)
{
    int i, misMatchNew, indexNew, currNew;

    if (noMism == m)
    {
        for (i=level; i<k; i++)
        {
            if (pTree->node[curr].ib.idx[index] > 0)
            {
                currNew = pTree->node[curr].ib.idx[index];

                if (i == k - 1)
                {
                    if (pTree->node[currNew].leaf == TRUE)
                        *kernelValue += pTree->node[currNew].value;
                    else
                    {
                        if (*printWarning)
                        {
                            Rprintf("Invalid leaf reached:\n");
                            Rprintf("    curr: %d, index: %d, pos:%d\n",
                                    curr, index, pos);
                        }
                    }
                    return;
                }
                else
                {
                    index = alphaInf->seqIndexMap[(int) x[++pos]];

                    if (index < 0)
                        return;

                    curr = currNew;
                }
            }
            else
                return;
        }
        return;
    }
    else
    {
        for (i = 0; i < alphaInf->numAlphabetChars; i++)
        {
            if (pTree->node[curr].ib.idx[i] > 0)
            {
                currNew = pTree->node[curr].ib.idx[i];
                
                if (level == k - 1)
                {
                    if (pTree->node[currNew].leaf == TRUE)
                        *kernelValue += pTree->node[currNew].value;
                    else
                    {
                        if (*printWarning)
                        {
                            Rprintf("Invalid leaf reached:\n");
                            Rprintf("    curr: %d, index: %d, pos: %d\n",
                                    currNew, index, pos);
                        }
                    }
                    continue;
                }

                if (pos < (length - 1))
                {
                    indexNew = alphaInf->seqIndexMap[(int) x[pos + 1]];

                    if (indexNew > -1)
                    {
                        misMatchNew = noMism;

                        if (index != i)
                            misMatchNew++;

                        traverseSubtree(x, length, indexNew, pos+1, currNew, k, m, level+1, misMatchNew,
                                        kernelValue, pTree, maxNoOfNodes, freeNode, unmapped, presence,
                                        printWarning, alphaInf);
                    }
                }
            }
        }
    }

    return;
}

RcppExport SEXP getMismatchKernelMatrix(NumericMatrix km, ByteStringVector x, ByteStringVector y, int sizeX,
                                        int sizeY, IntegerVector selX, IntegerVector selY, bool symmetric,
                                        int bioCharset, bool lowercase, bool unmapped, int k, int m,
                                        bool normalized, bool presence, int maxSeqLength,
                                        struct alphaInfo *alphaInf)
{
    int size1, size2, i, j, p, iX, iY, freeNode, index, maxNoOfNodes;
    uint64_t nodeLimit;
    double kernelVal, currVal, currValSqrt;
    bool printWarning = TRUE;
    bool reversed = FALSE;
    struct prefTree *pTree;
    struct indexBlock nullBlock;

    NumericVector kvd(0);
    NumericVector& kv = kvd;

    // setup init mask
    for (i=0; i< MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    // alloc mem for prefix tree
    // consider mismatch subtree
    nodeLimit = (pow(alphaInf->numAlphabetChars, k + 1) - 1) / (alphaInf->numAlphabetChars - 1);

    if (nodeLimit < (uint64_t) MAX_BLOCK)
        maxNoOfNodes = (int) nodeLimit;
    else
        maxNoOfNodes = MAX_BLOCK;

    pTree = (struct prefTree *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    if (symmetric)
    {
        for (i=0; i<sizeX; i++)
        {
            R_CheckUserInterrupt();

            iX = selX[i];
            freeNode = 1;

            km(i,i) = createMismatchTree((const char *) x.ptr[iX], x.nchar[iX], k, m, pTree, maxNoOfNodes,
                                         &freeNode, unmapped, presence, &printWarning, &nullBlock, alphaInf);

            if (km(i,i) == NA_REAL)
                return(createNAMatrix(sizeX, sizeY));

            if (normalized == TRUE)
            {
                currValSqrt = sqrt(km(i,i));
                km(i,i) = 1;
            }

            for (j=i+1; j<sizeX; j++)
            {
                iY = selX[j];
                kernelVal = 0;

                // search second sequence
                for (p=0; p<=x.nchar[iY] - k; p++)
                {
                    index = alphaInf->seqIndexMap[(int) x.ptr[iY][p]];
                    
                    if (index > -1)
                    {
                        traverseSubtree((const char *)&(x.ptr[iY])[p], x.nchar[iY] - p, index, 0, 0,
                                        k, m, 0, 0, &kernelVal, pTree, maxNoOfNodes, &freeNode,
                                        unmapped, presence, &printWarning, alphaInf);
                    }
                }

                if (normalized == TRUE)
                    km(i,j) = kernelVal / currValSqrt;
                else
                    km(i,j) = kernelVal;

                km(j,i) = km(i,j);
            }

            if (normalized == TRUE)
            {
                for (int v=0; v<i; v++)
                {
                    km(v,i) = km(v,i) / currValSqrt;
                    km(i,v) = km(v,i);
                }
            }
        }
    }
    else // not symmetric kernel matrix
    {
        ByteStringVector *seqs1 = &x;
        ByteStringVector *seqs2 = &y;
        size1 = sizeX;
        size2 = sizeY;

        // precompute kernel values for x or y sequences
        if (normalized == TRUE)
        {
            if (sizeX < sizeY)
            {
                size1 = sizeY;
                size2 = sizeX;
                seqs1 = &y;
                seqs2 = &x;
                reversed = TRUE;
            }

            NumericVector kvd(size2);
            kv = kvd;

            for (i=0; i<size2; i++)
            {
                if (reversed)
                    iX = selY[i];
                else
                    iX = selX[i];

                freeNode = 1;

                kv[i] = createMismatchTree((char *) seqs2->ptr[iX], seqs2->nchar[iX], k, m, pTree,
                                           maxNoOfNodes, &freeNode, unmapped, presence, &printWarning,
                                           &nullBlock, alphaInf);

                if (kv[i] != NA_REAL)
                    kv[i] = sqrt(kv[i]);
            }
        }

        for (i=0; i<size1; i++)
        {
            if (reversed)
                iX = selY[i];
            else
                iX = selX[i];

            freeNode = 1;

            //create new tree for shorter dimension element
            currVal = createMismatchTree((char *) seqs1->ptr[iX], seqs1->nchar[iX], k, m, pTree,
                                         maxNoOfNodes, &freeNode, unmapped, presence, &printWarning,
                                         &nullBlock, alphaInf);

            if (currVal == NA_REAL)
            {
                for (j=0; j<size2; j++)
                {
                    if (reversed)
                        km(j,i) = NA_REAL;
                    else
                        km(i,j) = NA_REAL;
                }
                break;
            }

            currValSqrt = sqrt(currVal);

            for (j=0; j<size2; j++)
            {
                if (reversed)
                    iY = selX[j];
                else
                    iY = selY[j];

                kernelVal = 0;

                // search second sequence
                for (int p=0; p<=seqs2->nchar[iY]-k; p++)
                {
                    index = alphaInf->seqIndexMap[(int) seqs2->ptr[iY][p]];

                    if (index > -1)
                    {
                        traverseSubtree(((const char *)&(seqs2->ptr[iY])[p]), seqs2->nchar[iY], index, 0, 0,
                                        k, m, 0, 0, &kernelVal, pTree, maxNoOfNodes, &freeNode, unmapped,
                                        presence, &printWarning, alphaInf);
                    }
                }

                if (normalized == TRUE)
                {
                    if (reversed)
                        km(j,i) = kernelVal / currValSqrt / kv[j];
                    else
                        km(i,j) = kernelVal / currValSqrt / kv[j];
                }
                else
                {
                    if (reversed)
                        km(j,i) = kernelVal;
                    else
                        km(i,j) = kernelVal;
                }
            }
        }
    }

    // dealloc of heap done by R after .Call and in error cases

    return(km);
}

void mutateFeatureIndexViaArray(uint64_t featureIndex, struct ifMutateFeature *intf)
{
    int i, numMismatches, pos, currStack, origIndex;
    uint32_t low, featureStored, stack[6*intf->k], currCount;
    uint64_t currFeatureIndex;
    bool finished;

    currStack = -1;
    finished = FALSE;
    currFeatureIndex = featureIndex;
    numMismatches = 0;
    pos = 0;
    i = 0;

    if (intf->countNonzeroFeatures)
        featureStored = intf->sampleIndex;
    else
        featureStored = 1;

    while (!finished)
    {
        if (numMismatches == intf->m || pos == intf->k)
        {
            if (intf->storeFeatures)
            {
                if (intf->featIndexMap != NULL)
                {
                    if (intf->featIndexMap[(uint32_t) currFeatureIndex] != MAXUINT32)
                    {
                        if (intf->presence)
                        {
                            if ((*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex]) != 1)
                            {
                                (*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex]) = 1;

                                if (intf->calcKernelValue)
                                    intf->kernelValue += 1;
                            }
                        }
                        else
                        {
                            (*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex])  =
                            (*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex])  + 1;

                            if (intf->calcKernelValue)
                            {
                                intf->kernelValue = intf->kernelValue -
                                ((*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex]) - 1) *
                                ((*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex])  - 1) +
                                (*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex])  *
                                (*(intf->pErd))(intf->sampleIndex, intf->featIndexMap[(uint32_t) currFeatureIndex]) ;
                            }
                        }
                    }
                }
                else
                {
                    if (intf->presence)
                    {
                        if ((*(intf->pErd))(intf->sampleIndex, currFeatureIndex) != 1)
                        {
                            (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) = 1;

                            if (intf->calcKernelValue)
                                intf->kernelValue += 1;
                        }
                    }
                    else
                    {
                        (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) =
                        (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) + 1;

                        if (intf->calcKernelValue)
                        {
                            intf->kernelValue = intf->kernelValue -
                            ((*(intf->pErd))(intf->sampleIndex, currFeatureIndex) - 1) *
                            ((*(intf->pErd))(intf->sampleIndex, currFeatureIndex) - 1) +
                            (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) *
                            (*(intf->pErd))(intf->sampleIndex, currFeatureIndex);
                        }
                    }
                }
            }
            else
            {
                if (intf->calcKernelValue)
                {
                    if (getBit(intf->featuresInSample, currFeatureIndex) != 0)
                    {
                        if (!intf->presence)
                        {
                            currCount = intf->featureCounts[currFeatureIndex];
                            intf->featureCounts[(uint32_t) currFeatureIndex] += 1;
                            intf->kernelValue = intf->kernelValue - currCount*currCount +
                                                (currCount + 1) * (currCount + 1);
                        }
                    }
                    else
                    {
                        setBit(intf->featuresInSample, currFeatureIndex);
                        intf->featureCounts[currFeatureIndex] = 1;
                        intf->kernelValue += 1;
                    }
                }

                // mark feature as used
                if (intf->featIndexMap[currFeatureIndex] < MAXUINT32 &&
                    intf->featIndexMap[currFeatureIndex] != featureStored)
                {
                    intf->featIndexMap[currFeatureIndex] = featureStored;

                    if (intf->countNonzeroFeatures)
                        intf->numNonzeroFeatures++;
                }
            }

            if (currStack > -1)
            {
                // retrieve previous state and continue
                origIndex = stack[currStack--];
                numMismatches = stack[currStack--];
                pos = stack[currStack--];
                i = stack[currStack--];
                low = stack[currStack--];
                currFeatureIndex = (((uint64_t) stack[currStack--]) << 32) + low;
            }
            else
                finished = TRUE;
        }

        if (i == 0)
        {
            origIndex = (currFeatureIndex / intf->powAlpha[pos]) % intf->powAlpha[1];
            currFeatureIndex = currFeatureIndex - (origIndex) * intf->powAlpha[pos];
        }

        if (i < (int) intf->powAlpha[1])
        {
            // push state to stack and descend
            stack[++currStack] = (uint32_t) ((currFeatureIndex + intf->powAlpha[pos]) >> 32);
            stack[++currStack] = (uint32_t) (currFeatureIndex + intf->powAlpha[pos]);
            stack[++currStack] = i + 1;
            stack[++currStack] = pos;
            stack[++currStack] = numMismatches;
            stack[++currStack] = origIndex;

            if (i != origIndex)
                numMismatches++;

            i = 0;
            pos++;
        }
        else
        {
            if (currStack > -1)
            {
                // retrieve previous state and continue
                origIndex = stack[currStack--];
                numMismatches = stack[currStack--];
                pos = stack[currStack--];
                i = stack[currStack--];
                low = stack[currStack--];
                currFeatureIndex = (((uint64_t) stack[currStack--]) << 32) + low;
            }
            else
                finished = TRUE;
        }
    }
}

void mutateFeatureIndex(uint64_t featureIndex, struct ifMutateFeature *intf)
{
    int currStack, result, numMismatches, pos;
    uint32_t i, stack[6*intf->k], low, featureStored, colIndex, currCount;
    uint64_t currFeatureIndex, origIndex;
    bool finished;
    khiter_t iter;

    currStack = -1;
    finished = FALSE;
    currFeatureIndex = featureIndex;
    numMismatches = 0;
    pos = 0;
    i = 0;

    if (intf->countNonzeroFeatures)
        featureStored = intf->sampleIndex;
    else
        featureStored = 1;

    while (!finished)
    {
        if (numMismatches == intf->m || pos == intf->k)
        {
            if (intf->pErd != NULL)
            {
                // count feature
                if (intf->hmap != NULL)
                {
                    iter = kh_get(fim, intf->hmap, currFeatureIndex);

                    if (iter == kh_end(intf->hmap))
                    {
                        Rprintf("Feature index %llu not found in index mapping\n", currFeatureIndex);
                        return;
                    }

                    colIndex = kh_value(intf->hmap, iter);

                    if (colIndex != MAXUINT32)
                    {
                        if (intf->presence)
                        {
                            if ((*(intf->pErd))(intf->sampleIndex, colIndex) != 1)
                            {
                                (*(intf->pErd))(intf->sampleIndex, colIndex) = 1;

                                if (intf->calcKernelValue)
                                    intf->kernelValue += 1;
                            }
                        }
                        else
                        {
                            (*(intf->pErd))(intf->sampleIndex, colIndex) =
                                (*(intf->pErd))(intf->sampleIndex, colIndex) + 1;

                            if (intf->calcKernelValue)
                            {
                                intf->kernelValue = intf->kernelValue -
                                    ((*(intf->pErd))(intf->sampleIndex, colIndex) - 1) *
                                    ((*(intf->pErd))(intf->sampleIndex, colIndex) - 1) +
                                    (*(intf->pErd))(intf->sampleIndex, colIndex) *
                                    (*(intf->pErd))(intf->sampleIndex, colIndex);
                            }
                        }
                    }
                }
                else
                {
                    if (intf->presence)
                    {
                        if ((*(intf->pErd))(intf->sampleIndex, currFeatureIndex) != 1)
                        {
                            (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) = 1;

                            if (intf->calcKernelValue)
                                intf->kernelValue += 1;
                        }
                    }
                    else
                    {
                        (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) =
                            (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) + 1;

                        if (intf->calcKernelValue)
                        {
                            intf->kernelValue = intf->kernelValue -
                                ((*(intf->pErd))(intf->sampleIndex, currFeatureIndex) - 1) *
                                ((*(intf->pErd))(intf->sampleIndex, currFeatureIndex) - 1) +
                                (*(intf->pErd))(intf->sampleIndex, currFeatureIndex) *
                                (*(intf->pErd))(intf->sampleIndex, currFeatureIndex);
                        }
                    }
                }
            }
            else
            {
                if (intf->calcKernelValue)
                {
                    iter = kh_get(fc, intf->fchmap, currFeatureIndex);

                    if (iter != kh_end(intf->fchmap))
                    {
                        currCount = kh_value(intf->fchmap, iter);
                        kh_value(intf->fchmap, iter) = currCount + 1;
                        intf->kernelValue = intf->kernelValue - currCount*currCount +
                        (currCount + 1) * (currCount + 1);
                    }
                    else
                    {
                        iter = kh_put(fc, intf->fchmap, currFeatureIndex, &result);

                        if (result == -1)
                        {
                            Rprintf("Storage of key %llu in feature count hashmap failed\n",
                                    currFeatureIndex);
                            return;
                        }

                        kh_value(intf->fchmap, iter) = 1;
                        intf->kernelValue += 1;
                    }
                }

                // mark feature as used
                if (intf->hmap != NULL)
                {
                    iter = kh_get(fim, intf->hmap, currFeatureIndex);

                    if (iter != kh_end(intf->hmap))
                    {
                        if (kh_value(intf->hmap, iter) != featureStored)
                        {
                            kh_value(intf->hmap, iter) = featureStored;

                            if (intf->countNonzeroFeatures)
                                intf->numNonzeroFeatures++;
                        }
                    }
                    else
                    {
                        if (!intf->populatedHashmap)
                        {
                            iter = kh_put(fim, intf->hmap, currFeatureIndex, &result);

                            if (result == -1)
                            {
                                Rprintf("Storage of key %llu in hashmap failed\n", currFeatureIndex);
                                return;
                            }

                            kh_value(intf->hmap, iter) = featureStored;

                            if (intf->countNonzeroFeatures)
                                intf->numNonzeroFeatures++;
                        }
                    }
                }
                else
                {
                    // $$$ TODO remove else branch
                    Rprintf("This branch should not be executed\n");
                    if (intf->countNonzeroFeatures)
                        intf->numNonzeroFeatures++;
                }
            }

            if (currStack > -1)
            {
                // retrieve previous state and continue
                origIndex = stack[currStack--];
                numMismatches = stack[currStack--];
                pos = stack[currStack--];
                i = stack[currStack--];
                low = stack[currStack--];
                currFeatureIndex = (((uint64_t) stack[currStack--]) << 32) + low;
            }
            else
                finished = TRUE;
        }

        if (i == 0)
        {
            origIndex = (currFeatureIndex / intf->powAlpha[pos]) % intf->powAlpha[1];
            currFeatureIndex = currFeatureIndex - (origIndex) * intf->powAlpha[pos];
        }

        if (i < intf->powAlpha[1])
        {
            // push state to stack and descend
            stack[++currStack] = (uint32_t) ((currFeatureIndex + intf->powAlpha[pos]) >> 32);
            stack[++currStack] = (uint32_t) ((currFeatureIndex + intf->powAlpha[pos]) & LOWER_LONG);
            stack[++currStack] = i + 1;
            stack[++currStack] = pos;
            stack[++currStack] = numMismatches;
            stack[++currStack] = origIndex;

            if (i != origIndex)
                numMismatches++;

            i = 0;
            pos++;
        }
        else
        {
            if (currStack > -1)
            {
                // retrieve previous state and continue
                origIndex = stack[currStack--];
                numMismatches = stack[currStack--];
                pos = stack[currStack--];
                i = stack[currStack--];
                low = stack[currStack--];
                currFeatureIndex = (((uint64_t) stack[currStack--]) << 32) + low;
            }
            else
                finished = TRUE;
        }
    }
}

bool getIndexMap(ByteStringVector x, int sizeX, IntegerVector selX, bool unmapped, int bioCharset,
                 int k, int m, bool normalized, bool presence, struct alphaInfo *alphaInf,
                 ByteStringVector features, uint64_t dimFeatureSpace, bool zeroFeatures, bool useHash,
                 void **indexMap, int *numUsedFeatures, bool countNonzeroFeatures, int *numNonzeroFeatures,
                 double **normValues)
{
    int i, j, iold, iX, index, patLength, result;
    uint32_t *featIndexMap, *featureCounts, *indexSet, featureTreated, maxValidIndex;
    uint64_t featureIndex, fIndex, prevIndex, *featuresInSample;
    bool calcKernelValue;
    const void *vmax;
    khiter_t iter;
    khash_t(fim) *hmap;
    khash_t(fc) *fchmap;
    khash_t(ind) *indset;
    struct ifMutateFeature intf;

    indset = NULL;
    indexSet = NULL;
    featureCounts = NULL;
    featuresInSample = NULL;
    vmax = NULL;

    // cyclic buffers for old index
    uint64_t *oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    uint64_t numAlphaPowK_1 = pow(alphaInf->numAlphabetChars, k - 1);
    uint64_t *powAlpha = (uint64_t *) R_alloc(k + 1, sizeof(uint64_t));

    for (i = 0; i <= k; i++)
        powAlpha[i] = pow(alphaInf->numAlphabetChars, i);

    if (normalized && !(zeroFeatures && features.length < 1))
    {
        calcKernelValue = TRUE;
        *normValues = (double *) Calloc(sizeX, double);
        pNormValues = *normValues;

        if (useHash)
        {
            fchmap = kh_init(fc);
            pFeatureCountsHMap = fchmap;
        }
        else
        {
            DEFINE_BITARRAY(usedFeatures, dimFeatureSpace);
            featuresInSample = usedFeatures;
            featureCounts = (uint32_t *) Calloc(dimFeatureSpace, uint32_t);
            pFeatureCounts = featureCounts;
        }
    }
    else
        calcKernelValue = FALSE;

    if (useHash)
    {
        hmap = kh_init(fim);    //alloc hash table for feature index mapping
        pFeatureHMap = hmap;
        *indexMap = (void *) hmap;
        featIndexMap = NULL;
        pfeatIndexMap = NULL;
        pIndexMap = NULL;

        if (!calcKernelValue)
        {
            indset = kh_init(ind);  //alloc hash table for index set
            pIndexHMap = indset;
        }
    }
    else
    {
        featIndexMap = (uint32_t *) Calloc(dimFeatureSpace, uint32_t);
        pfeatIndexMap = featIndexMap;
        *indexMap = (void *) featIndexMap;
        vmax = vmaxget();
        hmap = NULL;
        pFeatureHMap = NULL;
        pIndexHMap = NULL;

        if (!calcKernelValue)
        {
            indexSet = (uint32_t *) Calloc(dimFeatureSpace, uint32_t);
            pIndexMap = indexSet;
        }

        if (features.length < 1)
        {
            // no support for multibyte in memset because of endianess of multibyte integer
            for (i = 0; i < (int) dimFeatureSpace; i++)
                // init to MAXUINT32-1
                featIndexMap[i] = VALID_FEATURE;
        }
        else
            // init each element to MAXUINT32
            memset(featIndexMap, INDEX_MAP_INIT_FF, dimFeatureSpace*sizeof(uint32_t));

        R_CheckUserInterrupt();

        if (!calcKernelValue)
            memset(indexSet, INDEX_MAP_INIT_FF, dimFeatureSpace*sizeof(uint32_t));
    }

    intf.populatedHashmap = FALSE;

    if (features.length > 0)
    {
        intf.populatedHashmap = TRUE;

        // store requested features
        for (i=0; i < features.length; i++)
        {
            if (i % USER_INTERRUPT_LIMIT == 0)
                R_CheckUserInterrupt();

            featureIndex = 0;

            for (j=0; j < k; j++)
            {
                featureIndex = featureIndex * alphaInf->numAlphabetChars +
                               alphaInf->indexMap[(int) features.ptr[i][j]];
            }

            if (useHash)
            {
                iter = kh_get(fim, hmap, featureIndex);

                if (iter == kh_end(hmap))
                {
                    iter = kh_put(fim, hmap, featureIndex, &result);

                    if (result != -1)
                        kh_value(hmap, iter) = VALID_FEATURE;
                }
                else
                {
                    // $$$  TODO remove Rprintf
                    Rprintf("Features are not unique");
                    return(FALSE);
                }
            }
            else
                featIndexMap[featureIndex] = VALID_FEATURE;
        }
    }

    featureTreated = 1;

    intf.k = k;
    intf.m = m;
    intf.normalized = normalized;
    intf.presence = presence;
    intf.featIndexMap = featIndexMap;
    intf.hmap = hmap;
    intf.powAlpha = powAlpha;
    intf.numNonzeroFeatures = 0;
    intf.pErd = NULL;
    intf.storeFeatures = FALSE;
    intf.countNonzeroFeatures = countNonzeroFeatures;
    intf.calcKernelValue = calcKernelValue;
    intf.featureCounts = featureCounts;
    intf.fchmap = fchmap;
    intf.featuresInSample = featuresInSample;
    // intf.populatedHashmap is already assigned

    // determine used features
    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();

        if (normalized && !(zeroFeatures && features.length < 1))
        {
            if (useHash)
                kh_clear(fc, fchmap);
            else
                SET_BITARRAY(featuresInSample, dimFeatureSpace, 0);
        }

        featureIndex = 0;
        patLength = 0;
        iold = 0;
        iX = selX[i];
        intf.sampleIndex = i;
        intf.kernelValue = 0;

        if (numNonzeroFeatures != NULL)
            featureTreated = i;

        for (j = 0; j < x.nchar[iX]; j++)
        {
            index = alphaInf->seqIndexMap[(int) x.ptr[iX][j]];

            if (index > -1)
            {
                if (patLength < k)
                {
                    oldIndex[iold++] = index * numAlphaPowK_1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = featureIndex * alphaInf->numAlphabetChars + index;

                    patLength++;

                    if (patLength == k)
                    {
                        if (calcKernelValue)
                        {
                            if (useHash)
                                mutateFeatureIndex(featureIndex, &intf);
                            else
                                mutateFeatureIndexViaArray(featureIndex, &intf);
                        }
                        else
                        {
                            // mutate and store
                            if (useHash)
                            {
                                // check if new feature
                                iter = kh_get(ind, indset, featureIndex);

                                if (iter == kh_end(indset))
                                {
                                    iter = kh_put(ind, indset, featureIndex, &result);

                                    if (result != -1)
                                    {
                                        kh_value(indset, iter) = featureTreated;
                                        mutateFeatureIndex(featureIndex, &intf);
                                    }
                                    else
                                    {
                                        Rprintf("Feature index could not be stored\n");
                                        kh_destroy(ind, indset);
                                        pIndexHMap = NULL;
                                        return(false);
                                    }
                                }
                                else
                                {
                                    if (kh_value(indset, iter) != featureTreated)
                                    {
                                        kh_value(indset, iter) = featureTreated;
                                        mutateFeatureIndex(featureIndex, &intf);
                                    }
                                }
                            }
                            else
                            {
                                if (indexSet[featureIndex] != featureTreated)
                                {
                                    indexSet[featureIndex] = featureTreated;
                                    mutateFeatureIndexViaArray(featureIndex, &intf);
                                }
                            }
                        }
                    }
                }
                else
                {
                    prevIndex = oldIndex[iold];
                    oldIndex[iold++] = index * numAlphaPowK_1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = (featureIndex - prevIndex) *
                    alphaInf->numAlphabetChars + index;

                    if (calcKernelValue)
                    {
                        if (useHash)
                            mutateFeatureIndex(featureIndex, &intf);
                        else
                            mutateFeatureIndexViaArray(featureIndex, &intf);
                    }
                    else
                    {
                        // mutate and store
                        if (useHash)
                        {
                            // check if new feature
                            iter = kh_get(ind, indset, featureIndex);

                            if (iter == kh_end(indset))
                            {
                                iter = kh_put(ind, indset, featureIndex, &result);

                                if (result != -1)
                                {
                                    kh_value(indset, iter) = featureTreated;
                                    mutateFeatureIndex(featureIndex, &intf);
                                }
                                else
                                {
                                    Rprintf("Feature index could not be stored\n");
                                    kh_destroy(ind, indset);
                                    pIndexHMap = NULL;
                                    return(false);
                                }
                            }
                            else
                            {
                                if (kh_value(indset, iter) != featureTreated)
                                {
                                    kh_value(indset, iter) = featureTreated;
                                    mutateFeatureIndex(featureIndex, &intf);
                                }
                            }
                        }
                        else
                        {
                            if (indexSet[featureIndex] != featureTreated)
                            {
                                indexSet[featureIndex] = featureTreated;
                                mutateFeatureIndexViaArray(featureIndex, &intf);
                            }
                        }
                    }
                }
            }
            else // wrong character - restart kmer
            {
                featureIndex = 0;
                patLength = 0;
                iold = 0;
            }
        }

        if (calcKernelValue)
            (*normValues)[i] = sqrt(intf.kernelValue);

        // reset speed up array/hashmap for sparse ER after each sample
        if (countNonzeroFeatures && !calcKernelValue)
        {    if(useHash)
                kh_clear(ind, indset);
            else
                memset(indexSet, INDEX_MAP_INIT_FF, dimFeatureSpace*sizeof(uint32_t));
        }
    }

    if (!useHash)
        vmaxset(vmax);
    else
    {
        if (!calcKernelValue)
        {
            kh_destroy(ind, indset);
            pIndexHMap = NULL;
        }
    }

    if (countNonzeroFeatures)
        *numNonzeroFeatures = intf.numNonzeroFeatures;

    if (!zeroFeatures || features.length > 0)
    {
        if (!useHash)
        {
            fIndex = 0;

            if (zeroFeatures)
                maxValidIndex = MAXUINT32;
            else
                maxValidIndex = VALID_FEATURE;

            for (i=0; i < (int) dimFeatureSpace; i++)
            {
                if (i % USER_INTERRUPT_LIMIT == 0)
                    R_CheckUserInterrupt();

                if (featIndexMap[i] < maxValidIndex)
                    featIndexMap[i] = fIndex++;
                else
                {
                    if (featIndexMap[i] != MAXUINT32)
                        featIndexMap[i] = MAXUINT32;
                }
            }

            *numUsedFeatures = fIndex;
        }
        else
        {
            // sort entries
            int noOfEntries = kh_size(hmap);
            if (noOfEntries > 0)
            {
                uint64_t *featureIndices = (uint64_t *) malloc(noOfEntries * sizeof(uint64_t));

                i = 0;
                *numUsedFeatures = 0;

                for (iter = kh_begin(hmap); iter != kh_end(hmap); iter++)
                {
                    i++;

                    if (i % USER_INTERRUPT_LIMIT == 0)
                        R_CheckUserInterrupt();

                    if (kh_exist(hmap, iter))
                        featureIndices[(*numUsedFeatures)++] = kh_key(hmap, iter);
                }

                ks_mergesort(mism, noOfEntries, featureIndices, 0);

                // write col index to hash
                for (int i=0; i < *numUsedFeatures; i++)
                {
                    if (i % USER_INTERRUPT_LIMIT == 0)
                        R_CheckUserInterrupt();

                    iter = kh_get(fim, hmap, featureIndices[i]);

                    if (kh_exist(hmap, iter))
                    {
                        if (kh_value(hmap, iter) < VALID_FEATURE)
                            kh_value(hmap, iter) = i;
                    }
                    else
                    {
                        // $$$  TODO remove Rprintf
                        Rprintf("Internal error with hashing - entry not found\n");
                        return(false);
                    }
                }
            }
        }
    }

    return(true);
}

void getERDMismatch(ByteStringVector x, int sizeX, IntegerVector selX, bool unmapped, int bioCharset,
                    int k, int m, bool normalized, bool presence, struct alphaInfo *alphaInf,
                    ByteStringVector features, uint64_t dimFeatureSpace, bool zeroFeatures, bool useHash,
                    void *indexMap, int numUsedFeatures, NumericMatrix erd, double *normValues)
{
    int i, j, iold, iX, index, patLength;
    uint32_t *featIndexMap;
    uint64_t featureIndex, prevIndex;
    bool calcKernelValue;
    khash_t(fim) *hmap;
    struct ifMutateFeature intf;

    if (normalized && normValues == NULL)
    {
        normValues = (double *) Calloc(sizeX, double);
        pNormValues = normValues;
        calcKernelValue = TRUE;
    }
    else
        calcKernelValue = FALSE;

    if (useHash)
    {
        hmap = (khash_t(fim) *) indexMap;
        featIndexMap = NULL;
    }
    else
    {
        featIndexMap = (uint32_t *) indexMap;
        hmap = NULL;
    }

    // cyclic buffers for old index
    uint64_t *oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    uint64_t numAlphaPowK_1 = pow(alphaInf->numAlphabetChars, k - 1);
    uint64_t *powAlpha = (uint64_t *) R_alloc(k + 1, sizeof(uint64_t));

    for (i = 0; i <= k; i++)
        powAlpha[i] = pow(alphaInf->numAlphabetChars, i);

    intf.k = k;
    intf.m = m;
    intf.normalized = normalized;
    intf.presence = presence;
    intf.featIndexMap = featIndexMap;
    intf.hmap = hmap;
    intf.powAlpha = powAlpha;
    intf.numNonzeroFeatures = 0;
    intf.countNonzeroFeatures = FALSE;
    intf.populatedHashmap = FALSE;
    intf.storeFeatures = TRUE;
    intf.calcKernelValue = calcKernelValue;
    intf.pErd = &erd;

    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();

        featureIndex = 0;
        patLength = 0;
        iold = 0;
        iX = selX[i];
        intf.sampleIndex = i;
        intf.kernelValue = 0;

        for (j = 0; j < x.nchar[iX]; j++)
        {
            index = alphaInf->seqIndexMap[(int) x.ptr[iX][j]];

            if (index > -1)
            {
                if (patLength < k)
                {
                    oldIndex[iold++] = index * numAlphaPowK_1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = featureIndex * alphaInf->numAlphabetChars + index;

                    patLength++;

                    if (patLength == k)
                    {
                        // mutate and store
                        if (useHash)
                            mutateFeatureIndex(featureIndex, &intf);
                        else
                            mutateFeatureIndexViaArray(featureIndex, &intf);
                    }
                }
                else
                {
                    prevIndex = oldIndex[iold];
                    oldIndex[iold++] = index * numAlphaPowK_1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = (featureIndex - prevIndex) *
                                   alphaInf->numAlphabetChars + index;

                    if (useHash)
                        mutateFeatureIndex(featureIndex, &intf);
                    else
                        mutateFeatureIndexViaArray(featureIndex, &intf);
                }
            }
            else // wrong character - restart kmer
            {
                featureIndex = 0;
                patLength = 0;
                iold = 0;
            }
        }

        if (calcKernelValue)
            normValues[i] = sqrt(intf.kernelValue);
    }

    if (normalized)
    {
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            if (normValues[i] == 0)
                continue;

            for (j = 0; j < numUsedFeatures; j++)
            {
                if (erd(i,j) > 0)
                    erd(i,j) = erd(i,j) / normValues[i];
            }
        }
    }

    // memory release through mismatchFreeHeapC called in on.exit from R level

    return;
}

void storeFeatures(struct prefTree *pTree, SEXP slot_p, SEXP slot_j, SEXP slot_x, int *jIndex, int k,
                   bool normalized, double normValue, bool zeroFeatures, bool mapIndex, bool useHash,
                   void *indexMap, struct alphaInfo *alphaInf)
{
    int currBlock, currIndex, newBlock, maxBlockIndex, currStack, stack[2*k];
    uint32_t irrelevantFeature, *featIndexMap;
    uint64_t featureIndex, fIndex;
    bool skip;
    khash_t(fim) *hmap;
    khiter_t iter;

    if (useHash)
    {
        hmap = (khash_t(fim) *) indexMap;
        featIndexMap = NULL;
    }
    else
    {
        featIndexMap = (uint32_t *) indexMap;
        hmap = NULL;
    }

    // walk tree to store sparse feature vector to ers
    maxBlockIndex = alphaInf->maxAlphaIndex;
    currBlock = 0;
    currIndex = 0;
    currStack = -1;
    featureIndex = 0;
    skip = FALSE;

    if (zeroFeatures)
        irrelevantFeature = MAXUINT32;       // all valid
    else
        irrelevantFeature = VALID_FEATURE;   // used only

    while (currStack >= 0 || currIndex <= maxBlockIndex)
    {
        if (pTree->node[currBlock].ib.idx[currIndex] > 0)
        {
            newBlock = pTree->node[currBlock].ib.idx[currIndex];

            if (pTree->node[newBlock].leaf)
            {
                fIndex = featureIndex * alphaInf->numAlphabetChars + currIndex;

                if (mapIndex)
                {
                    if (useHash)
                    {
                        iter = kh_get(fim, hmap, featureIndex);

                        if (iter != kh_end(hmap))
                        {
                            if (kh_value(hmap, iter) < irrelevantFeature)
                            {
                                INTEGER(slot_j)[*jIndex] = kh_value(hmap, iter);
                            }
                            else
                                skip = TRUE;
                        }
                        else
                            skip = TRUE;
                    }
                    else
                    {
                        if (featIndexMap[fIndex] < irrelevantFeature)
                            INTEGER(slot_j)[*jIndex] = featIndexMap[fIndex];
                        else
                            skip = TRUE;
                    }
                }
                else
                    INTEGER(slot_j)[*jIndex] = fIndex;

                if (!skip)
                {
                    if (normalized)
                    {
                        if (pTree->node[newBlock].value == 0)
                            REAL(slot_x)[*jIndex] = 0;
                        else
                            REAL(slot_x)[*jIndex] = pTree->node[newBlock].value / normValue;
                    }
                    else
                        REAL(slot_x)[*jIndex] = pTree->node[newBlock].value;

                    (*jIndex)++;
                }
                else
                    skip = FALSE;

                currIndex++;

                if (currIndex > maxBlockIndex)
                {
                    if (currStack == -1)
                        continue;

                    while (currStack >= 0 && currIndex > maxBlockIndex)
                    {
                        currIndex = stack[currStack--];
                        currBlock = stack[currStack--];
                        featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                    }
                }
            }
            else
            {
                stack[++currStack] = currBlock;
                stack[++currStack] = currIndex + 1;
                featureIndex = featureIndex * alphaInf->numAlphabetChars + currIndex;
                currBlock = newBlock;
                currIndex = 0;
            }
        }
        else
        {
            currIndex++;

            if (currIndex > maxBlockIndex)
            {
                if (currStack == -1)
                    continue;

                while (currStack >= 0 && currIndex > maxBlockIndex)
                {
                    currIndex = stack[currStack--];
                    currBlock = stack[currStack--];
                    featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                }
            }
        }
    }
}

bool getERSMismatch(ByteStringVector x, int sizeX, IntegerVector selX, int maxSeqLength, bool unmapped,
                    int bioCharset, int k, int m, bool normalized, bool presence, struct alphaInfo *alphaInf,
                    ByteStringVector features, uint64_t dimFeatureSpace, bool zeroFeatures, bool mapIndex,
                    bool useHash, void *indexMap, int numUsedFeatures, SEXP slot_p, SEXP slot_j, SEXP slot_x,
                    double *normValues)
{
    int i, iX, freeNode, maxNoOfNodes, jIndex;
    uint64_t nodeLimit, maxNodesPerSequence;
    double kernelValue;
    bool printWarning = TRUE;
    struct prefTree *pTree;
    struct indexBlock nullBlock;

    // alloc mem for prefix tree
    // consider mismatch subtree
    nodeLimit = ((pow(alphaInf->numAlphabetChars, k + 1) - 1) / (alphaInf->numAlphabetChars - 1)) *
                pow(alphaInf->numAlphabetChars, k) + 1;

    maxNodesPerSequence = pow(alphaInf->numAlphabetChars, m) * k * (maxSeqLength - k + 1) + 1;

    if (maxNodesPerSequence < (uint64_t) nodeLimit)
        nodeLimit = (int) maxNodesPerSequence;

    maxNoOfNodes = MAX_BLOCK;

    if (nodeLimit < (uint64_t) MAX_BLOCK)
        maxNoOfNodes = (int) nodeLimit;

    pTree = (struct prefTree *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    // setup init mask
    for (i=0; i< MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    jIndex = 0;

    for (i=0; i<sizeX; i++)
    {
        R_CheckUserInterrupt();

        freeNode = 1;
        iX = selX[i];

        kernelValue = createMismatchTree((const char *) x.ptr[iX], x.nchar[iX], k, m, pTree, maxNoOfNodes,
                                         &freeNode, unmapped, presence, &printWarning, &nullBlock, alphaInf);

        if (kernelValue == NA_REAL)
        {
            Rprintf("Mismatch tree could not be created for sample %d\n", selX[i]);
            return(FALSE);
        }

        if (normalized)
            kernelValue = sqrt(kernelValue);

        INTEGER(slot_p)[i] = jIndex;

        storeFeatures(pTree, slot_p, slot_j, slot_x, &jIndex, k, normalized, kernelValue, zeroFeatures,
                      mapIndex, useHash, indexMap, alphaInf);
    }

    INTEGER(slot_p)[sizeX] = jIndex;

    // memory release through mismatchFreeHeapC called in on.exit from R level

    return(TRUE);
}

void assignFeatureNames(SEXP colnames, void *indexMap, int k, struct alphaInfo *alphaInf, uint64_t dimFeatureSpace,
                        bool mapIndex, bool useHash)
{
    int i, j, l;
    uint32_t *featIndexMap;
    khash_t(fim) *hmap;
    khiter_t iter;
    char kmer[k + 1]; // add char for \0

    kmer[k] = '\0';
    uint64_t *powAlpha = (uint64_t *) R_alloc(k + 1, sizeof(uint64_t));

    for (i = 0; i <= k; i++)
        powAlpha[i] = pow(alphaInf->numAlphabetChars, i);

    if (mapIndex)
    {
        if (useHash)
        {
            hmap = (khash_t(fim) *) indexMap;

            for (iter = kh_begin(hmap); iter != kh_end(hmap); iter++)
            {
                if (kh_exist(hmap, iter))
                {
                    for (j = 0; j < k; j++)
                    {
                        kmer[k - j - 1] =
                        alphaInf->reverseIndexMap[(int)((kh_key(hmap, iter) % powAlpha[j + 1]) / powAlpha[j])];
                    }

                    SET_STRING_ELT(colnames, kh_value(hmap, iter), Rf_mkChar(kmer));
                }
            }
        }
        else
        {
            featIndexMap = (uint32_t *) indexMap;
            l = 0;

            for (i = 0; i < (int)dimFeatureSpace; i++)
            {
                if (featIndexMap[i] != MAXUINT32)
                {
                    for (j = 0; j < k; j++)
                    {
                        kmer[k - j - 1] =
                            alphaInf->reverseIndexMap[(int)((i % (int) powAlpha[j + 1]) / (int) powAlpha[j])];
                    }

                    SET_STRING_ELT(colnames, l++, Rf_mkChar(kmer));
                }
            }
        }
    }
    else
    {
        for (i = 0; i < (int) dimFeatureSpace; i++)
        {
            if (i % USER_INTERRUPT_LIMIT == 0)
                R_CheckUserInterrupt();

            for (j = 0; j < k; j++)
                kmer[k - j - 1] =
                    alphaInf->reverseIndexMap[(int)((i % (int) powAlpha[j + 1]) / (int) powAlpha[j])];

            SET_STRING_ELT(colnames, i, Rf_mkChar(kmer));
        }
    }
}

RcppExport SEXP genExplRepMismatch(ByteStringVector x, int sizeX, IntegerVector selX, int maxSeqLength,
                                   int bioCharset, ByteStringVector features, int k, int m, bool presence,
                                   bool normalized, bool unmapped, bool lowercase, bool useRowNames,
                                   bool useColNames, bool zeroFeatures, bool sparse)
{
    int i, numProtect, numUsedFeatures, numNonzeroFeatures;
    uint64_t dimFeatureSpace;
    double *normValues;
    bool mapIndex, useHash;
    const void *vmax;
    void *indexMap;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    SEXP ers, dims, slot_p, slot_j, slot_x, rownames, colnames, dimnames;

    pfeatIndexMap = NULL;
    pFeatureHMap = NULL;
    pFeatureCounts = NULL;
    pFeatureCountsHMap = NULL;
    pIndexMap = NULL;
    pIndexHMap = NULL;
    hmap = NULL;
    indexMap = NULL;
    normValues = NULL;
    pNormValues = NULL;

    // check empty features vector
    if (features.length == 0)
        return(generateEmptyExplicitRep(sizeX, sparse));

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);
    dimFeatureSpace = pow(alphaInf.numAlphabetChars, k);

    // check if features space to large
    if (dimFeatureSpace > MAXINDEX32 && zeroFeatures && features.length == 0)
        return(generateEmptyExplicitRep(sizeX, sparse));

    if (dimFeatureSpace > FEATURE_SPACE_LIMIT)
    {
        Rprintf("feature space too large\n");
        return(generateEmptyExplicitRep(sizeX, sparse));
    }

    if (dimFeatureSpace <= HASH_MAP_LIMIT)
        useHash = FALSE;
    else
        useHash = TRUE;

    if (!zeroFeatures || features.length > 0)
        mapIndex = TRUE;
    else
    {
        mapIndex = FALSE;
        numUsedFeatures = dimFeatureSpace;
    }

    if (!sparse)
    {
        vmax = vmaxget();

        if (!zeroFeatures || features.length > 0)
        {
            if (!getIndexMap(x, sizeX, selX, unmapped, bioCharset, k, m, normalized, presence, &alphaInf,
                             features, dimFeatureSpace, zeroFeatures, useHash, &indexMap, &numUsedFeatures,
                             FALSE, NULL, &normValues))
            {
                vmaxset(vmax);
                return(generateEmptyExplicitRep(sizeX, sparse));
            }
        }

        // check if no features found
        if (numUsedFeatures < 1)
        {
            vmaxset(vmax);
            return(generateEmptyExplicitRep(sizeX, sparse));
        }

        // check if too many features
        if ((uint32_t) numUsedFeatures > MAX_FEATURES)
        {
            Rprintf("Too many features for explicit representation");
            vmaxset(vmax);
            return(generateEmptyExplicitRep(sizeX, sparse));
        }

        // free heap
        vmaxset(vmax);
        NumericMatrix erd(sizeX, numUsedFeatures);
        numProtect = 0;

        if (useColNames)
            PROTECT(colnames = Rf_allocVector(STRSXP, numUsedFeatures));
        else
            PROTECT(colnames = Rf_allocVector(STRSXP, 0));

        PROTECT(rownames = Rf_allocVector(STRSXP, 0));
        PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, rownames);
        SET_VECTOR_ELT(dimnames, 1, colnames);
        Rf_setAttrib(erd, R_DimNamesSymbol, dimnames);
        numProtect = 3;

        if (useColNames)
        {
            assignFeatureNames(colnames, indexMap, k, &alphaInf, dimFeatureSpace, mapIndex,
                               useHash);
        }

        vmax = vmaxget();

        getERDMismatch(x, sizeX, selX, unmapped, bioCharset, k, m, normalized, presence,
                       &alphaInf, features, dimFeatureSpace, zeroFeatures, useHash, indexMap,
                       numUsedFeatures, erd, normValues);

        vmaxset(vmax);

        if (numProtect > 0)
            UNPROTECT(numProtect);

        return(erd);
    }
    else
    {
        vmax = vmaxget();

        numNonzeroFeatures = 0;

        if (!getIndexMap(x, sizeX, selX, unmapped, bioCharset, k, m, normalized, presence,
                         &alphaInf, features, dimFeatureSpace, zeroFeatures, useHash, &indexMap,
                         &numUsedFeatures, TRUE, &numNonzeroFeatures, &normValues))
        {
            vmaxset(vmax);
            return(generateEmptyExplicitRep(sizeX, sparse));
        }

        // check if no features found
        if (numUsedFeatures < 1)
        {
            vmaxset(vmax);
            return(generateEmptyExplicitRep(sizeX, sparse));
        }

        // check if too many features
        if ((uint32_t) numUsedFeatures > MAX_FEATURES)
        {
            Rprintf("Too many features for explicit representation");
            vmaxset(vmax);
            return(generateEmptyExplicitRep(sizeX, sparse));
        }

        // free heap
        vmaxset(vmax);

        // allocate ExplicitRepSparse
        numProtect = 0;
        ers = PROTECT(NEW_OBJECT(MAKE_CLASS("ExplicitRepresentationSparse")));
        numProtect++;
        dims = PROTECT(Rf_allocVector(INTSXP, 2));
        numProtect++;
        SET_SLOT(ers, Rf_mkChar("Dim"), dims);
        INTEGER(dims)[0] = sizeX;
        INTEGER(dims)[1] = numUsedFeatures;
        slot_p = PROTECT(Rf_allocVector(INTSXP, sizeX + 1));
        numProtect++;
        SET_SLOT(ers, Rf_mkChar("p"), slot_p);

        if (useRowNames || useColNames)
        {
            PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
            PROTECT(rownames = Rf_allocVector(STRSXP, 0));
            numProtect += 2;

            if (useColNames && numUsedFeatures > 0)
            {
                PROTECT(colnames = Rf_allocVector(STRSXP, numUsedFeatures));
                numProtect++;
            }
            else
            {
                PROTECT(colnames = Rf_allocVector(STRSXP, 0));
                numProtect++;
            }

            SET_VECTOR_ELT(dimnames, 0, rownames);
            SET_VECTOR_ELT(dimnames, 1, colnames);
            SET_SLOT(ers, Rf_mkChar("Dimnames"), dimnames);
        }
        else
            colnames = NULL;

        if (numNonzeroFeatures == 0)
        {
            for (i = 0; i < (sizeX+1); i++)
                INTEGER(slot_p)[i] = 0;

            UNPROTECT(numProtect);
            return(ers);
        }
        else
        {
            slot_j = PROTECT(Rf_allocVector(INTSXP, numNonzeroFeatures));
            numProtect++;
            SET_SLOT(ers, Rf_mkChar("j"), slot_j);

            slot_x = PROTECT(Rf_allocVector(REALSXP, numNonzeroFeatures));
            numProtect++;
            SET_SLOT(ers, Rf_mkChar("x"), slot_x);
        }

        vmax = vmaxget();

        if (useColNames)
            assignFeatureNames(colnames, indexMap, k, &alphaInf, dimFeatureSpace, mapIndex, useHash);

        getERSMismatch(x, sizeX, selX, maxSeqLength, unmapped, bioCharset, k, m, normalized, presence,
                       &alphaInf, features, dimFeatureSpace, zeroFeatures, mapIndex, useHash, indexMap,
                       numUsedFeatures, slot_p, slot_j, slot_x, normValues);

        vmaxset(vmax);

        if (numProtect > 0)
            UNPROTECT(numProtect);

        return(ers);
    }
}

RcppExport SEXP mismatchKernelMatrixC(SEXP xR, SEXP yR, SEXP selXR, SEXP selYR, SEXP sizeXR, SEXP sizeYR,
                                      SEXP isXStringSetR, SEXP symmetricR, SEXP bioCharsetR, SEXP lowercaseR,
                                      SEXP unmappedR, SEXP maxSeqLengthR, SEXP kR, SEXP mR, SEXP normalizedR,
                                      SEXP presenceR)
{
    int sizeX = as<int>(sizeXR);
    int sizeY = as<int>(sizeYR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    bool symmetric = as<bool>(symmetricR);
    bool isXStringSet = as<bool>(isXStringSetR);
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    const void *vmax;

    if (symmetric)
        sizeY = sizeX;

    NumericMatrix km(sizeX, sizeY);

    vmax = vmaxget();

    ByteStringVector x, y;
    IntegerVector selX(selXR);
    IntegerVector selY(selYR);

    if (isXStringSet)
        x = XStringSet2ByteStringVec(xR);
    else
        x = charVector2ByteStringVec(xR);

    if (!Rf_isNull(yR))
    {
        if (isXStringSet)
            y = XStringSet2ByteStringVec(yR);
        else
            y = charVector2ByteStringVec(yR);
    }
    else
        y.length = 0;

    int k = as<int>(kR);
    int m = as<int>(mR);
    int bioCharset = as<int>(bioCharsetR);
    bool lowercase = as<bool>(lowercaseR);
    bool unmapped = as<bool>(unmappedR);
    bool normalized = as<bool>(normalizedR);
    bool presence = as<bool>(presenceR);

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    // no support for position or annotation specific kernel

    getMismatchKernelMatrix(km, x, y, sizeX, sizeY, selX, selY, symmetric, bioCharset, lowercase,
                            unmapped, k, m, normalized, presence, maxSeqLength, &alphaInf);
    vmaxset(vmax);

    return(km);
}

RcppExport void freeHeapPredProfMismatch()
{
    if (hmap != NULL)
    {
        kh_destroy(fw, hmap);
        hmap = NULL;
    }
}

RcppExport void featuresToHashmapMismatch(NumericMatrix featureWeights, int svmIndex, int k, int m,
                                          struct alphaInfo *alphaInf)
{
    int i, j, numFeatures, result;
    uint64_t featureIndex;
    const char *pattern;
    SEXP dimNames, colNames;
    khiter_t iter;
    struct hmData featureData;

    hmap = kh_init(fw);

    numFeatures = featureWeights.ncol();

    dimNames = Rf_getAttrib(featureWeights, R_DimNamesSymbol);
    colNames = VECTOR_ELT(dimNames, 1);

    featureData.unweightedPosIndex = MAXUINT32;

    for (i = 0; i < numFeatures; i++)
    {
        pattern = CHAR(STRING_ELT(colNames, i));
        featureIndex = 0;

        for (j = 0; j < k; j++)
            featureIndex = featureIndex * alphaInf->numAlphabetChars + alphaInf->indexMap[(int) pattern[j]];

        iter = kh_put(fw, hmap, featureIndex, &result);

        if (result == -1)
        {
            Rprintf("Storage of key %llu in hashmap failed\n", featureIndex);
            return;
        }

        featureData.featWeight = featureWeights(svmIndex, i);
        kh_value(hmap, iter) = featureData;
    }

    // dealloc of C heap needs to be done by on.exit hook

    return;
}

void mutateFeatureIndexPredProf(uint64_t featureIndex, int svmIndex, int sampleIndex, int currPos, int k, int m,
                                uint64_t *powAlpha, NumericMatrix pprof, NumericMatrix featureWeights, bool useHash,
                                bool calcKernelValue, khash_t(fc) *fchmap, double *kernelValue, bool presence)
{
    int i, l, currStack, numMismatches, pos, result;
    uint32_t stack[6*k], low, currCount;
    uint64_t currFeatureIndex, origIndex, lower;
    bool finished, mismatchPositions[k];
    struct hmData featureData;
    khiter_t iter;

    currStack = -1;
    finished = FALSE;
    currFeatureIndex = featureIndex;
    numMismatches = 0;
    origIndex = 0;
    pos = 0;
    i = 0;

    while (!finished)
    {
        if (numMismatches == m || pos == k)
        {
            // get mismatch positions
            lower = 1;
            memset(mismatchPositions, 0, k);

            for (i=0; i < k; i++)
            {
                if (((featureIndex / lower) % k) != ((currFeatureIndex / lower) % k))
                    mismatchPositions[i] = TRUE;

                lower *= k;
            }

            if (calcKernelValue)
            {
                iter = kh_get(fc, fchmap, featureIndex);

                if (iter != kh_end(fchmap))
                {
                    if (!presence)
                    {
                        currCount = kh_value(fchmap, iter);
                        kh_value(fchmap, iter) = currCount + 1;
                        *kernelValue = *kernelValue - currCount * currCount +
                                       (currCount + 1) * (currCount + 1);
                    }
                }
                else
                {
                    iter = kh_put(fc, fchmap, featureIndex, &result);

                    if (result == -1)
                    {
                        Rprintf("Storage of key %llu in feature count hashmap failed\n",
                                featureIndex);
                        return;
                    }

                    kh_value(fchmap, iter) = 1;
                    *kernelValue += 1;
                }
            }

            if (useHash)
            {
                // check if new feature
                iter = kh_get(fw, hmap, featureIndex);

                if (iter != kh_end(hmap))
                {
                    featureData = kh_value(hmap, iter);
                    featureData.featWeight /= (k - numMismatches);

                    for (l = 0; l < k; l++)
                    {
                        if (!mismatchPositions[l])
                            pprof(sampleIndex, currPos+1-k+l) += featureData.featWeight;
                    }
                }
            }
            else
            {
                for (l = 0; l < k; l++)
                {
                    if (!mismatchPositions[l])
                    {
                        pprof(sampleIndex, currPos+1-k+l) += featureWeights(svmIndex, featureIndex)
                                                         / (k - numMismatches);
                    }
                }
            }

            if (currStack > -1)
            {
                // retrieve previous state and continue
                origIndex = stack[currStack--];
                numMismatches = stack[currStack--];
                pos = stack[currStack--];
                i = stack[currStack--];
                low = stack[currStack--];
                currFeatureIndex = (((uint64_t) stack[currStack--]) << 32) + low;
            }
            else
                finished = TRUE;
        }

        if (i == 0)
        {
            origIndex = (currFeatureIndex / powAlpha[pos]) % powAlpha[1];
            currFeatureIndex = currFeatureIndex - (origIndex) * powAlpha[pos];
        }

        if (i < (int)powAlpha[1])
        {
            // push state to stack and descend
            stack[++currStack] = (uint32_t) ((currFeatureIndex + powAlpha[pos]) >> 32);
            stack[++currStack] = (uint32_t) ((currFeatureIndex + powAlpha[pos]) & LOWER_LONG);
            stack[++currStack] = i + 1;
            stack[++currStack] = pos;
            stack[++currStack] = numMismatches;
            stack[++currStack] = origIndex;

            if (i != (int) origIndex)
                numMismatches++;

            i = 0;
            pos++;
        }
        else
        {
            if (currStack > -1)
            {
                // retrieve previous state and continue
                origIndex = stack[currStack--];
                numMismatches = stack[currStack--];
                pos = stack[currStack--];
                i = stack[currStack--];
                low = stack[currStack--];
                currFeatureIndex = (((uint64_t) stack[currStack--]) << 32) + low;
            }
            else
                finished = TRUE;
        }
    }
}

void genPredProfileMismatch(NumericMatrix pprof, ByteStringVector x, IntegerVector selX, int numSamples,
                            ByteStringVector annCharset, ByteStringVector annX, int maxSeqLength, bool unmapped,
                            int kernelType, int k, int m, int bioCharset, NumericMatrix featureWeights,
                            int svmIndex, bool lowercase, bool normalized, bool presence)
{
    int i, j, iX, iold, patLength, index;
    bool useHash = FALSE, calcKernelValue=FALSE;
    uint64_t dimFeatureSpace, featureIndex, prevIndex;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    double kernelValue, *normValues;
    khash_t(fc) *fchmap;

    pfeatIndexMap = NULL;
    pFeatureHMap = NULL;
    pFeatureCounts = NULL;
    pFeatureCountsHMap = NULL;
    pIndexMap = NULL;
    pIndexHMap = NULL;
    normValues = NULL;
    pNormValues = NULL;
    hmap = NULL;
    fchmap = NULL;

    dimFeatureSpace = pow((uint64_t) alphaInf.numAlphabetChars, k);

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    // if not all features present use hash map
    if ((uint64_t) featureWeights.ncol() != dimFeatureSpace)
    {
        useHash = TRUE;
        featuresToHashmapMismatch(featureWeights, svmIndex, k, m, &alphaInf);
    }

    if (normalized)
    {
        normValues = (double *) Calloc(numSamples, double);
        pNormValues = normValues;
        fchmap = kh_init(fc);
        pFeatureCountsHMap = fchmap;
        calcKernelValue = TRUE;
    }

    // cyclic buffers for old index
    uint64_t *oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    uint64_t *powAlpha = (uint64_t *) R_alloc(k + 1, sizeof(uint64_t));

    for (i = 0; i <= k; i++)
        powAlpha[i] = pow(alphaInf.numAlphabetChars, i);

    for (i = 0 ; i < numSamples; i++)
    {
        featureIndex = 0;
        patLength = 0;
        iold = 0;
        iX = selX[i];
        kernelValue = 0;

        for (j = 0; j < x.nchar[iX]; j++) //
        {
            index = alphaInf.seqIndexMap[(int) x.ptr[iX][j]];

            if (index > -1)
            {
                if (patLength < k)
                {
                    oldIndex[iold++] = index * powAlpha[k-1];

                    if (iold == k)
                        iold = 0;

                    featureIndex = featureIndex * alphaInf.numAlphabetChars + index;

                    patLength++;

                    if (patLength == k)
                    {
                        mutateFeatureIndexPredProf(featureIndex, svmIndex, i, j, k, m, powAlpha,
                                                   pprof, featureWeights, useHash, calcKernelValue,
                                                   fchmap, &kernelValue, presence);
                    }
                }
                else
                {
                    prevIndex = oldIndex[iold];
                    oldIndex[iold++] = index * powAlpha[k-1];

                    if (iold == k)
                        iold = 0;

                    featureIndex = (featureIndex - prevIndex) *
                    alphaInf.numAlphabetChars + index;

                    mutateFeatureIndexPredProf(featureIndex, svmIndex, i, j, k, m, powAlpha,
                                               pprof, featureWeights, useHash, calcKernelValue,
                                               fchmap, &kernelValue, presence);
                }
            }
            else // wrong character - restart kmer
            {
                featureIndex = 0;
                patLength = 0;
                iold = 0;
            }
        }

        if (calcKernelValue)
        {
            kh_clear(fc, fchmap);
            normValues[i] = sqrt(kernelValue);
        }
    }

    if (normalized)
    {
        for (i = 0; i < numSamples; i++)
        {
            for (j = 0; j < pprof.ncol(); j++)
            {
                if (normValues[i] == 0)
                    continue;

                if (pprof(i,j) != 0)
                    pprof(i,j) /= normValues[i];
            }
        }
    }

    // dealloc of C heap done by on.exit hook

    return;
}

void freeHeapMismatch()
{
    if (pNormValues != NULL)
    {
        Free(pNormValues);
        pNormValues = NULL;
    }

    if (pfeatIndexMap != NULL)
    {
        Free(pfeatIndexMap);
        pfeatIndexMap = NULL;
    }

    if (pIndexMap != NULL)
    {
        Free(pIndexMap);
        pIndexMap = NULL;
    }

    if (pFeatureHMap != NULL)
    {
        kh_destroy(fim, pFeatureHMap);
        pFeatureHMap = NULL;
    }

    if (pIndexHMap != NULL)
    {
        kh_destroy(ind, pIndexHMap);
        pIndexHMap = NULL;
    }

    if (hmap != NULL)
    {
        kh_destroy(fw, hmap);
        hmap = NULL;
    }
}
