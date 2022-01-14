//
// C Routines for Gappy Pair Kernel
//
// Source : GappyPair.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014-2016   J o h a n n e s  P a l m e
//

#include "Kebabs.h"
#include "KernelUtils.h"

#define TF_NONE               0
#define TF_LEAF               1
#define TF_REV_COMPLEMENT     2

#define USER_INTERRUPT_LIMIT  100000
#define INDEX_MAP_INIT_FF     0xFF

extern "C"
{
    #include "BitArray.h"
    #include "IntegerPowers64.h"
    #include "khash.h"
    #include "ksort.h"
}

using namespace Rcpp;

struct hmData
{
    uint32_t featureIndex;
    uint32_t unweightedPosIndex;
};

#if __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

// hash map for mapping 64 bit feature indices to feature weights
KHASH_MAP_INIT_INT64(fw, struct hmData)
// hash map for mapping 64 bit feature indices to col indices
KHASH_MAP_INIT_INT64(fim, uint32_t)
// hash map for feature counting in case of normalization
KHASH_MAP_INIT_INT64(fc, uint32_t)
// sorting of 64 bit feature indices
KSORT_INIT(gappy, uint64_t, ks_lt_generic)


static khash_t(fw) *hmap;
static khash_t(fim) *pFeatureHMap;
static uint32_t *pFeatureMap;
static khash_t(fc) *pFeatureCountsHMap;
static uint32_t *pFeatureCounts;
static double *pNormValues;

#if __clang__
#pragma clang diagnostic pop
#endif

int findReverseComplementLeaf(const char *s, int slen, int currPos, const char* annptr, int k, int m, int g,
                              IntegerVector *annotationIndexMap, bool presence, struct prefTree *pTree,
                              int maxNoOfNodes, int *freeNode, struct indexBlock *nullBlock, bool *printWarning,
                              struct alphaInfo *alphaInf)
{

    int i, index, curr;

    curr = 0;

    if (annptr == NULL)
    {
        // create reverse path if necessary
        for (i = currPos + 2 * k + g - 1; i >= currPos + k + g; i--)
        {
            index = alphaInf->numAlphabetChars - alphaInf->seqIndexMap[(int)s[i]] - 1;

            // in this range sequence only contains valid characters
            if (pTree->node[curr].ib.idx[index] > 0)
                curr = pTree->node[curr].ib.idx[index];
            else
            {
                pTree->node[curr].ib.idx[index] = *freeNode;
                curr = *freeNode;

                if (curr >= maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }

                    return(0);
                }

                *freeNode = *freeNode + 1;
                pTree->node[curr].ib = *nullBlock;
                pTree->node[curr].leaf = TF_NONE;
            }
        }

        for (i = currPos + k - 1; i >= currPos; i--)
        {
            index = alphaInf->numAlphabetChars - alphaInf->seqIndexMap[(int)s[i]] - 1;

            // in this range sequence only contains valid characters
            if (pTree->node[curr].ib.idx[index] > 0)
                curr = pTree->node[curr].ib.idx[index];
            else
            {
                pTree->node[curr].ib.idx[index] = *freeNode;
                curr = *freeNode;

                if (curr >= maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }

                    return(0);
                }

                *freeNode = *freeNode + 1;

                if (i != currPos)
                {
                    pTree->node[curr].ib = *nullBlock;
                    pTree->node[curr].leaf = TF_NONE;
                }
                else
                {
                    pTree->node[curr].ib = *nullBlock;
                    pTree->node[curr].leaf = TF_LEAF;
                    pTree->node[curr].value = 0;
                }
            }
        }
    }
    else
    {
        for (i = currPos + 2 * k + g - 1; i >= currPos + k + g; i--)
        {
            index = alphaInf->numAlphabetChars - alphaInf->seqIndexMap[(int)s[i]] - 1;

            // in this range sequence only contains valid characters
            if (pTree->node[curr].ib.idx[index] > 0)
                curr = pTree->node[curr].ib.idx[index];
            else
            {
                pTree->node[curr].ib.idx[index] = *freeNode;
                curr = *freeNode;

                if (curr >= maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }

                    return(0);
                }

                *freeNode = *freeNode + 1;
                pTree->node[curr].ib = *nullBlock;
                pTree->node[curr].leaf = TF_NONE;
            }
        }

        for (i = currPos + k - 1; i >= currPos; i--)
        {
            index = alphaInf->numAlphabetChars - alphaInf->seqIndexMap[(int)s[i]] - 1;

            // in this range sequence only contains valid characters
            if (pTree->node[curr].ib.idx[index] > 0)
                curr = pTree->node[curr].ib.idx[index];
            else
            {
                pTree->node[curr].ib.idx[index] = *freeNode;
                curr = *freeNode;

                if (curr >= maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }

                    return(0);
                }

                *freeNode = *freeNode + 1;
                pTree->node[curr].ib = *nullBlock;
                pTree->node[curr].leaf = TF_NONE;
            }
        }

        for (i = currPos + 2 * k + g - 1; i >= currPos + k + g; i--)
        {
            index = (*annotationIndexMap)[annptr[i]];

            if (pTree->node[curr].ib.idx[index] > 0)
                curr = pTree->node[curr].ib.idx[index];
            else
            {
                pTree->node[curr].ib.idx[index] = *freeNode;
                curr = *freeNode;

                if (curr >= maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }

                    return(0);
                }

                *freeNode = *freeNode + 1;
                pTree->node[curr].ib = *nullBlock;
                pTree->node[curr].leaf = TF_NONE;
            }
        }

        for (i = currPos + k - 1; i >= currPos; i--)
        {
            index = (*annotationIndexMap)[annptr[i]];

            if (pTree->node[curr].ib.idx[index] > 0)
                curr = pTree->node[curr].ib.idx[index];
            else
            {
                pTree->node[curr].ib.idx[index] = *freeNode;
                curr = *freeNode;

                if (curr >= maxNoOfNodes)
                {
                    if (*printWarning)
                    {
                        Rprintf("Maximum number of nodes exceeded\n");
                        *printWarning = FALSE;
                    }

                    return(0);
                }

                *freeNode = *freeNode + 1;

                if (i != currPos)
                {
                    pTree->node[curr].ib = *nullBlock;
                    pTree->node[curr].leaf = TF_NONE;
                }
                else
                {
                    pTree->node[curr].ib = *nullBlock;
                    pTree->node[curr].leaf = TF_LEAF;
                    pTree->node[curr].value = 0;
                }
            }
        }
    }

    if (!(pTree->node[curr].leaf & TF_LEAF))
        return(0);

    return(curr);
}

double createTreeGappy(const char *s, int slen, const char* annptr, int k, int m,
                       IntegerVector *annotationIndexMap, bool presence,
                       bool reverseComplement, struct prefTree *pTree, int maxNoOfNodes,
                       int *freeNode, struct indexBlock *nullBlock, bool *printWarning,
                       struct alphaInfo *alphaInf)
{
    int i, j, g, l, offset, curr, curr2, currReverse, index;
    uint64_t forwardIndex, reverseIndex;
    double sum = 0;

    // init first node
    pTree->node[0].ib = *nullBlock;
    pTree->node[0].leaf = FALSE;

    if (annptr == NULL)
    {
        // build prefix tree for sequence
        for (i = 0; i <= slen - 2*k; i++)
        {
            curr = 0;

            for (j = 0; j < k; j++)
            {
                index = alphaInf->seqIndexMap[(int)s[i+j]];

                if (index > -1)
                {
                    if (pTree->node[curr].ib.idx[index] > 0)
                        curr = pTree->node[curr].ib.idx[index];
                    else
                    {
                        pTree->node[curr].ib.idx[index] = *freeNode;
                        curr = *freeNode;

                        if (curr >= maxNoOfNodes)
                        {
                            if (*printWarning)
                            {
                                Rprintf("Maximum number of nodes exceeded\n");
                                *printWarning = FALSE;
                            }
                            return(NA_REAL);
                        }

                        *freeNode = *freeNode + 1;

                        pTree->node[curr].ib = *nullBlock;
                        pTree->node[curr].leaf = FALSE;
                    }
                }
                else
                {
                    // interrupt kmer when unknown character found
                    curr = 0;
                    break;
                }
            }

            if (curr > 0)
            {
                curr2 = curr;

                for (g = 0; g <= m; g++)
                {
                    if (curr == 0 && index < 0)
                        break;

                    curr = curr2;
                    offset = i + k + g;

                    if (offset + k > slen)
                    {
                        curr = 0;
                        break;
                    }

                    for (j = 0; j < k; j++)
                    {
                        index = alphaInf->seqIndexMap[(int) s[offset+j]];

                        if (index > -1)
                        {
                            if (pTree->node[curr].ib.idx[index] > 0)
                            {
                                curr = pTree->node[curr].ib.idx[index];

                                if (j == k - 1)
                                {
                                    // $$$TODO Remove
                                    if (pTree->node[curr].leaf & TF_LEAF)
                                    {
                                        if (!((pTree->node[curr].leaf & TF_REV_COMPLEMENT)))
                                        {
                                            if (!presence)
                                            {
                                                pTree->node[curr].ib.idx[g] = pTree->node[curr].ib.idx[g] + 1;
                                                sum = sum - (double)((pTree->node[curr].ib.idx[g] - 1) *
                                                                     (pTree->node[curr].ib.idx[g] - 1))
                                                + (double)(pTree->node[curr].ib.idx[g] *
                                                           pTree->node[curr].ib.idx[g]);
                                            }
                                            else
                                            {
                                                if (pTree->node[curr].ib.idx[g] == 0)
                                                {
                                                    pTree->node[curr].ib.idx[g] = 1;
                                                    sum++;
                                                }
                                            }
                                            curr = 0;
                                        }
                                        else
                                        {
                                            currReverse = pTree->node[curr].value;

                                            if (currReverse == 0)
                                            {
                                                if (*printWarning)
                                                {
                                                    Rprintf("Leaf for reverse complement not found\n");
                                                    *printWarning = FALSE;
                                                }
                                                return(NA_REAL);
                                            }

                                            if (!presence)
                                            {
                                                pTree->node[currReverse].ib.idx[g] = pTree->node[currReverse].ib.idx[g] + 1;
                                                sum = sum - (double)((pTree->node[currReverse].ib.idx[g] - 1) *
                                                                     (pTree->node[currReverse].ib.idx[g] - 1))
                                                          + (double)(pTree->node[currReverse].ib.idx[g] *
                                                                     pTree->node[currReverse].ib.idx[g]);
                                            }
                                            else
                                            {
                                                if (pTree->node[currReverse].ib.idx[g] == 0)
                                                {
                                                    pTree->node[currReverse].ib.idx[g] = 1;
                                                    sum++;
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (*printWarning)
                                        {
                                            Rprintf("Invalid leaf reached:\n");
                                            Rprintf("    curr: %d, i: %d, j: %d, g: %d\n", curr, i, j, g);
                                        }
                                        return(NA_REAL);
                                    }

                                }
                            }
                            else
                            {
                                pTree->node[curr].ib.idx[index] = *freeNode;
                                curr = *freeNode;

                                if (curr >= maxNoOfNodes)
                                {
                                    if (*printWarning)
                                    {
                                        Rprintf("Maximum number of nodes exceeded\n");
                                        *printWarning = FALSE;
                                    }
                                    return(NA_REAL);
                                }

                                *freeNode = *freeNode + 1;

                                if (j == k - 1)
                                {
                                    pTree->node[curr].leaf = TF_LEAF;

                                    if (!reverseComplement)
                                    {
                                        pTree->node[curr].ib = *nullBlock;
                                        pTree->node[curr].ib.idx[g] = 1;
                                        sum += 1;
                                    }
                                    else
                                    {
                                        forwardIndex = 0;
                                        reverseIndex = 0;

                                        for (l = i; l < i + k; l++)
                                        {
                                            forwardIndex = forwardIndex * alphaInf->numAlphabetChars +
                                                           alphaInf->seqIndexMap[(int)s[l]];
                                        }

                                        for (l = i + k + g; l < i + 2*k + g; l++)
                                        {
                                            forwardIndex = forwardIndex * alphaInf->numAlphabetChars +
                                                           alphaInf->seqIndexMap[(int)s[l]];
                                        }

                                        for (l = i + 2*k + g - 1; l >= i + k + g; l--)
                                        {
                                            reverseIndex = (reverseIndex + 1) * alphaInf->numAlphabetChars -
                                            alphaInf->seqIndexMap[(int)s[l]] - 1;
                                        }

                                        for (l = i + k - 1; l >= i; l--)
                                        {
                                            reverseIndex = (reverseIndex + 1) * alphaInf->numAlphabetChars -
                                            alphaInf->seqIndexMap[(int)s[l]] - 1;
                                        }

                                        if (forwardIndex <= reverseIndex)
                                        {
                                            pTree->node[curr].ib = *nullBlock;
                                            pTree->node[curr].value = 1;
                                            sum += 1;
                                        }
                                        else
                                        {
                                            pTree->node[curr].leaf |= TF_REV_COMPLEMENT;

                                            currReverse = findReverseComplementLeaf(s, slen, i, annptr, k, m, g,
                                                                                    annotationIndexMap, presence,
                                                                                    pTree, maxNoOfNodes, freeNode,
                                                                                    nullBlock, printWarning, alphaInf);

                                            if (currReverse == 0)
                                            {
                                                if (*printWarning)
                                                {
                                                    Rprintf("Leaf for reverse complement not found\n");
                                                    *printWarning = FALSE;
                                                }
                                                return(NA_REAL);
                                            }

                                            pTree->node[curr].value = currReverse;

                                            if (!presence)
                                            {
                                                pTree->node[currReverse].ib.idx[g] = pTree->node[currReverse].ib.idx[g] + 1;
                                                sum = sum - (double)((pTree->node[currReverse].ib.idx[g] - 1) *
                                                                     (pTree->node[currReverse].ib.idx[g] - 1))
                                                          + (double)(pTree->node[currReverse].ib.idx[g] *
                                                                     pTree->node[currReverse].ib.idx[g]);
                                            }
                                            else
                                            {
                                                if (pTree->node[currReverse].ib.idx[g] == 0)
                                                {
                                                    pTree->node[currReverse].ib.idx[g] = 1;
                                                    sum++;
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    pTree->node[curr].ib = *nullBlock;
                                    pTree->node[curr].leaf = TF_NONE;
                                }
                            }
                        }
                        else
                        {
                            // interrupt kmer when unknown character found
                            curr = 0;
                            break;
                        }
                    }
                }
            }
        }
    }
    else
    {
        // build prefix tree for sequence
        for (int i=0; i<=slen-2*k; i++)
        {
            curr = 0;

            for (j = 0; j < k; j++)
            {
                index = alphaInf->seqIndexMap[(int) s[i+j]];

                if (index > -1)
                {
                    if (pTree->node[curr].ib.idx[index] > 0)
                        curr = pTree->node[curr].ib.idx[index];
                    else
                    {
                        pTree->node[curr].ib.idx[index] = *freeNode;
                        curr = *freeNode;

                        if (curr >= maxNoOfNodes)
                        {
                            if (*printWarning)
                            {
                                Rprintf("Maximum number of nodes exceeded\n");
                                *printWarning = FALSE;
                            }
                            return(NA_REAL);
                        }

                        *freeNode = *freeNode + 1;

                        pTree->node[curr].ib = *nullBlock;
                        pTree->node[curr].leaf = FALSE;
                    }
                }
                else
                {
                    // interrupt kmer when unknown character found
                    curr = 0;
                    break;
                }
            }

            if (curr > 0)
            {
                curr2 = curr;

                for (g = 0; g <= m; g++)
                {
                    if (curr == 0 && index < 0)
                        break;

                    curr = curr2;
                    offset = i + k + g;

                    if (offset + k > slen)
                    {
                        curr = 0;
                        break;
                    }

                    for (j = 0; j < k; j++)
                    {
                        index = alphaInf->seqIndexMap[(int) s[offset + j]];

                        if (index > -1)
                        {
                            if (pTree->node[curr].ib.idx[index] > 0)
                                curr = pTree->node[curr].ib.idx[index];
                            else
                            {
                                pTree->node[curr].ib.idx[index] = *freeNode;
                                curr = *freeNode;
                                if (curr >= maxNoOfNodes)
                                {
                                    if (*printWarning)
                                    {
                                        Rprintf("Maximum number of nodes exceeded\n");
                                        *printWarning = FALSE;
                                    }
                                    return(NA_REAL);
                                }

                                *freeNode = *freeNode + 1;

                                pTree->node[curr].ib = *nullBlock;
                                pTree->node[curr].leaf = FALSE;
                            }
                        }
                        else
                        {
                            // interrupt kmer when unknown character found
                            curr = 0;
                            break;
                        }
                    }

                    if (curr > 0)
                    {
                        for (j = 0; j < k; j++)
                        {
                            index = (*annotationIndexMap)[annptr[i+j]];

                            if (index > -1)
                            {
                                if (pTree->node[curr].ib.idx[index] > 0)
                                {
                                    curr = pTree->node[curr].ib.idx[index];
                                }
                                else
                                {
                                    pTree->node[curr].ib.idx[index] = *freeNode;
                                    curr = *freeNode;

                                    if (curr >= maxNoOfNodes)
                                    {
                                        if (*printWarning)
                                        {
                                            Rprintf("Maximum number of nodes exceeded\n");
                                            *printWarning = FALSE;
                                        }
                                        return(NA_REAL);
                                    }

                                    *freeNode = *freeNode + 1;

                                    pTree->node[curr].ib = *nullBlock;
                                    pTree->node[curr].leaf = FALSE;
                                }
                            }
                            else
                            {
                                // interrupt kmer when unknown character found
                                curr = 0;
                                break;
                            }
                        }

                        if (curr > 0)
                        {
                            for (j = 0; j < k; j++)
                            {
                                index = (*annotationIndexMap)[annptr[offset+j]];

                                if (index > -1)
                                {
                                    if (pTree->node[curr].ib.idx[index] > 0)
                                    {
                                        curr = pTree->node[curr].ib.idx[index];

                                        if (j == k - 1)
                                        {
                                            // $$$TODO Remove
                                            if ((pTree->node[curr].leaf & TF_LEAF))
                                            {
                                                if (!((pTree->node[curr].leaf & TF_REV_COMPLEMENT)))
                                                {
                                                    if (!presence)
                                                    {
                                                        pTree->node[curr].ib.idx[g] = pTree->node[curr].ib.idx[g] + 1;
                                                        sum = sum - (double)((pTree->node[curr].ib.idx[g] - 1) *
                                                                             (pTree->node[curr].ib.idx[g] - 1))
                                                                  + (double)(pTree->node[curr].ib.idx[g] *
                                                                             pTree->node[curr].ib.idx[g]);
                                                    }
                                                    else
                                                    {
                                                        if (pTree->node[curr].ib.idx[g] == 0)
                                                        {
                                                            pTree->node[curr].ib.idx[g] = 1;
                                                            sum++;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    currReverse = pTree->node[curr].value;

                                                    if (currReverse == 0)
                                                    {
                                                        if (*printWarning)
                                                        {
                                                            Rprintf("Leaf for reverse complement not found\n");
                                                            *printWarning = FALSE;
                                                        }
                                                        return(NA_REAL);
                                                    }

                                                    if (!presence)
                                                    {
                                                        pTree->node[currReverse].ib.idx[g] = pTree->node[currReverse].ib.idx[g] + 1;
                                                        sum = sum - (double)((pTree->node[currReverse].ib.idx[g] - 1) *
                                                                             (pTree->node[currReverse].ib.idx[g] - 1))
                                                                  + (double)(pTree->node[currReverse].ib.idx[g] *
                                                                             pTree->node[currReverse].ib.idx[g]);
                                                    }
                                                    else
                                                    {
                                                        if (pTree->node[currReverse].ib.idx[g] == 0)
                                                        {
                                                            pTree->node[currReverse].ib.idx[g] = 1;
                                                            sum++;
                                                        }
                                                    }
                                                }

                                                curr = 0;
                                            }
                                            else
                                            {
                                                if (*printWarning)
                                                {
                                                    Rprintf("Invalid leaf reached:\n");
                                                    Rprintf("    curr: %d, i: %d, j: %d, g: %d\n", curr, i, j, g);
                                                }
                                                return(NA_REAL);
                                            }
                                        }
                                    }
                                    else
                                    {
                                        pTree->node[curr].ib.idx[index] = *freeNode;
                                        curr = *freeNode;

                                        if (curr >= maxNoOfNodes)
                                        {
                                            if (*printWarning)
                                            {
                                                Rprintf("Maximum number of nodes exceeded\n");
                                                *printWarning = FALSE;
                                            }
                                            return(NA_REAL);
                                        }

                                        *freeNode = *freeNode + 1;

                                        if (j == k - 1)
                                        {
                                            pTree->node[curr].leaf = TF_LEAF;

                                            if (!reverseComplement)
                                            {
                                                pTree->node[curr].ib = *nullBlock;
                                                pTree->node[curr].ib.idx[g] = 1;
                                                sum += 1;
                                            }
                                            else
                                            {
                                                forwardIndex = 0;
                                                reverseIndex = 0;

                                                for (l = i; l < i + k; l++)
                                                {
                                                    forwardIndex = forwardIndex * alphaInf->numAlphabetChars +
                                                    alphaInf->seqIndexMap[(int)s[l]];
                                                }

                                                for (l = i + k + g; l < i + 2*k + g; l++)
                                                {
                                                    forwardIndex = forwardIndex * alphaInf->numAlphabetChars +
                                                    alphaInf->seqIndexMap[(int)s[l]];
                                                }

                                                for (l = i + 2*k + g - 1; l >= i + k + g; l--)
                                                {
                                                    reverseIndex = (reverseIndex + 1) * alphaInf->numAlphabetChars -
                                                    alphaInf->seqIndexMap[(int)s[l]] - 1;
                                                }

                                                for (l = i + k - 1; l >= i; l--)
                                                {
                                                    reverseIndex = (reverseIndex + 1) * alphaInf->numAlphabetChars -
                                                    alphaInf->seqIndexMap[(int)s[l]] - 1;
                                                }

                                                if (forwardIndex <= reverseIndex)
                                                {
                                                    pTree->node[curr].ib = *nullBlock;
                                                    pTree->node[curr].value = 1;
                                                    sum += 1;
                                                }
                                                else
                                                {
                                                    pTree->node[curr].leaf |= TF_REV_COMPLEMENT;

                                                    currReverse = findReverseComplementLeaf(s, slen, i, annptr, k, m, g,
                                                                                            annotationIndexMap, presence,
                                                                                            pTree, maxNoOfNodes, freeNode,
                                                                                            nullBlock, printWarning, alphaInf);

                                                    if (currReverse == 0)
                                                    {
                                                        if (*printWarning)
                                                        {
                                                            Rprintf("Leaf for reverse complement not found\n");
                                                            *printWarning = FALSE;
                                                        }
                                                        return(NA_REAL);
                                                    }

                                                    pTree->node[curr].value = currReverse;

                                                    if (!presence)
                                                    {
                                                        pTree->node[currReverse].ib.idx[g] = pTree->node[currReverse].ib.idx[g] + 1;
                                                        sum = sum - (double)((pTree->node[currReverse].ib.idx[g] - 1) *
                                                                             (pTree->node[currReverse].ib.idx[g] - 1))
                                                        + (double)(pTree->node[currReverse].ib.idx[g] *
                                                                   pTree->node[currReverse].ib.idx[g]);
                                                    }
                                                    else
                                                    {
                                                        if (pTree->node[currReverse].ib.idx[g] == 0)
                                                        {
                                                            pTree->node[currReverse].ib.idx[g] = 1;
                                                            sum++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else
                                        {
                                            pTree->node[curr].ib = *nullBlock;
                                            pTree->node[curr].leaf = TF_NONE;
                                        }
                                    }
                                }
                                else
                                {
                                    // interrupt kmer when unknown character found
                                    curr = 0;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return(sum);
}

template<typename T>
void getKMStdAnnGappy(T maxUnSignedIndex, NumericMatrix km, ByteStringVector x, ByteStringVector y,
                      int sizeX, int sizeY, IntegerVector selX, IntegerVector selY,
                      ByteStringVector annCharset, ByteStringVector annX, ByteStringVector annY,
                      int k, int m, bool normalized, bool symmetric, bool presence,
                      bool reverseComplement, int maxSeqLength, uint64_t dimFeatureSpace,
                      struct alphaInfo *alphaInf)
{
    T featureIndex, annotIndex;
    int i, j, l, iX, iY, elemIndex, noSamples = sizeX;
    int freeNode, nodeLimit, currBlock, currIndex, newBlock, maxBlockIndex;
    int maxNoOfNodes, seqnchar, currStack, stack[8*k];
    uint64_t maxNodesPerSequence, maxNumFeatures, numAnnPow2K = 0;
    double kv;
    bool printWarning = TRUE;
    const char *seqptr, *annptr;
    struct prefTree *pTree;
    struct indexBlock nullBlock;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);

    struct featVectorArrays {
        int32_t    *x;        // value array
        T          *idx;      // index array
    };

    // setup init mask
    for (i = 0; i < MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    if (!symmetric)
        noSamples += sizeY;

    // add one for the sentinel
    if (dimFeatureSpace < (uint64_t) maxSeqLength)
        maxNumFeatures = dimFeatureSpace + 1;
    else
        maxNumFeatures = (maxSeqLength - 2*k - 0.5*m + 1) * (m + 1) + 1;
    
    // allocate arrays for sparse feature vectors with 32 or 64 bit index
    // store only unnormalized k-mer count to avoid double space usage
    // add one for the sentinel
    featVectorArrays featVectors;
    featVectors.x   = (int32_t *) R_alloc(noSamples * maxNumFeatures, sizeof(int32_t));
    featVectors.idx = (T *) R_alloc(noSamples * maxNumFeatures, sizeof(T));

    double *normValues   = (double *) R_alloc(noSamples, sizeof(double));

    // alloc mem for prefix tree
    maxNoOfNodes = MAX_BLOCK;

    if (annX.length>0)
    {
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);

        nodeLimit = ((pow(alphaInf->numAlphabetChars, 2*k + 1) - 1) / (alphaInf->numAlphabetChars - 1)) +
                    pow(alphaInf->numAlphabetChars, 2*k) *
                          ((pow(annCharset.nchar[0], 2*k + 1) - 1) / (annCharset.nchar[0] - 1));
        maxNodesPerSequence = (3 * m + 4) * k * (maxSeqLength - 2 * k + 1) + 1;

        if (maxNodesPerSequence < (uint64_t) nodeLimit)
            nodeLimit = (int) maxNodesPerSequence;

        numAnnPow2K = ipow64(annCharset.nchar[0], 2*k);

        if (alphaInf->maxAlphaIndex > annCharset.nchar[0] - 1)
            maxBlockIndex = alphaInf->maxAlphaIndex;
        else
            maxBlockIndex = annCharset.nchar[0] - 1;
    }
    else
    {
        nodeLimit = (pow(alphaInf->numAlphabetChars, 2*k + 1) - 1) / (alphaInf->numAlphabetChars - 1);
        maxNodesPerSequence = (m + 2) * k * (maxSeqLength - 2 * k + 1) + 1;

        if (maxNodesPerSequence < (uint64_t) nodeLimit)
            nodeLimit = (int) maxNodesPerSequence;

        maxBlockIndex = alphaInf->maxAlphaIndex;
    }

    if (nodeLimit < maxNoOfNodes)
        maxNoOfNodes = nodeLimit;

    pTree = (struct prefTree *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    if (symmetric)
    {
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            iX = selX[i];
            annptr = NULL;

            if (annX.length > 0)
                annptr = annX.ptr[iX];

            freeNode = 1;

            //create new tree for sequence
            kv = createTreeGappy((const char *) x.ptr[iX], x.nchar[iX], annptr, k, m, &annotationIndexMap,
                                 presence, reverseComplement, pTree, maxNoOfNodes, &freeNode, &nullBlock,
                                 &printWarning, alphaInf);

            if (kv == NA_REAL)
            {
                for (j = 0; j < sizeX; j++)
                    for (l = 0; l < sizeX; l++)
                        km(j,l) = NA_REAL;

                return;
            }

            if (normalized)
                normValues[i] = sqrt(kv);
            else
                normValues[i] = kv;

            R_CheckUserInterrupt();

            // walk tree to generate sparse feature vector
            currBlock = 0;
            currIndex = 0;
            currStack = -1;
            elemIndex = i * maxNumFeatures;
            featureIndex = 0;
            annotIndex = 0;

            // for situation with no features set sentinel
            featVectors.idx[elemIndex] = maxUnSignedIndex;
            featVectors.x[elemIndex] = MAXINT32;

            while (currStack >= 0 || currIndex <= maxBlockIndex)
            {
                if (pTree->node[currBlock].ib.idx[currIndex] > 0)
                {
                    newBlock = pTree->node[currBlock].ib.idx[currIndex];

                    if (pTree->node[newBlock].leaf & TF_LEAF)
                    {
                        if (!(pTree->node[newBlock].leaf & TF_REV_COMPLEMENT))
                        {
                            for (j = 0; j <= m; j++)
                            {
                                if (pTree->node[newBlock].ib.idx[j] > 0)
                                {
                                    featVectors.x[elemIndex] = pTree->node[newBlock].ib.idx[j];

                                    if (currStack > 4 * k)
                                    {
                                        featVectors.idx[elemIndex++] =
                                        (featureIndex * numAnnPow2K + annotIndex * annCharset.nchar[0] + currIndex) * (m + 1) + j;
                                    }
                                    else
                                    {
                                        featVectors.idx[elemIndex++] =
                                        (featureIndex * alphaInf->numAlphabetChars + currIndex) * (m + 1) + j;
                                    }
                                }
                            }
                        }

                        currIndex++;

                        if (currIndex > maxBlockIndex)
                        {
                            if (currStack == -1)
                                continue;

                            while (currStack >= 0 && currIndex > maxBlockIndex)
                            {
                                currIndex = stack[currStack--];
                                currBlock = stack[currStack--];

                                if (currStack >= 4 * k - 2)
                                    annotIndex = (annotIndex - (currIndex - 1)) / annCharset.nchar[0];
                                else
                                    featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                            }
                        }
                    }
                    else
                    {
                        stack[++currStack] = currBlock;
                        stack[++currStack] = currIndex + 1;

                        // $$$ TODO remove check
                        if (currStack >= 8 * k)
                        {
                            Rprintf("Overflow of tree traversal stack\n");
                            return;
                        }

                        if (currStack > 4 * k)
                            annotIndex = annotIndex * annCharset.nchar[0] + currIndex;
                        else
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

                            if (currStack >= 4 * k - 2)
                                annotIndex = (annotIndex - (currIndex - 1)) / annCharset.nchar[0];
                            else
                                featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                        }
                    }
                }
            }

            featVectors.idx[elemIndex] = maxUnSignedIndex;
            featVectors.x[elemIndex] = MAXINT32;
        }

        computeKernelMatrix(maxUnSignedIndex, featVectors.idx, featVectors.x, km, normValues,
                            maxNumFeatures, sizeX, sizeY, normalized, symmetric);
    }
    else // not symmetric
    {
        int indexY = -1;

        // calculate kernel matrix
        for (i = 0; i < noSamples; i++)
        {
            R_CheckUserInterrupt();

            annptr = NULL;

            if (i < sizeX)
            {
                iX = selX[i];
                seqptr = x.ptr[iX];
                seqnchar = x.nchar[iX];

                if (annX.length > 0)
                    annptr = annX.ptr[iX];
            }
            else
            {
                indexY++;
                iY = selY[indexY];

                seqptr = y.ptr[iY];
                seqnchar = y.nchar[iY];

                if (annY.length > 0)
                    annptr = annY.ptr[iY];
            }

            freeNode = 1;

            //create new tree for sequence
            kv = createTreeGappy((const char *) seqptr, seqnchar, annptr, k, m, &annotationIndexMap,
                                 presence, reverseComplement, pTree, maxNoOfNodes, &freeNode,
                                 &nullBlock, &printWarning, alphaInf);

            if (kv == NA_REAL)
            {
                for (j = 0; j < sizeX; j++)
                    for (l = 0; l < sizeY; l++)
                        km(j,l) = NA_REAL;

                return;
            }

            if (normalized)
                normValues[i] = sqrt(kv);
            else
                normValues[i] = kv;

            R_CheckUserInterrupt();

            // walk tree to generate sparse feature vector
            currBlock = 0;
            currIndex = 0;
            currStack = -1;
            elemIndex = i * maxNumFeatures;
            featureIndex = 0;
            annotIndex = 0;

            // for situation with no features set sentinel
            featVectors.idx[elemIndex] = maxUnSignedIndex;
            featVectors.x[elemIndex] = MAXINT32;

            while (currStack >= 0 || currIndex <= maxBlockIndex)
            {
                if (pTree->node[currBlock].ib.idx[currIndex] > 0)
                {
                    newBlock = pTree->node[currBlock].ib.idx[currIndex];

                    if (pTree->node[newBlock].leaf & TF_LEAF)
                    {
                        if (!(pTree->node[newBlock].leaf & TF_REV_COMPLEMENT))
                        {
                            for (j = 0; j <= m; j++)
                            {
                                if (pTree->node[newBlock].ib.idx[j] > 0)
                                {
                                    featVectors.x[elemIndex] = pTree->node[newBlock].ib.idx[j];

                                    if (currStack > 4 * k)
                                    {
                                        featVectors.idx[elemIndex++] =
                                        (featureIndex * numAnnPow2K + annotIndex * annCharset.nchar[0] + currIndex) * (m + 1) + j;
                                    }
                                    else
                                    {
                                        featVectors.idx[elemIndex++] =
                                        (featureIndex * alphaInf->numAlphabetChars + currIndex) * (m + 1) + j;
                                    }
                                }
                            }
                        }

                        currIndex++;

                        if (currIndex > maxBlockIndex)
                        {
                            if (currStack == -1)
                                continue;

                            while (currStack >= 0 && currIndex > maxBlockIndex)
                            {
                                currIndex = stack[currStack--];
                                currBlock = stack[currStack--];

                                if (currStack >= 4 * k - 2)
                                    annotIndex = (annotIndex - (currIndex - 1)) / annCharset.nchar[0];
                                else
                                    featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                            }
                        }
                    }
                    else
                    {
                        stack[++currStack] = currBlock;
                        stack[++currStack] = currIndex + 1;

                        // $$$ TODO remove check
                        if (currStack >= 8 * k)
                        {
                            Rprintf("Overflow of tree traversal stack\n");
                            return;
                        }

                        if (currStack > 4 * k)
                            annotIndex = annotIndex * annCharset.nchar[0] + currIndex;
                        else
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

                            if (currStack >= 4 * k - 2)
                                annotIndex = (annotIndex - (currIndex - 1)) / annCharset.nchar[0];
                            else
                                featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                        }
                    }
                }
            }

            featVectors.idx[elemIndex] = maxUnSignedIndex;
            featVectors.x[elemIndex] = MAXINT32;
        }

        computeKernelMatrix(maxUnSignedIndex, featVectors.idx, featVectors.x, km, normValues,
                            maxNumFeatures, sizeX, sizeY, normalized, symmetric);
    }

    return;
}

template<typename T>
void getKMPosDistGappy(T maxUnSignedIndex, NumericMatrix km, ByteStringVector x, ByteStringVector y,
                       int sizeX, int sizeY, IntegerVector selX, IntegerVector selY,
                       IntegerVector offsetX, IntegerVector offsetY, int k, int m,
                       bool normalized,  bool symmetric, bool reverseComplement, bool posSpec,
                       NumericVector distWeight, int maxSeqLength, struct alphaInfo *alphaInf)
{
    T prevIndex;
    uint32_t km1, inew, elemIndex, tempIndex1, tempIndex2, fIndex;
    uint64_t featureIndex, fDimArray, invalidCharPos, numAlphaPowerK, numAlphaPowerK1, *featVectorsStart;
    int i, j, l, offset, maxNumPatternsPerPos, maxFeaturesPerSample, noSamples = sizeX;
    int patLength, seqnchar, index, iold, iX, iY, upper;
    const char *seqptr;
    double kv, *normValues;


    struct featVectorArrays
    {
        int32_t    *x;        // value array
        T          *idx;      // index array
    };

    if (!symmetric)
        noSamples += sizeY;

    km1 = k + m + 1;
    uint64_t *oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    uint64_t *newIndex = (uint64_t *) R_alloc(km1, sizeof(uint64_t));
    numAlphaPowerK1 = ipow64(alphaInf->numAlphabetChars, k - 1);
    numAlphaPowerK = numAlphaPowerK1 * alphaInf->numAlphabetChars;

    // allocate arrays for feature vectors with 8, 16, 32 or 64 bit index
    fDimArray = (maxSeqLength - 2*k - 0.5*m + 1) * (m + 1);
    featVectorArrays featVectors;

    featVectors.x   = (int32_t *) R_alloc(noSamples * fDimArray, sizeof(int32_t));
    featVectors.idx = (T *) R_alloc(noSamples * fDimArray, sizeof(T));
    featVectorsStart = (uint64_t *) R_alloc(noSamples + 1, sizeof(uint64_t));
    normValues = (double *) R_alloc(noSamples, sizeof(double));
    maxNumPatternsPerPos = m + 1;

    featVectorsStart[0] = 0;
    maxFeaturesPerSample = 0;
    elemIndex = 0;
    
    // walk along sequence and create sparse feature vector
    for (i = 0; i < noSamples; i++)
    {
        R_CheckUserInterrupt();

        offset = 0;

        if (i < sizeX)
        {
            iX = selX[i];
            seqptr = x.ptr[iX];
            seqnchar = x.nchar[iX];

            if (offsetX.length() > 0)
                offset = offsetX[iX];
        }
        else
        {
            iY = selY[i - sizeX];
            seqptr = y.ptr[iY];
            seqnchar = y.nchar[iY];

            if (offsetY.length() > 0)
                offset = offsetY[iY];
        }

        patLength = 0;
        featureIndex = 0;
        iold = 0;
        inew = 0;
        kv = 0;
        invalidCharPos = MAXUINT64;

        for (l = 0; l < seqnchar - 2 * k + 1; l++)
        {
            if (invalidCharPos != MAXUINT64 && l <= (int)invalidCharPos)
                continue;
            else
                invalidCharPos = MAXUINT64;

            // look ahead and fill newIndex buffer to process
            // each of the m + 1 patterns at a position immediately
            if (seqnchar - l >  2 * k + m - 1)
                upper = 2 * k + m - 1;
            else
                upper = seqnchar - l;

            if (patLength < 2*k + m - 1)
            {
                for (j = 0; j < upper; j++)
                {
                    if (invalidCharPos != MAXUINT64)
                    {
                        newIndex[inew++] = MAXUINT64;

                        if (inew == km1)
                            inew = 0;

                        continue;
                    }

                    index = alphaInf->seqIndexMap[(int) seqptr[l + j]];

                    if (index > -1)
                    {
                        if (patLength < k)
                        {
                            oldIndex[iold++] = index  * numAlphaPowerK1;

                            if (iold == k)
                                iold = 0;

                            featureIndex = featureIndex * alphaInf->numAlphabetChars + index;
                            patLength++;

                            if (patLength == k)
                            {
                                newIndex[inew++] = featureIndex;

                                if (inew == km1)
                                    inew = 0;
                            }
                        }
                        else
                        {
                            prevIndex = oldIndex[iold];
                            oldIndex[iold++] = index  * numAlphaPowerK1;

                            if (iold == k)
                                iold = 0;

                            featureIndex = (featureIndex - prevIndex) * alphaInf->numAlphabetChars + index;
                            patLength++;

                            newIndex[inew++] = featureIndex;

                            if (inew == km1)
                                inew = 0;
                        }
                    }
                    else
                    {
                        //store max index to buffer for not valid pattern
                        newIndex[inew++] = MAXUINT64;

                        if (inew == km1)
                            inew = 0;

                        invalidCharPos = l + j;
                        featureIndex = 0;
                        patLength = 0;
                        iold = 0;
                        break;
                    }
                }
            }

            if (l + 2 * k + m - 1 < seqnchar && invalidCharPos == MAXUINT64)
            {
                //calc next look ahead feature index and store to buffer
                index = alphaInf->seqIndexMap[(int) seqptr[l + 2 * k + m - 1]];

                if (index > -1)
                {
                    prevIndex = oldIndex[iold];
                    oldIndex[iold++] = index  * numAlphaPowerK1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = (featureIndex - prevIndex) * alphaInf->numAlphabetChars + index;

                    newIndex[inew++] = featureIndex;

                    if (inew == km1)
                        inew = 0;
                }
                else
                {
                    //store max index to buffer for not valid pattern
                    newIndex[inew++] = MAXUINT64;

                    if (inew == km1)
                        inew = 0;

                    invalidCharPos = l + 2 * k + m - 1;
                    featureIndex = 0;
                    patLength = 0;
                    iold = 0;
                }
            }
            else
            {
                inew++;

                if (inew == km1)
                    inew = 0;
            }

            if (seqnchar - 2 * k - l >= m )
                upper = m;
            else
                upper = seqnchar - 2 * k - l;

            if (invalidCharPos != MAXUINT64 && upper > (int) invalidCharPos - l - 2 * k)
                upper = invalidCharPos - l - 2 * k;

            for (j = 0; j <= upper; j++)
            {
                if (newIndex[inew] == MAXUINT64 ||
                    newIndex[(inew + k + j) % km1] == MAXUINT64)
                    continue;

                if (reverseComplement)
                {
                    tempIndex1 = newIndex[inew] * numAlphaPowerK + newIndex[(inew + k + j) % km1];
                    tempIndex2 = tempIndex1;
                    fIndex = 0;

                    for (l = 0; l < 2 * k; l++)
                    {
                        fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                 tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                        tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                    }

                    if (tempIndex1 < fIndex)
                        fIndex = tempIndex1;

                    fIndex = fIndex * (m + 1) + j;
                }
                else
                {
                    fIndex = (newIndex[inew] * numAlphaPowerK +
                              newIndex[(inew + k + j) % km1]) * (m + 1) + j;
                }

                featVectors.idx[elemIndex] = fIndex;
                featVectors.x[elemIndex++] = l - offset;
                kv++;
            }

            if (invalidCharPos != MAXUINT64)
                inew = 0;
        }

        featVectorsStart[i + 1] = elemIndex;

        if (maxFeaturesPerSample < (int) (featVectorsStart[i + 1] - featVectorsStart[i]))
            maxFeaturesPerSample = featVectorsStart[i + 1] - featVectorsStart[i];

        // for dist weighting the kernel value is determined in computeKernelMatrixPos
        if (distWeight.length() == 0)
        {
            if (normalized)
                normValues[i] = sqrt(kv);
            else
                normValues[i] = kv;
        }
    }

    computeKernelMatrixPos(maxUnSignedIndex, featVectors.idx, featVectors.x, featVectorsStart, km,
                           normValues, maxFeaturesPerSample, maxNumPatternsPerPos, sizeX, sizeY,
                           normalized, symmetric, FALSE, distWeight);

    return;
}

static bool getIndexMap(ByteStringVector x, int sizeX, IntegerVector selX, ByteStringVector annCharset,
                        ByteStringVector annX, IntegerVector annotationIndexMap, IntegerVector reverseAnnotationMap,
                        int k, int m, bool normalized, bool presence, bool reverseComplement,
                        struct alphaInfo *alphaInf, ByteStringVector features, uint64_t dimFeatureSpace,
                        bool zeroFeatures, bool useHash, void **indexMap, uint64_t *numUsedFeatures,
                        bool countNonzeroFeatures, uint64_t *numNonzeroFeatures, double **normValues)
{
    int i, j, j1, l, iold, inew, iX, patLength, index, indexAnn, result, km1, upper, patternLength;
    uint32_t *featIndexMap, *featureCounts, featureIndex32, currCount, upperLimit;
    uint64_t numAlphaPowK_1, numAlphaPowK, numAnnPowK_1, numAnnPowK, numAnnPow2K, *featuresInSample;
    uint64_t featureIndex, fIndex, *oldIndex, *newIndex, *oldAnnIndex, *newAnnIndex;
    uint64_t prevIndex, annotIndex, prevAnnIndex, tempIndex1, tempIndex2;
    double kernelValue;
    bool calcKernelValue;
    khiter_t iter;
    khash_t(fim) *hmap;
    khash_t(fc) *fchmap;

    // cyclic buffers for old index
    oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    newIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    numAlphaPowK_1 = ipow64(alphaInf->numAlphabetChars, k - 1);
    numAlphaPowK = numAlphaPowK_1 * alphaInf->numAlphabetChars;
    featuresInSample = NULL;
    featureCounts = NULL;
    featIndexMap = NULL;
    fchmap = NULL;
    hmap = NULL;

    if (annX.length > 0)
    {
        numAnnPowK_1 = ipow64(annCharset.nchar[0], k - 1);
        numAnnPowK = numAnnPowK_1 * annCharset.nchar[0];
        numAnnPow2K = numAnnPowK * numAnnPowK;
        oldAnnIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
        newAnnIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    }
    else
    {
        numAnnPowK_1 = 0;
        numAnnPowK = 0;
        numAnnPow2K = 0;
        oldAnnIndex = NULL;
        newAnnIndex = NULL;
    }

    if (!zeroFeatures || features.length > 0 || countNonzeroFeatures)
    {
        if (normalized && !(zeroFeatures && features.length < 1))
        {
            calcKernelValue = TRUE;
            *normValues = (double *) R_Calloc(sizeX, double);
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
                featureCounts = (uint32_t *) R_Calloc(dimFeatureSpace, uint32_t);
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

            if (features.length > 0)
            {
                for (i=0; i < features.length; i++)
                {
                    if (i % USER_INTERRUPT_LIMIT == 0)
                        R_CheckUserInterrupt();

                    featureIndex = 0;
                    j1 = 0;

                    patternLength = features.nchar[i];

                    if (annX.length > 0)
                        patternLength /= 2;

                    for (j=0; j < patternLength; j++)
                    {
                        if (features.ptr[i][j] == '.')
                            j1++;
                        else
                        {
                            featureIndex = featureIndex * alphaInf->numAlphabetChars +
                                           alphaInf->indexMap[(int) features.ptr[i][j]];
                        }
                    }

                    if (reverseComplement)
                    {
                        tempIndex2 = featureIndex;
                        fIndex = 0;

                        for (l = 0; l < 2 * k; l++)
                        {
                            fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                     tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                            tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                        }

                        if (fIndex < featureIndex)
                            featureIndex = fIndex;
                    }

                    if (annX.length > 0)
                    {
                        annotIndex = 0;

                        for (j=0; j < patternLength; j++)
                        {
                            if (features.ptr[i][patternLength + j] != '.')
                            {
                                annotIndex = annotIndex * annCharset.nchar[0] +
                                             annotationIndexMap[features.ptr[i][patternLength + j]];
                            }
                        }

                        featureIndex = featureIndex * numAnnPow2K + annotIndex;
                    }

                    featureIndex = featureIndex * (m + 1) + j1;
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
            }
        }
        else
        {
            featIndexMap = (uint32_t *) R_Calloc(dimFeatureSpace, uint32_t);
            pFeatureMap = featIndexMap;
            *indexMap = (void *) featIndexMap;

            if (features.length < 1)
            {
                for (i=0; i < (int) dimFeatureSpace; i++)
                {
                    if (i % 100000 == 0)
                        R_CheckUserInterrupt();

                    if (reverseComplement)
                    {
                        if (annX.length > 0)
                        {
                            featureIndex = (i - i % m) / (m + 1) / numAnnPowK / numAnnPow2K;
                            tempIndex2 = featureIndex;
                        }
                        else
                        {
                            featureIndex = (i - i % m) / (m + 1);
                            tempIndex2 = featureIndex;
                        }

                        fIndex = 0;

                        for (l = 0; l < 2 * k; l++)
                        {
                            fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                     tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                            tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                        }

                        if (fIndex < featureIndex)
                        {
                            featIndexMap[i] = MAXUINT32;
                            continue;
                        }
                    }

                    featIndexMap[i] = VALID_FEATURE;
                }
            }
            else
            {
                // init each element to MAXUINT32
                memset(featIndexMap, INDEX_MAP_INIT_FF, dimFeatureSpace*sizeof(uint32_t));

                for (i=0; i < features.length; i++)
                {
                    if (i % USER_INTERRUPT_LIMIT == 0)
                        R_CheckUserInterrupt();

                    featureIndex = 0;
                    j1 = 0;

                    patternLength = features.nchar[i];

                    if (annX.length > 0)
                        patternLength /= 2;

                    for (j=0; j < patternLength; j++)
                    {
                        if (features.ptr[i][j] == '.')
                            j1++;
                        else
                        {
                            featureIndex = featureIndex * alphaInf->numAlphabetChars +
                                           alphaInf->indexMap[(int) features.ptr[i][j]];
                        }
                    }

                    if (reverseComplement)
                    {
                        tempIndex2 = featureIndex;
                        fIndex = 0;

                        for (l = 0; l < 2 * k; l++)
                        {
                            fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                     tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                            tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                        }

                        if (fIndex < featureIndex)
                            continue;
                    }

                    if (annX.length > 0)
                    {
                        annotIndex = 0;

                        for (j=0; j < patternLength; j++)
                        {
                            if (features.ptr[i][patternLength + j] != '.')
                            {
                                annotIndex = annotIndex * annCharset.nchar[0] +
                                             annotationIndexMap[features.ptr[i][patternLength + j]];
                            }
                        }

                        featureIndex = featureIndex * numAnnPow2K + annotIndex;
                    }

                    featureIndex = featureIndex * (m + 1) + j1;
                    featIndexMap[featureIndex] = VALID_FEATURE;
                }
            }
        }

        km1 = k + m + 1;

        // run through all sequences once to get the used features
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
            annotIndex = 0;
            patLength = 0;
            iold = 0;
            inew = 0;
            iX = selX[i];
            kernelValue = 0;

            for (j = 0; j < x.nchar[iX]; j++)
            {
                index = alphaInf->seqIndexMap[(int) x.ptr[iX][j]];

                if (index > -1)
                {
                    if (annX.length == 0)
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
                                newIndex[inew++] = featureIndex;

                                if (inew == km1)
                                    inew = 0;
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

                            newIndex[inew++] = featureIndex;

                            if (inew == km1)
                                inew = 0;

                            patLength++;

                            if (patLength >= 2*k )
                            {
                                if (patLength >= 2*k + m)
                                    upper = m;
                                else
                                    upper = patLength - 2*k;

                                for (j1 = 0; j1 <= upper; j1++)
                                {
                                    if (reverseComplement)
                                    {
                                        tempIndex1 = newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                     newIndex[((inew-1+km1)%km1)];
                                        tempIndex2 = tempIndex1;
                                        fIndex = 0;

                                        for (l = 0; l < 2 * k; l++)
                                        {
                                            fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                                tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                            tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                        }

                                        if (tempIndex1 < fIndex)
                                            fIndex = tempIndex1;

                                        fIndex = fIndex * (m + 1) + j1;
                                    }
                                    else
                                    {
                                        fIndex = (newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                  newIndex[((inew-1+km1)%km1)])*(m+1) + j1;
                                    }

                                    if (useHash)
                                    {
                                        if (calcKernelValue)
                                        {
                                            iter = kh_get(fc, fchmap, fIndex);

                                            if (iter != kh_end(fchmap))
                                            {
                                                currCount = kh_value(fchmap, iter);
                                                kh_value(fchmap, iter) = currCount + 1;
                                                kernelValue = kernelValue - currCount * currCount +
                                                              (currCount + 1) * (currCount + 1);
                                            }
                                            else
                                            {
                                                iter = kh_put(fc, fchmap, fIndex, &result);

                                                if (result == -1)
                                                {
                                                    Rprintf("Storage of key %llu in feature count hashmap failed\n",
                                                            fIndex);
                                                    return (FALSE);
                                                }

                                                kh_value(fchmap, iter) = 1;
                                                kernelValue += 1;
                                            }
                                        }

                                        iter = kh_get(fim, hmap, fIndex);

                                        if (iter != kh_end(hmap))
                                        {
                                            if (kh_value(hmap, iter) != ((uint32_t) i))
                                            {
                                                kh_value(hmap, iter) = i;

                                                if (countNonzeroFeatures)
                                                    (*numNonzeroFeatures)++;
                                            }
                                        }
                                        else
                                        {
                                            if (features.length > 0)
                                                continue;

                                            iter = kh_put(fim, hmap, fIndex, &result);

                                            if (result == -1)
                                            {
                                                Rprintf("Storage of key %llu in hashmap failed\n", fIndex);
                                                return (FALSE);
                                            }

                                            kh_value(hmap, iter) = i;

                                            if (countNonzeroFeatures)
                                                (*numNonzeroFeatures)++;
                                        }
                                    }
                                    else
                                    {
                                        if (calcKernelValue)
                                        {
                                            if (getBit(featuresInSample, fIndex) != 0)
                                            {
                                                currCount = featureCounts[fIndex];
                                                featureCounts[fIndex] += 1;
                                                kernelValue = kernelValue - currCount * currCount +
                                                              (currCount + 1) * (currCount + 1);
                                            }
                                            else
                                            {
                                                setBit(featuresInSample, fIndex);
                                                featureCounts[fIndex] = 1;
                                                kernelValue += 1;
                                            }
                                        }

                                        if (featIndexMap[fIndex] <= VALID_FEATURE &&
                                            featIndexMap[fIndex] != ((uint32_t) i))
                                        {
                                            featIndexMap[fIndex] = i;

                                            if (countNonzeroFeatures)
                                                (*numNonzeroFeatures)++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else   // annX.length > 0
                    {
                        indexAnn = annotationIndexMap[annX.ptr[iX][j]];

                        if (patLength < k)
                        {
                            oldIndex[iold] = index * numAlphaPowK_1;
                            oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                            if (iold == k)
                                iold = 0;

                            featureIndex = featureIndex * alphaInf->numAlphabetChars + index;
                            annotIndex = annotIndex * annCharset.nchar[0] + indexAnn;
                            patLength++;

                            if (patLength == k)
                            {
                                newIndex[inew] = featureIndex;
                                newAnnIndex[inew++] = annotIndex;

                                if (inew == km1)
                                    inew = 0;
                            }
                        }
                        else
                        {
                            prevIndex = oldIndex[iold];
                            prevAnnIndex = oldAnnIndex[iold];
                            oldIndex[iold] = index * numAlphaPowK_1;
                            oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                            if (iold == k)
                                iold = 0;

                            featureIndex = (featureIndex - prevIndex) * alphaInf->numAlphabetChars + index;
                            annotIndex = (annotIndex - prevAnnIndex) * annCharset.nchar[0] + indexAnn;

                            newIndex[inew] = featureIndex;
                            newAnnIndex[inew++] = annotIndex;

                            if (inew == km1)
                                inew = 0;

                            patLength++;

                            if (patLength >= 2*k )
                            {
                                if (patLength >= 2*k + m)
                                    upper = m;
                                else
                                    upper = patLength - 2*k;

                                for (j1 = 0; j1 <= upper; j1++)
                                {
                                    if (reverseComplement)
                                    {
                                        tempIndex1 = newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                     newIndex[((inew-1+km1)%km1)];
                                        tempIndex2 = tempIndex1;
                                        fIndex = 0;

                                        for (l = 0; l < 2 * k; l++)
                                        {
                                            fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                                tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                            tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                        }

                                        if (tempIndex1 < fIndex)
                                            fIndex = tempIndex1;

                                        // annotation not reversed
                                        fIndex = (((fIndex / numAlphaPowK) * numAlphaPowK +
                                                   fIndex % numAlphaPowK) * numAnnPow2K +
                                                  newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                                  newAnnIndex[((inew-1+km1)%km1)]) * (m + 1) + j1;
                                    }
                                    else
                                    {
                                        fIndex = ((newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                   newIndex[((inew-1+km1)%km1)]) * numAnnPow2K +
                                                  newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                                  newAnnIndex[((inew-1+km1)%km1)]) * (m + 1) + j1;
                                    }

                                    if (useHash)
                                    {
                                        if (calcKernelValue)
                                        {
                                            iter = kh_get(fc, fchmap, fIndex);

                                            if (iter != kh_end(fchmap))
                                            {
                                                currCount = kh_value(fchmap, iter);
                                                kh_value(fchmap, iter) = currCount + 1;
                                                kernelValue = kernelValue - currCount * currCount +
                                                              (currCount + 1) * (currCount + 1);
                                            }
                                            else
                                            {
                                                iter = kh_put(fc, fchmap, fIndex, &result);

                                                if (result == -1)
                                                {
                                                    Rprintf("Storage of key %llu in feature count hashmap failed\n",
                                                            fIndex);
                                                    return (FALSE);
                                                }

                                                kh_value(fchmap, iter) = 1;
                                                kernelValue += 1;
                                            }
                                        }

                                        iter = kh_get(fim, hmap, fIndex);

                                        if (iter != kh_end(hmap))
                                        {
                                            if (kh_value(hmap, iter) != ((uint32_t) i))
                                            {
                                                kh_value(hmap, iter) = i;

                                                if (countNonzeroFeatures)
                                                    (*numNonzeroFeatures)++;
                                            }
                                        }
                                        else
                                        {
                                            if (features.length > 0)
                                                continue;

                                            iter = kh_put(fim, hmap, fIndex, &result);

                                            if (result == -1)
                                            {
                                                Rprintf("Storage of key %llu in hashmap failed\n", fIndex);
                                                return (FALSE);
                                            }

                                            kh_value(hmap, iter) = i;

                                            if (countNonzeroFeatures)
                                                (*numNonzeroFeatures)++;
                                        }
                                    }
                                    else
                                    {
                                        if (calcKernelValue)
                                        {
                                            if (getBit(featuresInSample, fIndex) != 0)
                                            {
                                                currCount = featureCounts[fIndex];
                                                featureCounts[fIndex] += 1;
                                                kernelValue = kernelValue - currCount * currCount +
                                                (currCount + 1) * (currCount + 1);
                                            }
                                            else
                                            {
                                                setBit(featuresInSample, fIndex);
                                                featureCounts[fIndex] = 1;
                                                kernelValue += 1;
                                            }
                                        }

                                        if (featIndexMap[fIndex] <= VALID_FEATURE &&
                                            featIndexMap[fIndex] != ((uint32_t) i))
                                        {
                                            featIndexMap[fIndex] = i;

                                            if (countNonzeroFeatures)
                                                (*numNonzeroFeatures)++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    patLength = 0;
                    featureIndex = 0;
                    annotIndex = 0;
                    iold = 0;
                    inew = 0;
                }
            }

            if (calcKernelValue)
                (*normValues)[i] = sqrt(kernelValue);
        }
    }

    if (!zeroFeatures || features.length > 0)
    {
        if (!useHash)
        {
            featureIndex32 = 0;

            if (zeroFeatures)
                upperLimit = MAXUINT32;
            else
                upperLimit = VALID_FEATURE;

            for (i=0; i < (int) dimFeatureSpace; i++)
            {
                if (i % USER_INTERRUPT_LIMIT == 0)
                    R_CheckUserInterrupt();

                if (featIndexMap[i] < upperLimit)
                    featIndexMap[i] = featureIndex32++;
                else
                {
                    if (featIndexMap[i] != MAXUINT32)
                        featIndexMap[i] = MAXUINT32;
                }
            }

            *numUsedFeatures = featureIndex32;
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

                ks_mergesort(gappy, noOfEntries, featureIndices, 0);

                // write col index to hash
                for (int i=0; i < (int) *numUsedFeatures; i++)
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
                        return(FALSE);
                    }
                }
            }
        }
    }

    return(TRUE);
}

void getERDGappy(NumericMatrix erd, ByteStringVector x, int sizeX, IntegerVector selX,
                 ByteStringVector annCharset, ByteStringVector annX, int k, int m,
                 bool normalized, bool presence, bool reverseComplement, struct alphaInfo *alphaInf,
                 ByteStringVector features, uint64_t dimFeatureSpace, bool zeroFeatures, bool useHash,
                 bool mapIndex, void *indexMap, uint64_t numUsedFeatures, double *normValues)
{
    uint32_t *featIndexMap;
    uint64_t prevIndex, featureIndex, prevAnnIndex, annotIndex, fIndex, tempIndex1, tempIndex2;
    uint64_t *oldIndex, *newIndex, *oldAnnIndex, *newAnnIndex;
    uint64_t numAlphaPowK, numAlphaPowK_1, numAnnPowK, numAnnPowK_1, numAnnPow2K;
    int i, j, j1, l, patLength, iold, inew, index, indexAnn, iX, km1, upper;
    double kv, temp;
    bool calcKernelValue;
    khiter_t iter;
    khash_t(fim) *hmap;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);

    if (normalized && normValues == NULL)
    {
        normValues = (double *) R_Calloc(sizeX, double);
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

    if (annX.length > 0)
    {
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);
        numAnnPowK_1 = ipow64(annCharset.nchar[0], k - 1);
        numAnnPowK = numAnnPowK_1 * annCharset.nchar[0];
        numAnnPow2K = numAnnPowK * numAnnPowK;
        oldAnnIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
        newAnnIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    }
    else
    {
        numAnnPowK_1 = 0;
        numAnnPowK = 0;
        numAnnPow2K = 0;
        oldAnnIndex = NULL;
        newAnnIndex = NULL;
    }

    // cyclic buffers for old index
    oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    newIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    numAlphaPowK_1 = ipow64(alphaInf->numAlphabetChars, k - 1);
    numAlphaPowK = numAlphaPowK_1 * alphaInf->numAlphabetChars;
    km1 = k + m + 1;

    // generate ER content
    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();

        featureIndex = 0;
        annotIndex = 0;
        patLength = 0;
        iold = 0;
        inew = 0;
        iX = selX[i];
        kv = 0;

        for (j = 0; j < x.nchar[iX]; j++)
        {
            index = alphaInf->seqIndexMap[(int) x.ptr[iX][j]];

            if (index > -1)
            {
                if (annX.length == 0)
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
                            newIndex[inew++] = featureIndex;

                            if (inew == km1)
                                inew = 0;
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

                        newIndex[inew++] = featureIndex;

                        if (inew == km1)
                            inew = 0;

                        patLength++;

                        if (patLength >= 2 * k )
                        {
                            if (patLength >= 2 * k + m)
                                upper = m;
                            else
                                upper = patLength - 2 * k;

                            for (j1 = 0; j1 <= upper; j1++)
                            {
                                if (reverseComplement)
                                {
                                    tempIndex1 = newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                 newIndex[((inew-1+km1)%km1)];
                                    tempIndex2 = tempIndex1;
                                    fIndex = 0;

                                    for (l = 0; l < 2 * k; l++)
                                    {
                                        fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                                 tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                        tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                    }

                                    if (tempIndex1 < fIndex)
                                        fIndex = tempIndex1;

                                    fIndex = fIndex * (m + 1) + j1;
                                }
                                else
                                {
                                    fIndex = (newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                              newIndex[((inew-1+km1)%km1)])*(m+1) + j1;
                                }

                                if (mapIndex)
                                {
                                    if (useHash)
                                    {
                                        iter = kh_get(fim, hmap, fIndex);

                                        if (iter != kh_end(hmap))
                                            fIndex = kh_value(hmap, iter);
                                        else
                                            continue;
                                    }
                                    else
                                    {
                                        fIndex = featIndexMap[fIndex];

                                        if (fIndex == MAXUINT32)
                                            continue;
                                    }
                                }

                                if (presence)
                                {
                                    if (calcKernelValue)
                                        kv = kv - erd(i, fIndex) + 1;

                                    erd(i, fIndex) = 1;
                                }
                                else
                                {
                                    if (calcKernelValue)
                                    {
                                        temp = erd(i, fIndex);
                                        kv = kv - temp * temp + (temp + 1) * (temp + 1);
                                    }

                                    erd(i, fIndex) += 1;
                                }
                            }
                        }
                    }
                }
                else
                {
                    indexAnn = annotationIndexMap[annX.ptr[iX][j]];

                    if (patLength < k)
                    {
                        oldIndex[iold] = index * numAlphaPowK_1;
                        oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = featureIndex * alphaInf->numAlphabetChars + index;
                        annotIndex = annotIndex * annCharset.nchar[0] + indexAnn;
                        patLength++;

                        if (patLength == k)
                        {
                            newIndex[inew] = featureIndex;
                            newAnnIndex[inew++] = annotIndex;

                            if (inew == km1)
                                inew = 0;
                        }
                    }
                    else
                    {
                        prevIndex = oldIndex[iold];
                        prevAnnIndex = oldAnnIndex[iold];
                        oldIndex[iold] = index * numAlphaPowK_1;
                        oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = (featureIndex - prevIndex) *
                                       alphaInf->numAlphabetChars + index;
                        annotIndex = (annotIndex - prevAnnIndex) *
                                     annCharset.nchar[0] + indexAnn;

                        newIndex[inew] = featureIndex;
                        newAnnIndex[inew++] = annotIndex;

                        if (inew == km1)
                            inew = 0;

                        patLength++;


                        if (patLength >= 2*k )
                        {
                            if (patLength >= 2*k + m)
                                upper = m;
                            else
                                upper = patLength - 2*k;

                            for (j1 = 0; j1 <= upper; j1++)
                            {
                                if (reverseComplement)
                                {
                                    tempIndex1 = newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                 newIndex[((inew-1+km1)%km1)];
                                    tempIndex2 = tempIndex1;
                                    fIndex = 0;

                                    for (l = 0; l < 2 * k; l++)
                                    {
                                        fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                                 tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                        tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                    }

                                    if (tempIndex1 < fIndex)
                                        fIndex = tempIndex1;

                                    // annotation not reversed
                                    fIndex = (((fIndex / numAlphaPowK) * numAlphaPowK +
                                              fIndex % numAlphaPowK) * numAnnPow2K +
                                              newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                              newAnnIndex[((inew-1+km1)%km1)]) * (m + 1) + j1;
                                }
                                else
                                {
                                    fIndex = ((newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                               newIndex[((inew-1+km1)%km1)]) * numAnnPow2K +
                                              newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                              newAnnIndex[((inew-1+km1)%km1)]) * (m + 1) + j1;
                                }

                                if (mapIndex)
                                {
                                    if (useHash)
                                    {
                                        iter = kh_get(fim, hmap, fIndex);

                                        if (iter != kh_end(hmap))
                                            fIndex = kh_value(hmap, iter);
                                        else
                                            continue;
                                    }
                                    else
                                    {
                                        fIndex = featIndexMap[fIndex];

                                        if (fIndex == MAXUINT32)
                                            continue;
                                    }
                                }

                                if (presence)
                                {
                                    if (calcKernelValue)
                                        kv = kv - erd(i, fIndex) + 1;

                                    erd(i, fIndex) = 1;
                                }
                                else
                                {
                                    if (calcKernelValue)
                                    {
                                        temp = erd(i, fIndex);
                                        kv = kv - temp * temp + (temp + 1) * (temp + 1);
                                    }

                                    erd(i, fIndex) += 1;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                patLength = 0;
                featureIndex = 0;
                annotIndex = 0;
                iold = 0;
                inew = 0;
            }
        }

        if (calcKernelValue)
            normValues[i] = sqrt(kv);
    }

    if (normalized)
    {

        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            iX = selX[i];

            for (j = 0; j < (int) numUsedFeatures; j++)
            {
                if (erd(i,j) > 0)
                    erd(i,j) = erd(i,j) / normValues[i];
            }
        }
    }

    return;
}

bool getERSGappy(ByteStringVector x, int sizeX, IntegerVector selX, ByteStringVector annCharset,
                 ByteStringVector annX, int maxSeqLength, int k, int m, bool normalized,
                 bool presence, bool reverseComplement, struct alphaInfo *alphaInf, ByteStringVector features,
                 uint64_t dimFeatureSpace, bool zeroFeatures, bool useHash, bool mapIndex, void *indexMap,
                 uint64_t numUsedFeatures,SEXP slot_p, SEXP slot_j, SEXP slot_x, double *normValues)
{
    uint32_t *featIndexMap;
    uint64_t featureIndex, annotIndex, fIndex, jIdx, nodeLimit;
    uint64_t numAnnPowK, numAnnPow2K, maxNodesPerSequence;
    int i, j, iX, currStack, stack[8*k];
    int freeNode, currBlock, currIndex, newBlock, maxBlockIndex, maxNoOfNodes;
    double kv;
    bool saveKernelValue, printWarning = TRUE;
    const char *annptr;
    khiter_t iter;
    khash_t(fim) *hmap;
    struct prefTree *pTree;
    struct indexBlock nullBlock;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);



    if (normalized && normValues == NULL)
    {
        normValues = (double *) R_Calloc(sizeX, double);
        pNormValues = normValues;
        saveKernelValue = TRUE;
    }
    else
        saveKernelValue = FALSE;

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

    // setup init mask
    for (i = 0; i < MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    if (annX.length > 0)
    {
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);
        numAnnPowK = ipow64(annCharset.nchar[0], k);
        numAnnPow2K = numAnnPowK * numAnnPowK;
        nodeLimit = ((pow(alphaInf->numAlphabetChars, 2 * k + 1) - 1) / (alphaInf->numAlphabetChars - 1)) +
                    pow(alphaInf->numAlphabetChars, 2 * k) *
                    ((pow(annCharset.nchar[0], 2 * k + 1) - 1) / (annCharset.nchar[0] - 1));
        maxNodesPerSequence = (3 * m + 4) * k * (maxSeqLength - 2 * k + 1) + 1;

        if (maxNodesPerSequence < (uint64_t) nodeLimit)
            nodeLimit = (int) maxNodesPerSequence;

        if (alphaInf->maxAlphaIndex > annCharset.nchar[0] - 1)
            maxBlockIndex = alphaInf->maxAlphaIndex;
        else
            maxBlockIndex = annCharset.nchar[0] - 1;

    }
    else
    {
        numAnnPowK = 0;
        numAnnPow2K = 0;
        nodeLimit = (pow(alphaInf->numAlphabetChars, 2 * k + 1) - 1) / (alphaInf->numAlphabetChars - 1);
        maxNodesPerSequence = (m + 2) * k * (maxSeqLength - 2 * k + 1) + 1;

        if (maxNodesPerSequence < (uint64_t) nodeLimit)
            nodeLimit = (int) maxNodesPerSequence;

        maxBlockIndex = alphaInf->maxAlphaIndex;
    }

    // alloc mem for prefix tree
    maxNoOfNodes = MAX_BLOCK;

    if (nodeLimit < (uint64_t) maxNoOfNodes)
        maxNoOfNodes = nodeLimit;

    pTree = (struct prefTree *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    jIdx = 0;


    // fill ER sparse
    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();

        INTEGER(slot_p)[i] = jIdx;

        iX = selX[i];
        annptr = NULL;

        if (annX.length > 0)
            annptr = annX.ptr[iX];

        freeNode = 1;

        //create new tree for sequence
        kv = createTreeGappy((const char *) x.ptr[iX], x.nchar[iX], annptr, k, m, &annotationIndexMap,
                             presence, reverseComplement, pTree, maxNoOfNodes, &freeNode, &nullBlock,
                             &printWarning, alphaInf);

        if (kv == NA_REAL)
            return(FALSE);

        if (saveKernelValue)
            normValues[i] = sqrt(kv);

        R_CheckUserInterrupt();

        // walk tree to generate sparse feature vector and kernel value
        currBlock = 0;
        currIndex = 0;
        currStack = -1;
        featureIndex = 0;
        annotIndex = 0;

        while (currStack >= 0 || currIndex <= maxBlockIndex)
        {
            if (pTree->node[currBlock].ib.idx[currIndex] > 0)
            {
                newBlock = pTree->node[currBlock].ib.idx[currIndex];

                if (pTree->node[newBlock].leaf & TF_LEAF)
                {
                    if (!(pTree->node[newBlock].leaf & TF_REV_COMPLEMENT))
                    {
                        for (j = 0; j <= m; j++)
                        {
                            if (pTree->node[newBlock].ib.idx[j] > 0)
                            {
                                if (currStack > 4 * k)
                                    fIndex = (featureIndex * numAnnPow2K + annotIndex * annCharset.nchar[0] + currIndex) * (m + 1) + j;
                                else
                                    fIndex = (featureIndex * alphaInf->numAlphabetChars + currIndex) * (m + 1) + j;

                                if (features.length > 0)
                                {
                                    if (useHash)
                                    {
                                        iter = kh_get(fim, hmap, fIndex);

                                        if (iter != kh_end(hmap))
                                            INTEGER(slot_j)[jIdx] = kh_value(hmap, iter);
                                        else
                                            continue;
                                    }
                                    else
                                    {
                                        if (featIndexMap[fIndex] < VALID_FEATURE)
                                            INTEGER(slot_j)[jIdx] = featIndexMap[fIndex];
                                        else
                                            continue;
                                    }
                                }
                                else
                                {
                                    if (mapIndex)
                                    {
                                        if (useHash)
                                        {
                                            iter = kh_get(fim, hmap, fIndex);

                                            if (iter != kh_end(hmap))
                                                INTEGER(slot_j)[jIdx] = kh_value(hmap, iter);
                                            else
                                                continue;
                                        }
                                        else
                                            INTEGER(slot_j)[jIdx] = featIndexMap[fIndex];
                                    }
                                    else
                                        INTEGER(slot_j)[jIdx] = fIndex;
                                }

                                if (normalized)
                                {
                                    if (pTree->node[newBlock].ib.idx[j] == 0)
                                        REAL(slot_x)[jIdx] = pTree->node[newBlock].ib.idx[j];
                                    else
                                        REAL(slot_x)[jIdx] = pTree->node[newBlock].ib.idx[j] / normValues[i];
                                }
                                else
                                    REAL(slot_x)[jIdx] = pTree->node[newBlock].ib.idx[j];

                                jIdx++;
                            }
                        }
                    }

                    currIndex++;

                    if (currIndex > maxBlockIndex)
                    {
                        while (currStack >= 0 && currIndex > maxBlockIndex)
                        {
                            currIndex = stack[currStack--];
                            currBlock = stack[currStack--];

                            if (currStack >= 4 * k - 2)
                                annotIndex = (annotIndex - (currIndex - 1)) / annCharset.nchar[0];
                            else
                                featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;

                        }
                    }
                }
                else
                {
                    stack[++currStack] = currBlock;
                    stack[++currStack] = currIndex + 1;

                    // $$$ TODO remove check
                    if (currStack >= 8 * k)
                    {
                        Rprintf("Overflow of tree traversal stack\n");
                        return(FALSE);
                    }

                    if (currStack >= 4 * k)
                        annotIndex = annotIndex * annCharset.nchar[0] + currIndex;
                    else
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
                    while (currStack >= 0 && currIndex > maxBlockIndex)
                    {
                        currIndex = stack[currStack--];
                        currBlock = stack[currStack--];

                        if (currStack >= 4 * k - 2)
                            annotIndex = (annotIndex - (currIndex - 1)) / annCharset.nchar[0];
                        else
                            featureIndex = (featureIndex - (currIndex - 1)) / alphaInf->numAlphabetChars;
                    }
                }
            }
        }
    }

    INTEGER(slot_p)[sizeX] = jIdx;

    return(TRUE);
}

static void assignFeatureNames(SEXP colnames, void *indexMap, int k, int m, struct alphaInfo *alphaInf,
                               uint64_t dimFeatureSpace, ByteStringVector annX, ByteStringVector annCharset,
                               IntegerVector reverseAnnotationMap, bool mapIndex, bool useHash)
{
    int i, j, j1, index, kmerOffset, offset;
    uint32_t *featIndexMap;
    uint64_t *powAnnot, featIndex, currIndex;
    khash_t(fim) *hmap;
    khiter_t iter;

    uint64_t *powAlpha = (uint64_t *) R_alloc(2 * k + 1, sizeof(uint64_t));

    for (i = 0; i <= 2 * k; i++)
        powAlpha[i] = ipow64(alphaInf->numAlphabetChars, i);

    if (annX.length > 0)
    {
        kmerOffset = 2 * k + m;
        powAnnot = (uint64_t *) R_alloc(2 * k + 1, sizeof(uint64_t));

        for (i = 0; i <= 2 * k; i++)
            powAnnot[i] = ipow64(annCharset.nchar[0], i);
    }
    else
    {
        powAnnot = NULL;
        kmerOffset = 0;
    }

    char kmer[2 * k + m + 1 + kmerOffset]; // add char for \0

    if (mapIndex)
    {
        if (useHash)
        {
            hmap = (khash_t(fim) *) indexMap;

            for (iter = kh_begin(hmap); iter != kh_end(hmap); iter++)
            {
                if (kh_exist(hmap, iter))
                {
                    if (kh_value(hmap, iter) < VALID_FEATURE)
                    {
                        currIndex = kh_key(hmap, iter);
                        featIndex = currIndex / (m + 1);
                        offset = currIndex % (m + 1);

                        if (annX.length > 0)
                        {
                            index = 4 * k + 2 * offset - 1;
                            kmer[4 * k + 2 * offset] = '\0';

                            for (j = 0; j < 2 * k; j++)
                            {
                                if (j == k)
                                {
                                    for (j1 = 1; j1 <= offset; j1++)
                                        kmer[index--] = '.';
                                }

                                kmer[index--] =
                                    reverseAnnotationMap[(int)((featIndex % powAnnot[j + 1]) / powAnnot[j])];
                            }

                            featIndex = featIndex / powAnnot[2 * k];
                        }
                        else
                        {
                            index = 2 * k + offset - 1;
                            kmer[2 * k + offset] = '\0';
                        }

                        for (j = 0; j < 2 * k; j++)
                        {
                            if (j == k)
                            {
                                for (j1 = 1; j1 <= offset; j1++)
                                    kmer[index--] = '.';
                            }

                            kmer[index--] =
                            alphaInf->reverseIndexMap[(int)((featIndex % powAlpha[j + 1]) / powAlpha[j])];
                        }

                        SET_STRING_ELT(colnames, kh_value(hmap, iter), Rf_mkChar(kmer));
                    }
                }
            }
        }
        else
        {
            featIndexMap = (uint32_t *) indexMap;

            for (i = 0; i < (int)dimFeatureSpace; i++)
            {
                currIndex = i;

                if (currIndex % USER_INTERRUPT_LIMIT == 0)
                    R_CheckUserInterrupt();

                if (featIndexMap[currIndex] != MAXUINT32)
                {
                    featIndex = currIndex / (m + 1);
                    offset = currIndex % (m + 1);

                    if (annX.length > 0)
                    {
                        index = 4 * k + 2 * offset - 1;
                        kmer[4 * k + 2 * offset] = '\0';

                        for (j = 0; j < 2 * k; j++)
                        {
                            if (j == k)
                            {
                                for (j1 = 1; j1 <= offset; j1++)
                                    kmer[index--] = '.';
                            }

                            kmer[index--] =
                                reverseAnnotationMap[(int)((featIndex % powAnnot[j + 1]) / powAnnot[j])];
                        }

                        featIndex = featIndex / powAnnot[2*k];
                    }
                    else
                    {
                        index = 2 * k + offset - 1;
                        kmer[2 * k + offset] = '\0';
                    }

                    for (j = 0; j < 2 * k; j++)
                    {
                        if (j == k)
                        {
                            for (j1 = 1; j1 <= offset; j1++)
                                kmer[index--] = '.';
                        }

                        kmer[index--] =
                            alphaInf->reverseIndexMap[(int)((featIndex % powAlpha[j + 1]) / powAlpha[j])];
                    }

                    index = featIndexMap[i];
                    SET_STRING_ELT(colnames, index, Rf_mkChar(kmer));
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

            featIndex = i / (m + 1);
            offset = i % (m + 1);

            if (annX.length > 0)
            {
                index = 4 * k + 2 * offset - 1;
                kmer[4 * k + 2 * offset] = '\0';

                for (j = 0; j < 2 * k; j++)
                {
                    if (j == k)
                    {
                        for (j1 = 1; j1 <= offset; j1++)
                            kmer[index--] = '.';
                    }

                    kmer[index--] =
                        reverseAnnotationMap[(int)((featIndex % powAnnot[j + 1]) / powAnnot[j])];
                }

                featIndex = featIndex / powAnnot[2*k];
            }
            else
            {
                index = 2 * k + offset - 1;
                kmer[2 * k + offset] = '\0';
            }

            for (j = 0; j < 2 * k; j++)
            {
                if (j == k)
                {
                    for (j1 = 1; j1 <= offset; j1++)
                        kmer[index--] = '.';
                }

                kmer[index--] =
                    alphaInf->reverseIndexMap[(int)((featIndex % powAlpha[j + 1]) / powAlpha[j])];
            }

            SET_STRING_ELT(colnames, i, Rf_mkChar(kmer));
        }
    }

    return;
}

RcppExport SEXP genExplRepGappyPair(ByteStringVector x, int sizeX, IntegerVector selX,
                                    ByteStringVector annCharset, ByteStringVector annX,
                                    int maxSeqLength, int bioCharset, ByteStringVector features,
                                    int k, int m, bool presence, bool reverseComplement,
                                    bool normalized, bool unmapped, bool lowercase,
                                    bool useRowNames, bool useColNames, bool zeroFeatures,
                                    bool sparse)
{
    int i, numProtect;
    uint64_t dimFeatureSpace, numUsedFeatures, numNonzeroFeatures;
    double *normValues;
    bool mapIndex, useHash;
    const void *vmax;
    void *indexMap;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    SEXP ers, dims, slot_p, slot_j, slot_x, rownames, colnames, dimnames;

    pNormValues = NULL;
    normValues = NULL;
    colnames = NULL;
    pFeatureCounts = NULL;
    pFeatureCountsHMap = NULL;
    pFeatureMap = NULL;
    pFeatureHMap = NULL;
    hmap = NULL;

    // if no features in feature subsetting left
    if (features.length == 0)
        return(generateEmptyExplicitRep(sizeX, sparse));

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    dimFeatureSpace = (m + 1) * ipow64(alphaInf.numAlphabetChars, 2 * k);

    if (annX.length > 0)
    {
        dimFeatureSpace = dimFeatureSpace * ipow64(annCharset.nchar[0], 2 * k);
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);
    }

    // check if features space to large with zero features
    if (dimFeatureSpace > MAXINDEX32 && zeroFeatures && features.length == 0)
    {
        Rprintf("feature space too large with zero features\n");
        return(generateEmptyExplicitRep(sizeX, sparse));
    }

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
            if (!getIndexMap(x, sizeX, selX, annCharset, annX, annotationIndexMap, reverseAnnotationMap,
                             k, m, normalized, presence, reverseComplement, &alphaInf, features, dimFeatureSpace,
                             zeroFeatures, useHash, &indexMap, &numUsedFeatures, FALSE, NULL, &normValues))
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

        vmax = vmaxget();

        if (useColNames)
        {
            assignFeatureNames(colnames, indexMap, k, m, &alphaInf, dimFeatureSpace, annX, annCharset,
                               reverseAnnotationMap, mapIndex, useHash);
        }

        getERDGappy(erd, x, sizeX, selX, annCharset, annX, k, m, normalized, presence,
                    reverseComplement, &alphaInf, features, dimFeatureSpace, zeroFeatures, useHash,
                    mapIndex, indexMap, numUsedFeatures, normValues);

        vmaxset(vmax);

        if (numProtect > 0)
            UNPROTECT(numProtect);

        return(erd);
    }
    else
    {
        numNonzeroFeatures = 0;

        vmax = vmaxget();

        if (!getIndexMap(x, sizeX, selX, annCharset, annX, annotationIndexMap, reverseAnnotationMap,
                         k, m, normalized, presence, reverseComplement, &alphaInf, features, dimFeatureSpace,
                         zeroFeatures, useHash, &indexMap, &numUsedFeatures, TRUE, &numNonzeroFeatures, &normValues))
        {
            vmaxset(vmax);
            return(generateEmptyExplicitRep(sizeX, sparse));
        }

        if (zeroFeatures && features.length == 0)
        {
            indexMap = NULL;
            numUsedFeatures = dimFeatureSpace;
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
        {
            assignFeatureNames(colnames, indexMap, k, m, &alphaInf, dimFeatureSpace, annX,
                               annCharset, reverseAnnotationMap, mapIndex, useHash);
        }

        getERSGappy(x, sizeX, selX, annCharset, annX, maxSeqLength, k, m, normalized,
                    presence, reverseComplement, &alphaInf, features, dimFeatureSpace, zeroFeatures,
                    useHash, mapIndex, indexMap, numUsedFeatures, slot_p, slot_j, slot_x,
                    normValues);

        vmaxset(vmax);

        if (numProtect > 0)
            UNPROTECT(numProtect);

        return(ers);
    }
}

RcppExport SEXP gappyPairKernelMatrixC(SEXP xR, SEXP yR, SEXP selXR, SEXP selYR, SEXP sizeXR,
                    SEXP sizeYR, SEXP isXStringSetR, SEXP symmetricR, SEXP offsetXR,
                    SEXP offsetYR, SEXP annCharsetR, SEXP annXR, SEXP annYR, SEXP bioCharsetR,
                    SEXP ignoreLowerR, SEXP unmappedR, SEXP maxSeqLengthR, SEXP kR, SEXP mR,
                    SEXP posSpecR, SEXP distWeightR, SEXP normalizedR, SEXP presenceR,
                    SEXP reverseComplementR)
{
    int indexSize, sizeX = as<int>(sizeXR);
    int sizeY = as<int>(sizeYR);
    uint8_t maxUIndex8 = MAXUINT8;
    uint16_t maxUIndex16 = MAXUINT16;
    uint32_t maxUIndex32 = MAXUINT32;
    uint64_t dimFeatureSpace, temp, maxUIndex64 = MAXUINT64;
    bool symmetric = as<bool>(symmetricR);
    bool isXStringSet = as<bool>(isXStringSetR);
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    const void *vmax;

    if (symmetric)
        sizeY = sizeX;

    NumericMatrix km(sizeX, sizeY);

    vmax = vmaxget();

    // all R related allocations after this statement
    // to ensure release with vmaxset
    ByteStringVector x, y, annX, annY, annCharset;
    IntegerVector selX(selXR);
    IntegerVector selY(selYR);
    IntegerVector offsetX(offsetXR);
    IntegerVector offsetY(offsetYR);
    NumericVector distWeight(distWeightR);

    if (isXStringSet)
        x = XStringSet2ByteStringVec(xR);
    else
        x = charVector2ByteStringVec(xR);

    annCharset.length = 0;
    annX.length = 0;
    annY.length = 0;

    if (!Rf_isNull(yR))
    {
        if (isXStringSet)
            y = XStringSet2ByteStringVec(yR);
        else
            y = charVector2ByteStringVec(yR);

        if (!Rf_isNull(annYR))
            annY = charVector2ByteStringVec(annYR);
    }
    else
        y.length = 0;

    if (!Rf_isNull(annXR))
    {
        annCharset = charVector2ByteStringVec(annCharsetR);
        annX = charVector2ByteStringVec(annXR);
    }

    int k = as<int>(kR);
    int m = as<int>(mR);
    int bioCharset = as<int>(bioCharsetR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    bool lowercase = !as<bool>(ignoreLowerR);
    bool posSpec = as<bool>(posSpecR);
    bool unmapped = as<bool>(unmappedR);
    bool normalized = as<bool>(normalizedR);
    bool presence = as<bool>(presenceR);
    bool reverseComplement = as<bool>(reverseComplementR);

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    dimFeatureSpace = (m + 1) * ipow64(alphaInf.numAlphabetChars, 2 * k);

    if (annX.length > 0)
        dimFeatureSpace = dimFeatureSpace * ipow64(annCharset.nchar[0], 2 * k);

    indexSize = 1;
    temp = dimFeatureSpace - 1;

    while ((temp >>= 8) > 0)
        indexSize++;

    switch (indexSize)
    {
        case 1:
        {
            if (posSpec || distWeight.length() > 0)
            {
                getKMPosDistGappy(maxUIndex8, km, x, y, sizeX, sizeY, selX, selY, offsetX, offsetY,
                                  k, m, normalized, symmetric, reverseComplement, posSpec,
                                  distWeight, maxSeqLength, &alphaInf);
            }
            else
            {
                getKMStdAnnGappy(maxUIndex8, km, x, y, sizeX, sizeY, selX, selY, annCharset, annX, annY,
                                 k, m, normalized, symmetric, presence, reverseComplement,
                                 maxSeqLength, dimFeatureSpace, &alphaInf);
            }

            break;
        }

        case 2:
        {
            if (posSpec || distWeight.length() > 0)
            {
                getKMPosDistGappy(maxUIndex16, km, x, y, sizeX, sizeY, selX, selY, offsetX, offsetY,
                                  k, m, normalized, symmetric, reverseComplement, posSpec,
                                  distWeight, maxSeqLength, &alphaInf);
            }
            else
            {
                getKMStdAnnGappy(maxUIndex16, km, x, y, sizeX, sizeY, selX, selY, annCharset, annX, annY,
                                 k, m, normalized, symmetric, presence, reverseComplement,
                                 maxSeqLength, dimFeatureSpace, &alphaInf);
            }

            break;
        }

        case 3:
        case 4:
        {
            if (posSpec || distWeight.length() > 0)
            {
                getKMPosDistGappy(maxUIndex32, km, x, y, sizeX, sizeY, selX, selY, offsetX, offsetY,
                                  k, m, normalized, symmetric, reverseComplement, posSpec,
                                  distWeight, maxSeqLength, &alphaInf);
            }
            else
            {
                getKMStdAnnGappy(maxUIndex32, km, x, y, sizeX, sizeY, selX, selY, annCharset, annX, annY,
                                 k, m, normalized, symmetric, presence, reverseComplement,
                                 maxSeqLength, dimFeatureSpace, &alphaInf);
            }

            break;
        }

        default:
        {
            if (posSpec || distWeight.length() > 0)
            {
                getKMPosDistGappy(maxUIndex64, km, x, y, sizeX, sizeY, selX, selY, offsetX, offsetY,
                                  k, m, normalized, symmetric, reverseComplement, posSpec,
                                  distWeight, maxSeqLength, &alphaInf);
            }
            else
            {
                getKMStdAnnGappy(maxUIndex64, km, x, y, sizeX, sizeY, selX, selY, annCharset, annX, annY,
                                 k, m, normalized, symmetric, presence, reverseComplement,
                                 maxSeqLength, dimFeatureSpace, &alphaInf);
            }

            break;
        }
    }

    vmaxset(vmax);

    return(km);
}

uint64_t * featureNamesToIndexGappyPair(SEXP featureNames, int numFeatures, ByteStringVector annCharset,
                                        IntegerVector annotationIndexMap, int k, int m, bool reverseComplement,
                                        struct alphaInfo *alphaInf)
{
    int i, j, l, numDots;
    uint64_t featureIndex, *featIndex, annIndex, tempIndex, fIndex, numAnnPow2K = 0;
    const char *pattern;

    featIndex = (uint64_t *) R_alloc(numFeatures, sizeof(uint64_t));

    if (annCharset.length > 0)
        numAnnPow2K = (uint64_t) ipow64(annCharset.nchar[0], 2 * k);

    for (i = 0; i < numFeatures; i++)
    {
        featureIndex = 0;
        annIndex = 0;

        pattern = CHAR(STRING_ELT(featureNames, i));

        for (j = 0; j < k; j++)
            featureIndex = featureIndex * alphaInf->numAlphabetChars + alphaInf->indexMap[(int) pattern[j]];

        numDots = 0;
        while(pattern[k + numDots] == '.')
            numDots++;

        for (j = k + numDots; j < 2 * k + numDots; j++)
            featureIndex = featureIndex * alphaInf->numAlphabetChars + alphaInf->indexMap[(int) pattern[j]];

        if (reverseComplement)
        {
            tempIndex = featureIndex;
            fIndex = 0;

            for (l = 0; l < 2 * k; l++)
            {
                fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                         tempIndex % (uint64_t) alphaInf->numAlphabetChars - 1;
                tempIndex /= (uint64_t) alphaInf->numAlphabetChars;
            }

            if (fIndex < featureIndex)
                featureIndex = fIndex;
        }

        if (annCharset.length > 0)
        {
            for (j = 2 * k + numDots; j < 3 * k + numDots; j++)
                annIndex = annIndex * annCharset.nchar[0] + annotationIndexMap[pattern[j]];

            for (j = 3 * k + 2 * numDots; j < 4 * k + 2 * numDots; j++)
                annIndex = annIndex * annCharset.nchar[0] + annotationIndexMap[pattern[j]];

            featIndex[i] = (featureIndex * numAnnPow2K + annIndex) * (m + 1) + numDots;

        }
        else
            featIndex[i] = featureIndex * (m + 1) + numDots;
    }

    return(featIndex);
}

RcppExport void featuresToHashmapGappyPair(NumericMatrix featureWeights, int svmIndex, int k, int m,
                                           struct alphaInfo *alphaInf, bool annSpec, ByteStringVector annCharset,
                                           IntegerVector annotationIndexMap)
{
    int i, j, numFeatures, result, numDots;
    uint64_t featureIndex, annotIndex, numAnnPow2K = 0;
    const char *pattern;
    SEXP dimNames, colNames;
    khiter_t iter;
    struct hmData featureData;

    hmap = kh_init(fw);

    numFeatures = featureWeights.ncol();

    if (annSpec)
        numAnnPow2K = (uint64_t) ipow64(annCharset.nchar[0], 2 * k);

    dimNames = Rf_getAttrib(featureWeights, R_DimNamesSymbol);
    colNames = VECTOR_ELT(dimNames, 1);

    for (i = 0; i < numFeatures; i++)
    {
        pattern = CHAR(STRING_ELT(colNames, i));
        featureIndex = 0;
        numDots = 0;

        for (j = 0; j < k; j++)
            featureIndex = featureIndex * alphaInf->numAlphabetChars + alphaInf->indexMap[(int) pattern[j]];

        j = k;

        while (pattern[j] == '.' && j <= (k + m))
        {
            numDots++;
            j++;
        }

        for (j = k + numDots; j < (2*k + numDots); j++)
            featureIndex = featureIndex * alphaInf->numAlphabetChars + alphaInf->indexMap[(int) pattern[j]];

        if (annSpec)
        {
            annotIndex = 0;

            for (j = 2 * k + numDots; j < 3 * k + numDots; j++)
                annotIndex = annotIndex * annCharset.nchar[0] + annotationIndexMap[(int) pattern[j]];

            for (j = 3 * k + 2 * numDots; j < 4 * k + 2 * numDots; j++)
                annotIndex = annotIndex * annCharset.nchar[0] + annotationIndexMap[(int) pattern[j]];

            featureIndex = (featureIndex * numAnnPow2K + annotIndex) * (m + 1) + numDots;
        }
        else
            featureIndex = featureIndex * (m + 1) + numDots;

        iter = kh_put(fw, hmap, featureIndex, &result);

        if (result == -1)
        {
            Rprintf("Storage of key %llu in hashmap failed\n", featureIndex);
            return;
        }

        featureData.unweightedPosIndex = numDots;
        featureData.featureIndex = i;
        kh_value(hmap, iter) = featureData;
    }

    // dealloc of C heap needs to be done by on.exit hook

    return;
}

// delivers the unnormalized feature vectors in compressed format together with kernel values
template<typename T>
void genFeatVectorsPosDepGappyPairT(T maxUnSignedIndex, ByteStringVector x, int sizeX, IntegerVector selX,
                                    IntegerVector offsetX, ByteStringVector annX, ByteStringVector annCharset,
                                    int maxSeqLength, int k, int m, struct alphaInfo *alphaInf, bool presence,
                                    bool normalized, bool reverseComplement, bool posSpecific, NumericVector distWeight,
                                    int sortType, uint64_t **startIndex, T **featVectorIndex,
                                    int32_t **featVectorValue, double **kernelValue)
{
    T prevIndex, featureIndex, tempIndex1, tempIndex2, fIndex;
    int i, j, j1, l, index, indexAnn, iold, inew, patternLength, offset, km1, upper;
    uint64_t elemIndex, numAlphaPowK_1, numAlphaPowK, numAnnPowK_1, numAnnPowK, numAnnPow2K;
    uint64_t *oldIndex, *newIndex, *oldAnnIndex, *newAnnIndex, prevAnnIndex, annotIndex;
    double kv;
    char *seqptr, *annptr;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    IntegerVector selCurr(1), selY(0), offsetY(0);
    NumericMatrix kmOne(1,1);
    ByteStringVector y;
    
    y.length = 0;
    *featVectorIndex = (T *) R_alloc(sizeX * maxSeqLength * (m + 1), sizeof(T));
    *featVectorValue = (int32_t *) R_alloc(sizeX * maxSeqLength * (m + 1), sizeof(int32_t));
    *startIndex = (uint64_t *) R_alloc(sizeX + 1, sizeof(uint64_t));

    if (normalized)
        *kernelValue = (double *) R_alloc(sizeX, sizeof(double));

    if (annX.length > 0)
    {
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);
        numAnnPowK_1 = ipow64(annCharset.nchar[0], k - 1);
        numAnnPowK = numAnnPowK_1 * annCharset.nchar[0];
        numAnnPow2K = numAnnPowK * numAnnPowK;
        oldAnnIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
        newAnnIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    }
    else
    {
        numAnnPowK_1 = 0;
        numAnnPowK = 0;
        numAnnPow2K = 0;
        oldAnnIndex = NULL;
        newAnnIndex = NULL;
        annptr = NULL;
    }

    // cyclic buffers for old index and new index
    oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    newIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    numAlphaPowK_1 = ipow64(alphaInf->numAlphabetChars, k - 1);
    numAlphaPowK = numAlphaPowK_1 * alphaInf->numAlphabetChars;
    km1 = k + m + 1;

    elemIndex = 0;

    for (i = 0; i < sizeX; i++)
    {
        (*startIndex)[i] = elemIndex;
        seqptr = x.ptr[selX[i]];

        if (annX.length > 0)
            annptr = annX.ptr[selX[i]];

        patternLength = 0;
        featureIndex = 0;
        annotIndex = 0;
        iold = 0;
        inew = 0;
        kv = 0;

        if (offsetX.length() > 0)
            offset = offsetX[selX[i]];
        else
            offset = 0;

        // $$$ TODO rewrite loop to avoid sorting of feature vector according to position
        //          by generating the features per position
        for (j = 0; j < x.nchar[selX[i]]; j++)
        {
            index = alphaInf->seqIndexMap[(int)seqptr[j]];

            if (index > -1)
            {
                if (annX.length == 0)
                {
                    if (patternLength < k)
                    {
                        oldIndex[iold++] = index * numAlphaPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = featureIndex * alphaInf->numAlphabetChars + index;
                        patternLength++;

                        if (patternLength == k)
                        {
                            newIndex[inew++] = featureIndex;

                            if (inew == km1)
                                inew = 0;
                        }
                    }
                    else
                    {
                        prevIndex = oldIndex[iold];
                        oldIndex[iold++] = index * numAlphaPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = (featureIndex - prevIndex) * alphaInf->numAlphabetChars + index;
                        newIndex[inew++] = featureIndex;

                        if (inew == km1)
                            inew = 0;

                        patternLength++;

                        if (patternLength >= 2 * k )
                        {
                            if (patternLength >= 2 * k + m)
                                upper = m;
                            else
                                upper = patternLength - 2 * k;

                            for (j1 = upper; j1 >= 0; j1--)
                            {
                                if (reverseComplement)
                                {
                                    tempIndex1 = (newIndex[(inew-1-k-j1+km1)%km1] * numAlphaPowK +
                                                  newIndex[(inew-1+km1)%km1]);
                                    tempIndex2 = tempIndex1;
                                    fIndex = 0;

                                    for (l = 0; l < 2 * k; l++)
                                    {
                                        fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                                 tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                        tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                    }

                                    if (tempIndex1 < fIndex)
                                        fIndex = tempIndex1;

                                    fIndex = fIndex * (m + 1) + j1;
                                }
                                else
                                {
                                    fIndex = (newIndex[(inew-1-k-j1+km1)%km1] * numAlphaPowK +
                                              newIndex[(inew-1+km1)%km1])*(m+1) + j1;
                                }

                                (*featVectorIndex)[elemIndex] = fIndex;
                                (*featVectorValue)[elemIndex++] = j + 1 - 2 * k - j1 - offset + 1;
                                kv++;
                            }
                        }
                    }
                }
                else
                {
                    indexAnn = annotationIndexMap[annptr[j]];

                    if (patternLength < k)
                    {
                        oldIndex[iold] = index * numAlphaPowK_1;
                        oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = featureIndex * alphaInf->numAlphabetChars + index;
                        annotIndex = annotIndex * annCharset.nchar[0] + indexAnn;
                        patternLength++;

                        if (patternLength == k)
                        {
                            newIndex[inew] = featureIndex;
                            newAnnIndex[inew++] = annotIndex;

                            if (inew == km1)
                                inew = 0;
                        }
                    }
                    else
                    {
                        prevIndex = oldIndex[iold];
                        prevAnnIndex = oldAnnIndex[iold];
                        oldIndex[iold] = index * numAlphaPowK_1;
                        oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = (featureIndex - prevIndex) * alphaInf->numAlphabetChars + index;
                        annotIndex = (annotIndex - prevAnnIndex) * annCharset.nchar[0] + indexAnn;

                        newIndex[inew] = featureIndex;
                        newAnnIndex[inew++] = annotIndex;

                        if (inew == km1)
                            inew = 0;

                        patternLength++;

                        if (patternLength >= 2*k )
                        {
                            if (patternLength >= 2*k + m)
                                upper = m;
                            else
                                upper = patternLength - 2*k;

                            for (j1 = upper; j1 >= 0; j1--)
                            {
                                if (reverseComplement)
                                {
                                    tempIndex1 = newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                 newIndex[((inew-1+km1)%km1)];
                                    tempIndex2 = tempIndex1;
                                    fIndex = 0;

                                    for (l = 0; l < 2 * k; l++)
                                    {
                                        fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                                 tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                        tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                    }

                                    if (tempIndex1 < fIndex)
                                        fIndex = tempIndex1;

                                    // annotation not reversed
                                    fIndex = (((fIndex / numAlphaPowK) * numAlphaPowK +
                                               fIndex % numAlphaPowK) * numAnnPow2K +
                                              newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                              newAnnIndex[((inew-1+km1)%km1)]) * (m + 1) + j1;
                                }
                                else
                                {
                                    fIndex = ((newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                               newIndex[((inew-1+km1)%km1)]) * numAnnPow2K +
                                              newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                              newAnnIndex[((inew-1+km1)%km1)]) * (m + 1) + j1;
                                }

                                (*featVectorIndex)[elemIndex] = fIndex;
                                (*featVectorValue)[elemIndex++] = j + 1 - 2 * k - j1 - offset + 1;
                                kv++;
                            }
                        }
                    }
                }
            }
            else
            {
                patternLength = 0;
                featureIndex = 0;
                annotIndex = 0;
                iold = 0;
                inew = 0;
            }
        }

        if (normalized)
        {
            if (distWeight.length() == 0)
                (*kernelValue)[i] = kv;
            else
            {
                selCurr[0] = selX[i];
                int currSeqLength = x.nchar[selX[i]];
                
                getKMPosDistGappy(maxUnSignedIndex, kmOne, x, y, 1, 1, selCurr, selY, offsetX, offsetY,
                                  k, m, FALSE, TRUE, reverseComplement, posSpecific, distWeight,
                                  currSeqLength, alphaInf);
                (*kernelValue)[i] = kmOne(0,0);
                
            }
        }
    }

    (*startIndex)[sizeX] = elemIndex;

    if (sortType == KBS_SORT_BY_POSITION)
    {
        sort2Arrays((int32_t) (MAXUINT32 >> 1), *featVectorValue, *featVectorIndex, 1,
                    (int) elemIndex, *startIndex);
    }

    return;
}

void genFeatVectorsGappyPair(ByteStringVector x, int sizeX, IntegerVector selX, IntegerVector offsetX,
                             ByteStringVector annX, ByteStringVector annCharset, int maxSeqLength,
                             int k, int m, struct alphaInfo *alphaInf, uint64_t dimFeatureSpace,
                             bool presence, bool normalized, bool reverseComplement, bool posSpecific,
                             NumericVector distWeight, int sortType, int numPositions, uint64_t **startIndex,
                             void **featVectorIndex, int32_t **featVectorValue, double **kernelValue,
                             int *indexSize)
{
    uint8_t maxUIndex8 = MAXUINT8;
    uint16_t maxUIndex16 = MAXUINT16;
    uint32_t maxUIndex32 = MAXUINT32;
    uint64_t maxUIndex64 = MAXUINT64;
    uint64_t temp;

    if (sortType == KBS_UNSORTED)
        temp = dimFeatureSpace;
    else
        temp = dimFeatureSpace * numPositions - 1;

    *indexSize = 1;

    while ((temp >>= 8) > 0)
        (*indexSize)++;

    switch (*indexSize)
    {
        case 1:
        {
            genFeatVectorsPosDepGappyPairT(maxUIndex8, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
                                           k, m, alphaInf, presence, normalized, reverseComplement, posSpecific,
                                           distWeight, sortType, startIndex, (uint8_t **) featVectorIndex,
                                           featVectorValue, kernelValue);
            return;
        }

        case 2:
        {
            genFeatVectorsPosDepGappyPairT(maxUIndex16, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
                                           k, m, alphaInf, presence, normalized, reverseComplement, posSpecific,
                                           distWeight, sortType, startIndex, (uint16_t **) featVectorIndex,
                                           featVectorValue, kernelValue);
            return;
        }

        case 3:
        case 4:
        {
            genFeatVectorsPosDepGappyPairT(maxUIndex32, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
                                           k, m, alphaInf, presence, normalized, reverseComplement, posSpecific,
                                           distWeight, sortType, startIndex, (uint32_t **) featVectorIndex,
                                           featVectorValue, kernelValue);
            return;
        }

        default:
        {
            genFeatVectorsPosDepGappyPairT(maxUIndex64, x, sizeX, selX, offsetX, annX, annCharset, maxSeqLength,
                                           k, m, alphaInf, presence, normalized, reverseComplement, posSpecific,
                                           distWeight, sortType, startIndex, (uint64_t **) featVectorIndex,
                                           featVectorValue, kernelValue);
            return;
        }
    }
}

template<typename T>
bool getSVWeightsFeatGappyPair(T maxUnSignedIndex, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, ByteStringVector x,
                               int sizeX, IntegerVector selX, IntegerVector offsetX, int maxSeqLength, NumericVector coefs,
                               bool reverseComplement, bool posSpecific, NumericVector distWeight, double weightLimit,
                               int k, int m, int minPos, int maxPos, struct alphaInfo *alphaInf, bool normalized,
                               uint64_t *numKeys, T **keys)
{
    T prevIndex, featureIndex, tempIndex1, tempIndex2;
    int i, j, j1, l, iX, index, iold, inew, offset, result, patternLength, upper, km1, currSeqLength;
    uint64_t numAlphaPowK_1, numAlphaPowK, numAlphaPow2K, numEntries, key, fIndex, *oldIndex, *newIndex;
    double kv, limit, normFactor;
    char *seqptr;
    khiter_t iter;
    IntegerVector selCurr(1), selY(0), offsetY(0);
    NumericMatrix kmOne(1,1);
    ByteStringVector y;
    numAlphaPowK_1 = ipow64(alphaInf->numAlphabetChars, k - 1);
    numAlphaPowK = numAlphaPowK_1 * alphaInf->numAlphabetChars;
    numAlphaPow2K = numAlphaPowK * numAlphaPowK;
    oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    newIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));

    km1 = k + m + 1;
    normFactor = 1;
    y.length = 0;

    for (i = 0; i < sizeX; i++)
    {
        if (i % USER_INTERRUPT_LIMIT == 0)
            R_CheckUserInterrupt();

        iX = selX[i];
        seqptr = x.ptr[iX];

        if (offsetX.length() > 0)
            offset = offsetX[iX];
        else
            offset = 0;

        kv = 0;
        patternLength = 0;

        // get kernel value
        if (normalized)
        {
            if (distWeight.length() == 0)
            {
                for (j = 0; j < x.nchar[iX]; j++)
                {
                    index = alphaInf->seqIndexMap[(int) seqptr[j]];
                    
                    if (index > -1)
                    {
                        if (patternLength < 2*k + m)
                            patternLength++;
                        
                        if (patternLength >= 2 * k)
                        {
                            if (patternLength == 2 * k + m)
                                upper = m;
                            else
                                upper = patternLength - 2 * k;
                            
                            kv += upper + 1;
                        }
                    }
                    else
                        patternLength = 0;
                }
            }
            else
            {
                selCurr[0] = iX;
                currSeqLength = x.nchar[iX];
                
                getKMPosDistGappy(maxUnSignedIndex, kmOne, x, y, 1, 1, selCurr, selY, offsetX, offsetY, k, m,
                                  FALSE, TRUE, reverseComplement, posSpecific, distWeight, currSeqLength, alphaInf);
                kv = kmOne(0,0);
            }
            
            normFactor = 1.0 / sqrt(kv);
        }

        patternLength = 0;
        featureIndex = 0;
        iold = 0;
        inew = 0;

        for (j = 0; j < x.nchar[iX]; j++)
        {
            index = alphaInf->seqIndexMap[(int) seqptr[j]];

            if (index > -1)
            {
                if (patternLength < k)
                {
                    oldIndex[iold++] = index * numAlphaPowK_1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = featureIndex * alphaInf->numAlphabetChars + index;
                    patternLength++;

                    if (patternLength == k)
                    {
                        newIndex[inew++] = featureIndex;

                        if (inew == km1)
                            inew = 0;
                    }
                }
                else
                {
                    prevIndex = oldIndex[iold];
                    oldIndex[iold++] = index * numAlphaPowK_1;

                    if (iold == k)
                        iold = 0;

                    featureIndex = (featureIndex - prevIndex) * alphaInf->numAlphabetChars + index;
                    newIndex[inew++] = featureIndex;

                    if (inew == km1)
                        inew = 0;

                    patternLength++;

                    if (patternLength >= 2 * k )
                    {
                        if (patternLength >= 2 * k + m)
                            upper = m;
                        else
                            upper = patternLength - 2 * k;

                        for (j1 = 0; j1 <= upper; j1++)
                        {
                            if (reverseComplement)
                            {
                                tempIndex1 = (newIndex[(inew-1-k-j1+km1)%km1] * numAlphaPowK +
                                              newIndex[(inew-1+km1)%km1]);
                                tempIndex2 = tempIndex1;
                                fIndex = 0;

                                for (l = 0; l < 2 * k; l++)
                                {
                                    fIndex = (fIndex + 1) * alphaInf->numAlphabetChars -
                                             tempIndex2 % (uint64_t) alphaInf->numAlphabetChars - 1;
                                    tempIndex2 /= (uint64_t) alphaInf->numAlphabetChars;
                                }

                                if (tempIndex1 < fIndex)
                                    fIndex = tempIndex1;

                                fIndex = fIndex * (m + 1) + j1;
                            }
                            else
                            {
                                fIndex = (newIndex[(inew-1-k-j1+km1)%km1] * numAlphaPowK +
                                          newIndex[(inew-1+km1)%km1])*(m+1) + j1;
                            }

                            iter = kh_put(pdfi, pdfimap, fIndex, &result);
                            key = (j - 2 * k + 1 - j1 - offset - minPos + 1) * (m + 1) * numAlphaPow2K + fIndex;
                            iter = kh_put(pdfw, pdfwmap, key, &result);

                            if (result)
                                kh_value(pdfwmap, iter) = normFactor * coefs[iX];
                            else
                                kh_value(pdfwmap, iter) = kh_value(pdfwmap, iter) + normFactor * coefs[iX];
                        }
                    }
                }
            }
            else
            {
                patternLength = 0;
                featureIndex = 0;
                iold = 0;
                inew = 0;
            }
        }
    }

    *numKeys = kh_size(pdfwmap);

    if (kh_size(pdfwmap) == 0)
        return(TRUE);

    // create mapping to index
    *keys = (T *) R_Calloc(kh_size(pdfimap) + 1, T);
    numEntries = 0;

    for (iter = kh_begin(pdfimap); iter != kh_end(pdfimap); iter++)
    {
        if (kh_exist(pdfimap, iter))
            (*keys)[numEntries++] = kh_key(pdfimap, iter);
    }

    sortArray(maxUnSignedIndex, *keys, 1, numEntries);

    for (i = 0; i < (int) numEntries; i++)
    {
        iter = kh_get(pdfi, pdfimap, (*keys)[i]);
        kh_value(pdfimap, iter) = i;
    }

    // perform weight pruning for pos specific kernels
    // for pos dependent kernels weight pruning is done at prediction or
    // profile generation when the actual weights are calculated
    if (posSpecific)
        limit = 0;
    else
        limit = weightLimit;

    *keys = (T *) R_Calloc(kh_size(pdfwmap), T);
    numEntries = 0;

    for (iter = kh_begin(pdfwmap); iter != kh_end(pdfwmap); iter++)
    {
        if (numEntries % USER_INTERRUPT_LIMIT == 0)
            R_CheckUserInterrupt();

        if (kh_exist(pdfwmap, iter) && (fabs(kh_value(pdfwmap, iter)) > limit))
            (*keys)[numEntries++] = kh_key(pdfwmap, iter);
    }

    // Realloc to shrink size
    if (*numKeys != numEntries)
    {
        *numKeys = numEntries;
        *keys = (T *) R_Realloc(*keys, *numKeys, T);
    }

    // sort keys according to position and feature index
    sortArray(maxUnSignedIndex, *keys, 1, *numKeys);

    return(TRUE);
}

template<typename T>
void getWeightedFeatOfSVGappyPair(T maxUnSignedIndex, SEXP **pdFeatWeights, khash_t(pdfw) *pdfwmap,
                                  khash_t(pdfi) *pdfimap, ByteStringVector x, int sizeX, IntegerVector selX,
                                  IntegerVector offsetX, int maxSeqLength, NumericVector coefs, bool posSpecific,
                                  NumericVector distWeight, bool reverseComplement, double weightLimit, int k,
                                  int m, int minPos, int maxPos, struct alphaInfo *alphaInf, bool normalized,
                                  uint64_t *numKeys, T **keys)
{
    int i, j, j1, index, row, numProtect, numDots;
    char kmer[k + 1], position[12];
    uint64_t featIndex, *powAlpha;
    khiter_t iter;
    SEXP rownames, colnames, dimnames, slot_p, slot_i, slot_x, dims;

    if (!getSVWeightsFeatGappyPair(maxUnSignedIndex, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength, coefs,
                                   reverseComplement, posSpecific, distWeight, weightLimit, k, m, minPos, maxPos, alphaInf,
                                   normalized, numKeys, keys))
    {
        pdFeatWeights = NULL;
        return;
    }

    // alloc dgCMatrix for position dependent feature weights
    numProtect = 0;
    **pdFeatWeights = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    numProtect++;
    dims = PROTECT(Rf_allocVector(INTSXP, 2));
    numProtect++;
    SET_SLOT(**pdFeatWeights, Rf_mkChar(MATRIX_DIM_SLOT), dims);
    INTEGER(dims)[0] = kh_size(pdfimap);
    INTEGER(dims)[1] = maxPos - minPos + 1;

    PROTECT(slot_p = Rf_allocVector(INTSXP, maxPos - minPos + 2));
    PROTECT(slot_i = Rf_allocVector(INTSXP, *numKeys));
    PROTECT(slot_x = Rf_allocVector(REALSXP, *numKeys));
    numProtect += 3;
    SET_SLOT(**pdFeatWeights, Rf_mkChar(MATRIX_I_SLOT), slot_i);
    SET_SLOT(**pdFeatWeights, Rf_mkChar(MATRIX_P_SLOT), slot_p);
    SET_SLOT(**pdFeatWeights, Rf_mkChar(MATRIX_X_SLOT), slot_x);

    PROTECT(rownames = Rf_allocVector(STRSXP, kh_size(pdfimap)));
    PROTECT(colnames = Rf_allocVector(STRSXP, maxPos - minPos + 1));
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    numProtect += 3;
    SET_VECTOR_ELT(dimnames, 0, rownames);
    SET_VECTOR_ELT(dimnames, 1, colnames);
    SET_SLOT(**pdFeatWeights, Rf_mkChar(MATRIX_DIMNAMES_SLOT), dimnames);

    powAlpha = (uint64_t *) R_alloc(2 * k + 1, sizeof(uint64_t));
    powAlpha[0] = 1;

    for (i = 0; i <= 2*k; i++)
        powAlpha[i+1] = powAlpha[i] * alphaInf->numAlphabetChars;

    for (i = minPos; i <= maxPos; i++)
    {
        snprintf(position, 12, "%d", i);
        SET_STRING_ELT(colnames, i - minPos, Rf_mkChar(position));
    }

    kmer[k] = '\0';

    for (iter = kh_begin(pdfimap); iter != kh_end(pdfimap); iter++)
    {
        if (kh_exist(pdfimap, iter))
        {
            row = kh_value(pdfimap, iter);
            featIndex = kh_key(pdfimap, iter);
            numDots = featIndex % (m + 1);
            featIndex = featIndex / (m + 1);
            kmer[2 * k + numDots] = '\0';
            index = 2 * k + numDots - 1;

            for (j = 0; j < 2 * k; j++)
            {
                if (j == k)
                {
                    for (j1 = 1; j1 <= numDots; j1++)
                        kmer[index--] = '.';
                }

                kmer[index--] =
                    alphaInf->reverseIndexMap[(int)((featIndex % powAlpha[j + 1]) / powAlpha[j])];
            }

            SET_STRING_ELT(rownames, row, Rf_mkChar(kmer));
        }
    }

    UNPROTECT(numProtect);

    return;
}

void getFeaturesOfSVGappyPair(SEXP **pdFeatWeights, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, ByteStringVector x,
                              int sizeX, IntegerVector selX, IntegerVector offsetX, int maxSeqLength, NumericVector coefs,
                              bool reverseComplement, bool posSpecific, NumericVector distWeight, double weightLimit, int k,
                              int m, int minPos, int maxPos, uint64_t dimFeatureSpace, struct alphaInfo *alphaInf,
                              bool normalized, int featIndexSize, uint64_t *numKeys, void **keys)
{
    uint8_t maxUIndex8 = MAXUINT8;
    uint16_t maxUIndex16 = MAXUINT16;
    uint32_t maxUIndex32 = MAXUINT32;
    uint64_t maxUIndex64 = MAXUINT64;


    switch (featIndexSize)
    {
        case 1:
        {
            getWeightedFeatOfSVGappyPair(maxUIndex8, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                         coefs, posSpecific, distWeight, reverseComplement, weightLimit, k, m, minPos, maxPos,
                                         alphaInf, normalized, numKeys, (uint8_t **) keys);
            break;
        }

        case 2:
        {
            getWeightedFeatOfSVGappyPair(maxUIndex16, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                         coefs, posSpecific, distWeight, reverseComplement, weightLimit, k, m, minPos, maxPos,
                                         alphaInf, normalized, numKeys, (uint16_t **) keys);
            break;
        }

        case 3:
        case 4:
        {
            getWeightedFeatOfSVGappyPair(maxUIndex32, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                         coefs, posSpecific, distWeight, reverseComplement, weightLimit, k, m, minPos, maxPos,
                                         alphaInf, normalized, numKeys, (uint32_t **) keys);
            break;
        }

        default:
        {
            getWeightedFeatOfSVGappyPair(maxUIndex64, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                         coefs, posSpecific, distWeight, reverseComplement, weightLimit, k, m, minPos, maxPos,
                                         alphaInf, normalized, numKeys, (uint64_t **) keys);
            break;
        }
    }

    return;
}

void genPredProfileGappyPair(NumericMatrix pprof, ByteStringVector x, IntegerVector selX, int numSamples,
                             ByteStringVector annCharset, ByteStringVector annX, int maxSeqLength, bool unmapped,
                             bool reverseComplement, int kernelType, int k, int m, int bioCharset,
                             NumericMatrix featureWeights, int svmIndex, bool lowercase, bool normalized,
                             bool presence)
{
    int i, j, l, j1, km1, iX, iold, inew, patLength, index, indexAnn, upper, result;
    uint32_t currCount;
    uint64_t dimFeatureSpace, featureIndex, featureIndex1, prevIndex, tempIndex1, tempIndex2, fIndex;
    uint64_t annotIndex, prevAnnIndex, *oldAnnIndex, *newAnnIndex, numAnnPowK, numAnnPowK_1, numAnnPow2K;
    double kernelValue, *normValues, featureWeight;
    bool useHash = FALSE, calcKernelValue = FALSE;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    struct hmData featureData;
    khiter_t iter;
    khash_t(fc) *fchmap;

    pNormValues = NULL;
    normValues = NULL;
    pFeatureCounts = NULL;
    pFeatureCountsHMap = NULL;
    pFeatureMap = NULL;
    pFeatureHMap = NULL;
    hmap = NULL;
    fchmap = NULL;

    dimFeatureSpace = (m + 1) * ipow64((uint64_t) alphaInf.numAlphabetChars, 2 * k);

    if (annX.length > 0)
    {
        dimFeatureSpace *= ipow64((uint64_t) annCharset.nchar[0], 2 * k);
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);
        numAnnPowK_1 = ipow64(annCharset.nchar[0], k - 1);
        numAnnPowK = numAnnPowK_1 * annCharset.nchar[0];
        numAnnPow2K = numAnnPowK * numAnnPowK;
        oldAnnIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
        newAnnIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    }
    else
    {
        numAnnPowK_1 = 0;
        numAnnPowK = 0;
        numAnnPow2K = 0;
        oldAnnIndex = NULL;
        newAnnIndex = NULL;
    }

    if (normalized)
    {
        normValues = (double *) R_Calloc(numSamples, double);
        pNormValues = normValues;
        fchmap = kh_init(fc);
        pFeatureCountsHMap = fchmap;
        calcKernelValue = TRUE;
    }

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    // if not all features present use hash map
    if ((uint64_t) featureWeights.ncol() != dimFeatureSpace)
    {
        useHash = TRUE;
        featuresToHashmapGappyPair(featureWeights, svmIndex, k, m, &alphaInf, annX.length > 0,
                                   annCharset, annotationIndexMap);
    }

    // cyclic buffers for old index
    uint64_t *oldIndex = (uint64_t *) R_alloc(k, sizeof(uint64_t));
    uint64_t *newIndex = (uint64_t *) R_alloc(k + m + 1, sizeof(uint64_t));
    uint64_t numAlphaPowK_1 = ipow64(alphaInf.numAlphabetChars, k - 1);
    uint64_t numAlphaPowK = numAlphaPowK_1 * alphaInf.numAlphabetChars;
    km1 = k + m + 1;

    for (i = 0 ; i < numSamples; i++)
    {
        kernelValue = 0;
        featureIndex = 0;
        annotIndex = 0;
        patLength = 0;
        iold = 0;
        inew = 0;
        iX = selX[i];

        for (j = 0; j < x.nchar[iX]; j++)
        {
            index = alphaInf.seqIndexMap[(int) x.ptr[iX][j]];

            if (index > -1)
            {
                if (annX.length == 0)
                {
                    if (patLength < k)
                    {
                        oldIndex[iold++] = index * numAlphaPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = featureIndex * alphaInf.numAlphabetChars + index;
                        patLength++;

                        if (patLength == k)
                        {
                            newIndex[inew++] = featureIndex;

                            if (inew == km1)
                                inew = 0;
                        }
                    }
                    else
                    {
                        prevIndex = oldIndex[iold];
                        oldIndex[iold++] = index * numAlphaPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = (featureIndex - prevIndex) *
                                       alphaInf.numAlphabetChars + index;

                        newIndex[inew++] = featureIndex;

                        if (inew == km1)
                            inew = 0;

                        patLength++;

                        if (patLength >= 2*k )
                        {
                            if (patLength >= 2*k + m)
                                upper = m;
                            else
                                upper = patLength - 2*k;

                            for (j1 = 0; j1 <= upper; j1++)
                            {
                                if (reverseComplement)
                                {
                                    tempIndex1 = newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                    newIndex[((inew-1+km1)%km1)];
                                    tempIndex2 = tempIndex1;
                                    fIndex = 0;

                                    for (l = 0; l < 2 * k; l++)
                                    {
                                        fIndex = (fIndex + 1) * alphaInf.numAlphabetChars -
                                                 tempIndex2 % (uint64_t) alphaInf.numAlphabetChars - 1;
                                        tempIndex2 /= (uint64_t) alphaInf.numAlphabetChars;
                                    }

                                    if (tempIndex1 < fIndex)
                                        fIndex = tempIndex1;

                                    fIndex = fIndex * (m + 1) + j1;
                                }
                                else
                                {
                                    fIndex = (newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                              newIndex[((inew-1+km1)%km1)])*(m+1) + j1;
                                }

                                if (calcKernelValue)
                                {
                                    iter = kh_get(fc, fchmap, fIndex);

                                    if (iter != kh_end(fchmap))
                                    {
                                        if (!presence)
                                        {
                                            currCount = kh_value(fchmap, iter);
                                            kh_value(fchmap, iter) = currCount + 1;
                                            kernelValue = kernelValue - currCount * currCount +
                                                          (currCount + 1) * (currCount + 1);
                                        }
                                    }
                                    else
                                    {
                                        iter = kh_put(fc, fchmap, fIndex, &result);

                                        if (result == -1)
                                        {
                                            Rprintf("Storage of key %llu in feature count hashmap failed\n",
                                                    fIndex);
                                            return;
                                        }

                                        kh_value(fchmap, iter) = 1;
                                        kernelValue += 1;
                                    }
                                }

                                if (useHash)
                                {
                                    // check if new feature
                                    iter = kh_get(fw, hmap, fIndex);

                                    if (iter != kh_end(hmap))
                                    {
                                        featureData = kh_value(hmap, iter);
                                        featureWeight = featureWeights(svmIndex, featureData.featureIndex) / (2 * k);

                                        for (l = j+1-2*k-j1; l < j+1-k-j1; l++)
                                        {
                                            pprof(i, l) += featureWeight;
                                            pprof(i, l + k + j1) += featureWeight;
                                        }
                                    }
                                }
                                else
                                {
                                    for (l = j+1-2*k-j1; l < j+1-k-j1; l++)
                                    {
                                        pprof(i, l) += featureWeights(svmIndex, fIndex) / (2 * k);
                                        pprof(i, l + k + j1) += featureWeights(svmIndex, fIndex) / (2 * k);
                                    }
                                }
                            }
                        }

                    }
                }
                else
                {
                    indexAnn = annotationIndexMap[annX.ptr[iX][j]];

                    if (patLength < k)
                    {
                        oldIndex[iold] = index * numAlphaPowK_1;
                        oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = featureIndex * alphaInf.numAlphabetChars + index;
                        annotIndex = annotIndex * annCharset.nchar[0] + indexAnn;
                        patLength++;

                        if (patLength == k)
                        {
                            newIndex[inew] = featureIndex;
                            newAnnIndex[inew++] = annotIndex;

                            if (inew == km1)
                                inew = 0;
                        }
                    }
                    else
                    {
                        prevIndex = oldIndex[iold];
                        prevAnnIndex = oldAnnIndex[iold];
                        oldIndex[iold] = index * numAlphaPowK_1;
                        oldAnnIndex[iold++] = indexAnn * numAnnPowK_1;

                        if (iold == k)
                            iold = 0;

                        featureIndex = (featureIndex - prevIndex) *
                                       alphaInf.numAlphabetChars + index;

                        annotIndex = (annotIndex - prevAnnIndex) *
                                     annCharset.nchar[0] + indexAnn;

                        newIndex[inew] = featureIndex;
                        newAnnIndex[inew++] = annotIndex;

                        if (inew == km1)
                            inew = 0;

                        patLength++;

                        if (patLength >= 2*k )
                        {
                            if (patLength >= 2*k + m)
                                upper = m;
                            else
                                upper = patLength - 2*k;

                            for (j1 = 0; j1 <= upper; j1++)
                            {
                                featureIndex1 = ((newIndex[((inew-1-k-j1+km1)%km1)] * numAlphaPowK +
                                                  newIndex[((inew-1+km1)%km1)]) * numAnnPow2K +
                                                 newAnnIndex[((inew-1-k-j1+km1)%km1)] * numAnnPowK +
                                                 newAnnIndex[((inew-1+km1)%km1)]) * (m+1) + j1;

                                if (calcKernelValue)
                                {
                                    iter = kh_get(fc, fchmap, featureIndex1);

                                    if (iter != kh_end(fchmap))
                                    {
                                        if (!presence)
                                        {
                                            currCount = kh_value(fchmap, iter);
                                            kh_value(fchmap, iter) = currCount + 1;
                                            kernelValue = kernelValue - currCount * currCount +
                                                          (currCount + 1) * (currCount + 1);
                                        }
                                    }
                                    else
                                    {
                                        iter = kh_put(fc, fchmap, featureIndex1, &result);

                                        if (result == -1)
                                        {
                                            Rprintf("Storage of key %llu in feature count hashmap failed\n",
                                                    featureIndex1);
                                            return;
                                        }

                                        kh_value(fchmap, iter) = 1;
                                        kernelValue += 1;
                                    }
                                }

                                if (useHash)
                                {
                                    // check if new feature
                                    iter = kh_get(fw, hmap, featureIndex1);

                                    if (iter != kh_end(hmap))
                                    {
                                        featureData = kh_value(hmap, iter);
                                        featureWeight = featureWeights(svmIndex, featureData.featureIndex) / (2 * k);

                                        for (l = j+1-2*k-j1; l < j+1-k-j1; l++)
                                        {
                                            pprof(i, l) += featureWeight;
                                            pprof(i, l + k + j1) += featureWeight;
                                        }
                                    }
                                }
                                else
                                {
                                    for (l = j+1-2*k-j1; l < j+1-k-j1; l++)
                                    {
                                        pprof(i, l) += featureWeights(svmIndex, featureIndex1) / (2 * k);
                                        pprof(i, l + k + j1) += featureWeights(svmIndex, featureIndex1) / (2 * k);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                featureIndex = 0;
                annotIndex = 0;
                patLength = 0;
                iold = 0;
                inew = 0;
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

void freeHeapGappyPair()
{
    if (pNormValues != NULL)
    {
        R_Free(pNormValues);
        pNormValues = NULL;
    }

    if (pFeatureCounts != NULL)
    {
        R_Free(pFeatureCounts);
        pFeatureCounts = NULL;
    }

    if (pFeatureCountsHMap != NULL)
    {
        kh_destroy(fc, pFeatureCountsHMap);
        pFeatureCountsHMap = NULL;
    }

    if (pFeatureMap != NULL)
    {
        R_Free(pFeatureMap);
        pFeatureMap = NULL;
    }

    if (pFeatureHMap != NULL)
    {
        kh_destroy(fim, pFeatureHMap);
        pFeatureHMap = NULL;
    }

    if (hmap != NULL)
    {
        kh_destroy(fw, hmap);
        hmap = NULL;
    }
}
