//
// C Routines for Motif Kernel
//
// Source : Motif.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014-2016   J o h a n n e s  P a l m e
//


#include "Kebabs.h"
#include "KernelUtils.h"
#include "khash.h"
#include "ksort.h"

#define TF_NONE               0
#define TF_LEAF               1
#define TF_PATTERN_BLOCK      2
#define TF_USED_MOTIF         4
#define TF_MOTIF_IN_SAMPLE    8
#define TF_ANNOTATED_LEAF     16

#define TF_CHAR_MATCH         0
#define TF_GROUP_MATCH        1
#define TF_NEG_GROUP_MATCH    2

#define TF_RESET_FLAG         0
#define TF_RESET_ANNOT_ROOT   1

#define MAX_MOTIF_LENGTH      1000
#define INIT_POOL_SIZE        67108864    // 64MB
#define USER_INTERRUPT_LIMIT  100000
#define MAX_FEAT_VEC_LENGTH   1073741823

using namespace Rcpp;

#if __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

struct fmData
{
    uint32_t counter;
    uint32_t featureIndex;
    uint32_t sampleIndex;
    double    weight;
};

// hash map for mapping annotated feature to column index
KHASH_MAP_INIT_STR(fim, struct fmData);
KSORT_INIT(str, ksstr_t, ks_lt_str)

struct treeNodeMotif {
    struct indexBlock ib;
    uint32_t          value;
    uint8_t           flags;
};

struct prefTreeMotif {
    struct treeNodeMotif node[MAX_BLOCK];
};

struct intfFindMotifs {
    const char *seqptr;
    const char *annptr;
    uint32_t seqnchar;
    struct alphaInfo *alphaInf;
    struct prefTreeMotif *pTree;
    struct indexBlock *nullBlock;
    struct allIndMaps *allIndexMaps;
    int maxNoOfNodes;
    int *freeNode;
    int maxMotifLength;
    int offset;
    int rowIndex;
    int numUsedMotifs;
    uint64_t numNonzeroFeatures;
    int featuresPerSample;
    int svmIndex;
    char *keyPool;
    uint64_t *keyPoolSize;
    uint64_t *poolNextFree;
    uint64_t currFeatVecLength;
    int32_t *pFeatVecValue;
    uint32_t *pFeatVecIndex;
    uint32_t *unweightedPos;
    uint32_t fDim;
    uint32_t elemIndex;
    double kernelValue;
    bool getKernelValue;
    bool presence;
    bool zeroFeatures;
    bool markUsedOnly;
    bool markMotifsInSample;
    bool *printWarning;
    khash_t(fim) *featMap;
    IntegerVector *motifLengths;
    IntegerVector *annotationIndexMap;
    IntegerVector *unweightedPosStart;
    NumericMatrix *pErd;
    NumericMatrix *pProf;
    NumericMatrix *featureWeights;
};

struct intfStorePattern {
    char *pattern;
    int patternLength;
    int weightIndex;
    char *keyPool;
    double weight;
    uint64_t *keyPoolSize;
    uint64_t *poolNextFree;
    bool zeroFeatures;
    bool annSpec;
    struct prefTreeMotif *pTree;
    struct alphaInfo *alphaInf;
    khash_t(fim) *featMap;
    IntegerVector *motifLengths;
    IntegerVector *annotationIndexMap;
    // output
    int leafBlock;
    int motifLength;
    int motifPatternLength;
};

// for C heap release in on.exit hook
static khash_t(fim) *pFeatureHMap;
static uint32_t *pUnweightedPos;
static char *pKeypool;
static char **pKeys;

#if __clang__
#pragma clang diagnostic pop
#endif

void descendMotif(struct prefTreeMotif *pTree, int index, char kmer[], int level, int k,
                  struct alphaInfo *alphaInf)
{
    char kmern[MAX_MOTIF_LENGTH];
    uint32_t currMask, saveLevel;

    if (pTree->node[index].flags & TF_LEAF)
    {
        //visit node
        Rprintf("%s\n", kmer);
    }

    strcpy(kmern, kmer);

    saveLevel = level;

    // descend in alphabetical order
    for (int i = 0; i < (MAX_ALPHA_SIZE - 2); i++)
    {
        int next = (int) pTree->node[index].ib.idx[i];

        if (next > 0)
        {
            // append next char and pass copy
            strcpy(kmern, kmer);

            if (i == alphaInf->numAlphabetChars + 1)
            {
                kmern[level++] = '.';
                kmern[level] = '\0';

                descendMotif(pTree, next, kmern, level, k, alphaInf);

                level = saveLevel;
                kmern[level] = '\0';

            }
            else if (i != alphaInf->numAlphabetChars)
            {
                kmern[level++] = alphaInf->reverseIndexMap[i];
                kmern[level] = '\0';

                descendMotif(pTree, next, kmern, level, k, alphaInf);

                level = saveLevel;
                kmern[level] = '\0';
            }
            else
            {
                int saveLevelPattern;
                bool invertedGroup;
                saveLevelPattern = saveLevel;

                // loop through substitution group blocks
                int patternBlock = (int) pTree->node[index].ib.idx[alphaInf->numAlphabetChars];

                while (patternBlock != 0)
                {
                    for (uint32_t j = 0 ; j < pTree->node[patternBlock].value; j++)
                    {
                        invertedGroup = pTree->node[patternBlock].ib.idx[2 * j + 1] & (1UL << 31);
                        kmern[level++] = '[';

                        if (invertedGroup)
                            kmern[level++] = '^';

                        for (int g = 0; g < alphaInf->numAlphabetChars; g++)
                        {
                            currMask = 1UL << g;

                            if (invertedGroup)
                            {
                                if ((~pTree->node[patternBlock].ib.idx[2 * j + 1] & currMask) > 0)
                                    kmern[level++] = alphaInf->reverseIndexMap[g];
                            }
                            else
                            {
                                if ((pTree->node[patternBlock].ib.idx[2 * j + 1] & currMask) > 0)
                                    kmern[level++] = alphaInf->reverseIndexMap[g];
                            }
                        }

                        kmern[level++] = ']';
                        kmern[level] = '\0';

                        descendMotif(pTree, pTree->node[patternBlock].ib.idx[2 * j], kmern, level, k, alphaInf);

                        level = saveLevelPattern;
                        kmern[level] = '\0';
                    }

                    patternBlock = (int) pTree->node[patternBlock].ib.idx[MAX_ALPHA_SIZE - 2];
                }

                level = saveLevel;
                kmern[level] = '\0';
            }
        }
    }
}

void listTreeMotif(struct prefTreeMotif *pTree, int k, struct alphaInfo *alphaInf)
{
    char kmer[k];

    kmer[k] = '\0';
    descendMotif(pTree, 0, kmer, 0, k, alphaInf);
}


// $$$ TODO ensure that addressing via feature index is correct
// once the feature index is known it is used for addressing
// make sure features are in the correct order - either original and never sorted or always sorted
// or user flag that specifies whether they should be sorted

bool createMotifTree(ByteStringVector motifs, int maxMotifLength, struct prefTreeMotif *pTree,
                     int maxNoOfNodes, int *freeNode, struct indexBlock *nullBlock, bool *printWarning,
                     struct alphaInfo *alphaInf, bool setFeatureIndex)
{
    int i, j, motifChar, curr, index;
    uint32_t g, groupBits;
    bool substitutionGroup, invertedGroup, groupFound;

    // init first node
    pTree->node[0].ib = *nullBlock;
    pTree->node[0].flags = TF_NONE;
    *freeNode = 1;
    groupFound = FALSE;
    invertedGroup = FALSE;

    groupBits = 0;

    // build prefix tree for motifs
    for (i = 0; i < motifs.length; i++)
    {
        curr = 0;
        substitutionGroup = FALSE;

        for (j = 0; j < motifs.nchar[i]; j++)
        {
            motifChar = motifs.ptr[i][j];

            switch (motifChar)
            {
                case '[':
                {
                    if (substitutionGroup)
                    {
                        Rprintf("Substitution group within substitution group\n"
                                "        not allowed\n");
                        return(FALSE);
                    }

                    substitutionGroup = TRUE;
                    invertedGroup =FALSE;
                    groupFound =FALSE;
                    groupBits = 0;

                    if (pTree->node[curr].ib.idx[alphaInf->numAlphabetChars] != 0)
                    {
                        // set curr to first pattern block
                        curr = pTree->node[curr].ib.idx[alphaInf->numAlphabetChars];
                    }
                    else
                    {
                        //link first pattern block
                        pTree->node[curr].ib.idx[alphaInf->numAlphabetChars] = *freeNode;
                        curr = *freeNode;

                        if (curr>= maxNoOfNodes)
                        {
                            if (*printWarning)
                            {
                                Rprintf("Maximum number of nodes exceeded\n");
                                *printWarning = FALSE;
                            }
                            return(FALSE);
                        }

                        pTree->node[curr].flags = TF_PATTERN_BLOCK;
                        pTree->node[curr].value = 0; // init pattern count
                        pTree->node[curr].ib = *nullBlock;
                        *freeNode = *freeNode + 1;
                    }

                    break;
                }

                case '^':
                {
                    if (!substitutionGroup)
                    {
                        Rprintf("'^' only allowed at start of substitution group\n");
                        return(FALSE);
                    }

                    invertedGroup = TRUE;

                    break;
                }

                case ']':
                {
                    substitutionGroup = FALSE;

                    if (invertedGroup)
                        groupBits = ~groupBits;

                    // check if pattern exists already
                    while (pTree->node[curr].flags & TF_PATTERN_BLOCK)
                    {
                        for (g = 0; g < (2 * pTree->node[curr].value); g = g + 2)
                        {
                            if (pTree->node[curr].ib.idx[g + 1] == groupBits)
                            {
                                groupFound = TRUE;
                                curr = pTree->node[curr].ib.idx[g];
                                break;
                            }
                        }

                        if (groupFound || pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2] == 0)
                            break;
                        else
                            curr = pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2];
                    }

                    if (!groupFound)
                    {
                        // check if still space in pattern block for new substitution group
                        if (pTree->node[curr].value >= (MAX_ALPHA_SIZE / 2 - 1))
                        {
                            //link new parallel pattern block
                            pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2] = *freeNode;
                            curr = *freeNode;

                            if (curr>= maxNoOfNodes)
                            {
                                if (*printWarning)
                                {
                                    Rprintf("Maximum number of nodes exceeded\n");
                                    *printWarning = FALSE;
                                }
                                return(FALSE);
                            }

                            pTree->node[curr].flags = TF_PATTERN_BLOCK;
                            pTree->node[curr].value = 0; // init pattern count
                            pTree->node[curr].ib = *nullBlock;
                            *freeNode = *freeNode + 1;
                        }

                        // add substitution group to pattern block
                        pTree->node[curr].ib.idx[2 * pTree->node[curr].value] = *freeNode;
                        pTree->node[curr].ib.idx[2 * pTree->node[curr].value + 1] = groupBits;
                        pTree->node[curr].value++;
                        curr = *freeNode;

                        if (curr>= maxNoOfNodes)
                        {
                            if (*printWarning)
                            {
                                Rprintf("Maximum number of nodes exceeded\n");
                                *printWarning = FALSE;
                            }
                            return(FALSE);
                        }

                        pTree->node[curr].flags = TF_NONE;
                        pTree->node[curr].ib = *nullBlock;
                        *freeNode = *freeNode + 1;
                    }

                    if (j == motifs.nchar[i] - 1)
                    {
                        pTree->node[curr].flags |= TF_LEAF;
                        pTree->node[curr].value = 0;
                        pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2] = i;   // save motif index

                        if (setFeatureIndex)
                            pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1] = i;
                    }

                    break;
                }

                case '.':
                {
                    if (substitutionGroup)
                    {
                        Rprintf("Wildcard '.' not allowed in substitution group\n");
                        return(FALSE);
                    }

                    if (pTree->node[curr].ib.idx[alphaInf->numAlphabetChars + 1] > 0)
                        curr = pTree->node[curr].ib.idx[alphaInf->numAlphabetChars + 1];
                    else
                    {
                        // set continuation index for wildcard
                        pTree->node[curr].ib.idx[alphaInf->numAlphabetChars + 1] = *freeNode;
                        curr = *freeNode;

                        if (curr >= maxNoOfNodes)
                        {
                            if (*printWarning)
                            {
                                Rprintf("Maximum number of nodes exceeded\n");
                                *printWarning = FALSE;
                            }
                            return(FALSE);
                        }

                        pTree->node[curr].flags = TF_NONE;
                        pTree->node[curr].ib = *nullBlock;
                        *freeNode = *freeNode + 1;
                    }

                    if (j == motifs.nchar[i] - 1)
                    {
                        pTree->node[curr].flags |= TF_LEAF;
                        pTree->node[curr].value = 0;
                        pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2] = i;   // save motif index

                        if (setFeatureIndex)
                            pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1] = i;
                    }

                    break;
                }

                default:
                {
                    index = alphaInf->indexMap[motifChar];

                    if (index > -1)
                    {
                        if (substitutionGroup)
                        {
                            // add char to group
                            groupBits |= (1UL << index);
                        }
                        else
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
                                    return(FALSE);
                                }

                                pTree->node[curr].flags = TF_NONE;
                                pTree->node[curr].ib = *nullBlock;
                                *freeNode = *freeNode + 1;
                            }

                            if (j == motifs.nchar[i] - 1)
                            {
                                pTree->node[curr].flags |= TF_LEAF;
                                pTree->node[curr].value = 0;
                                pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2] = i;   // save motif index

                                if (setFeatureIndex)
                                    pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1] = i;
                            }
                        }
                    }
                    else
                    {
                        // interrupt when unknown character found
                        curr = 0;
                        return(FALSE);
                    }

                    break;
                }
            }
        }
    }

//    listTreeMotif(pTree, maxMotifLength, alphaInf);

    return(TRUE);
}

bool getLeaf(struct intfStorePattern *intf)
{
    int i, g, index, motifLength, currLength;
    uint32_t currBlock, group;
    uint8_t leafFlag;
    char motifChar;
    bool groupFound, substitutionGroup, invertedGroup;

    if (!intf->zeroFeatures)
        leafFlag = TF_USED_MOTIF;
    else
        leafFlag = TF_LEAF;

    currBlock = 0;
    motifLength = 0;
    substitutionGroup = FALSE;
    invertedGroup = FALSE;
    group = 0;

    for (i = 0; i < intf->patternLength; i++)
    {
        motifChar = intf->pattern[i];

        switch (motifChar)
        {
            case '.':
            {
                if (substitutionGroup)
                    return(FALSE);

                if (intf->pTree->node[currBlock].ib.idx[intf->alphaInf->numAlphabetChars + 1] == 0)
                {
                    intf->leafBlock = 0;
                    return(TRUE);
                }

                currBlock = intf->pTree->node[currBlock].ib.idx[intf->alphaInf->numAlphabetChars + 1];
                motifLength++;

                if (intf->pTree->node[currBlock].flags & leafFlag)
                    break;
                else
                    continue;
            }

            case '^':
            {
                if (!substitutionGroup)
                    return(FALSE);

                invertedGroup = TRUE;
                continue;
            }

            case '[':
            {
                if (intf->pTree->node[currBlock].ib.idx[intf->alphaInf->numAlphabetChars] == 0)
                {
                    intf->leafBlock = 0;
                    return(TRUE);
                }

                substitutionGroup = TRUE;
                invertedGroup = FALSE;
                group = 0;
                continue;
            }

            case ']':
            {
                substitutionGroup = FALSE;

                if (invertedGroup)
                    group = ~group;

                currBlock = intf->pTree->node[currBlock].ib.idx[intf->alphaInf->numAlphabetChars];
                groupFound = FALSE;

                while (intf->pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
                {
                    for (g = 0; g < (int) (2 * intf->pTree->node[currBlock].value); g = g + 2)
                    {
                        if (intf->pTree->node[currBlock].ib.idx[g + 1] == group)
                        {
                            groupFound = TRUE;
                            currBlock = intf->pTree->node[currBlock].ib.idx[g];
                            break;
                        }
                    }

                    if (groupFound)
                        break;

                    if (intf->pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] == 0)
                    {
                        Rprintf("Continuation block for motif not found\n");
                        return(FALSE);
                    }
                    else
                        currBlock = intf->pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2];
                }

                motifLength++;

                if (intf->pTree->node[currBlock].flags & leafFlag)
                    break;
                else
                    continue;
            }

            default:
            {
                index = intf->alphaInf->indexMap[(int) motifChar];

                if (index > -1)
                {
                    if (substitutionGroup)
                    {
                        group |= (1UL << index);
                        continue;
                    }

                    if (intf->pTree->node[currBlock].ib.idx[index] == 0)
                    {
                        intf->leafBlock = 0;
                        return(TRUE);
                    }

                    currBlock = intf->pTree->node[currBlock].ib.idx[index];

                    motifLength++;

                    if (intf->pTree->node[currBlock].flags & leafFlag)
                        break;
                    else
                        continue;
                }
                else
                {
                    Rprintf("Error in finding pattern in motif tree\n");
                    return(FALSE);
                }
            }
        }

        currLength = i + 1;

        if (intf->annSpec)
            currLength += motifLength;

        // verify if it is correct leaf (= deepest leaf) for current pattern
        if (currLength != intf->patternLength)
            continue;

        intf->leafBlock = currBlock;
        intf->motifLength = motifLength;
        intf->motifPatternLength = i + 1;

        return(TRUE);
    }

    intf->leafBlock = 0;
    return(TRUE);
}

void setFeatureIndex(struct prefTreeMotif *pTree, int maxMotifLength, int maxPatternLength, struct alphaInfo *alphaInf,
                     bool assignFeatureName, bool zeroFeatures, ByteStringVector motifs, IntegerVector *motifLengths,
                     bool annSpec, SEXP colnames, khash_t(fim) *featMap, char **keyPool, uint64_t *keyPoolSize,
                     uint64_t *poolNextFree, bool motifWithAnnotation, bool setUsedMotif)
{
    int i, j, g, curr, pos, numEntries, featureIndex, result;
    uint32_t currIndex, motifIndex;
    char *pKey, localKey[1024], digit, currMotif[2 * maxPatternLength + 1]; // including annotation chars
    struct intfStorePattern intf;
    struct fmData entryData;
    khiter_t iter;

    featureIndex = 0;

    if (!annSpec || motifWithAnnotation)
    {
        intf.pTree = pTree;
        intf.alphaInf = alphaInf;
        intf.annSpec = annSpec;
        intf.zeroFeatures = zeroFeatures;

        for (i = 0; i < motifs.length; i++)
        {
            intf.pattern = motifs.ptr[i];
            intf.patternLength = motifs.nchar[i];

            if (!getLeaf(&intf))
                return;

            // motif not found
            if (intf.leafBlock == 0)
                continue;

            if (!motifWithAnnotation)
            {
                pTree->node[intf.leafBlock].ib.idx[MAX_ALPHA_SIZE - 1] = featureIndex++;
                
                if (setUsedMotif)
                    pTree->node[intf.leafBlock].flags |= TF_USED_MOTIF;

                if (assignFeatureName)
                {
                    currMotif[motifs.nchar[i]] = '\0';

                    for (g = 0; g < motifs.nchar[i]; g++)
                        currMotif[g] = motifs.ptr[i][g];

                    SET_STRING_ELT(colnames, featureIndex - 1, Rf_mkChar(currMotif));
                }
            }
            else
            {
                motifIndex = pTree->node[intf.leafBlock].ib.idx[MAX_ALPHA_SIZE - 2];
                pos = 0;

                for (j = 0; j < 10; j++)
                {
                    localKey[9 - pos++] = '0' + motifIndex % 10;
                    motifIndex = motifIndex / 10;

                }

                localKey[pos++] = '_';
                memcpy(&(localKey[pos]), &(intf.pattern[intf.motifPatternLength]), intf.motifLength);
                pos += intf.motifLength;
                localKey[pos++] = '\0';

                // new key
                if (*keyPoolSize < *poolNextFree + 10 + intf.motifLength)
                {
                    // realloc pool
                    (*keyPoolSize) *= 2;
                    pKeypool = (char *) R_Realloc(*keyPool, *keyPoolSize, char);
                    *keyPool = pKeypool;
                }

                pKey = &((*keyPool)[*poolNextFree]);
                memcpy(&((*keyPool)[*poolNextFree]), localKey, pos);
                (*poolNextFree) += pos;

                iter = kh_put(fim, featMap, pKey, &result);

                if (result != -1)
                {
                    // create annotation entry
                    entryData.sampleIndex = MAXUINT32;
                    entryData.featureIndex = i;
                    entryData.counter = 0;

                    kh_value(featMap, iter) = entryData;
                }
                else
                {
                    // $$$  TODO remove Rprintf
                    Rprintf("Annotated motif could not be stored in hash map\n");
                    return;
                }
            }
        }
    }
    else
    {
        numEntries = kh_size(featMap);
        curr = 0;

        entryData.sampleIndex = MAXUINT32;
        entryData.counter = 0;

        if (numEntries > 0)
        {
            // get key pointers
            char **keys = (char **) R_Calloc(numEntries, char *);
            pKeys = keys;
            curr = 0;

            for (iter = kh_begin(featMap); iter != kh_end(featMap); iter++)
            {
                if (kh_exist(featMap, iter))
                    keys[curr++] = (char *)kh_key(featMap, iter);
            }

            // sort keys
            ks_mergesort(str, numEntries, (const char **) keys, NULL);

            // loop through keys in sorted order and set increasing col index
            for (i = 0; i < numEntries; i++)
            {
                iter = kh_get(fim, featMap, keys[i]);

                if (iter != kh_end(featMap))
                {
                    if (assignFeatureName)
                    {
                        currIndex = 0;
                        curr = 0;

                        digit = kh_key(featMap, iter)[curr++];

                        while (digit != '_')
                        {
                            currIndex = currIndex*10 + digit - '0';
                            digit = kh_key(featMap, iter)[curr++];
                        }

                        currMotif[motifs.nchar[currIndex] + (*motifLengths)[currIndex]] = '\0';

                        for (j = 0; j < motifs.nchar[currIndex]; j++)
                            currMotif[j] = motifs.ptr[currIndex][j];

                        memcpy(&(currMotif[motifs.nchar[currIndex]]), &(keys[i][curr]), (*motifLengths)[currIndex]);
                        SET_STRING_ELT(colnames, i, Rf_mkChar(currMotif));
                    }

                    entryData.featureIndex = featureIndex++;
                    kh_value(featMap, iter) = entryData;
                }
            }

            R_Free(pKeys);
            pKeys = NULL;
        }
    }

    return;
}

void resetInfoInTree(uint8_t task, uint8_t flagToReset, struct prefTreeMotif *pTree,
                     int maxMotifLength, struct alphaInfo *alphaInf, bool zeroFeatures)
{
    uint32_t currBlock, blockStack[4 * (maxMotifLength + 1)];
    int currIndex, maxBlockIndex, currStack;

    maxBlockIndex = MAX_ALPHA_SIZE - 3; // highest valid index for motif
    currBlock = 0;
    currIndex = 0;
    currStack = -1;

    while (currStack >= 0 || currIndex <= maxBlockIndex)
    {
        if (currIndex == 0 && (pTree->node[currBlock].flags & TF_LEAF))
        {
            if (task == TF_RESET_FLAG)
                pTree->node[currBlock].flags &= ~flagToReset;
            else if (task == TF_RESET_ANNOT_ROOT)
                pTree->node[currBlock].value = 0;
            else
                return;
        }

        if (pTree->node[currBlock].ib.idx[currIndex] != 0)
        {
            if (pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
            {
                if ((currIndex + 2 <= (int) (2 * pTree->node[currBlock].value - 2)) ||
                    ((currIndex + 2 == MAX_ALPHA_SIZE - 2) &&
                     (pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] != 0)))

                {
                    blockStack[++currStack] = currIndex + 2;
                    blockStack[++currStack] = currBlock;
                }
            }
            else
            {
                if (currIndex + 1 <= alphaInf->numAlphabetChars + 1)
                {
                    blockStack[++currStack] = currIndex + 1;
                    blockStack[++currStack] = currBlock;
                }
            }

            currBlock = pTree->node[currBlock].ib.idx[currIndex];
            currIndex = 0;
        }
        else
        {
            if (pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
            {
                currIndex = currIndex + 2;

                if (currIndex < (int) (2 * pTree->node[currBlock].value))
                    continue;

                if (pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] != 0)
                {
                    // current block finished afterwards - nothing put on stack
                    currBlock = pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2];
                    currIndex = 0;
                    continue;
                }
            }
            else
            {
                currIndex++;
                while ((pTree->node[currBlock].ib.idx[currIndex] == 0) &&
                       (currIndex <= alphaInf->numAlphabetChars + 1))
                    currIndex++;

                if (currIndex <= alphaInf->numAlphabetChars + 1)
                    continue;
            }

            if (currStack == -1)
                continue;

            currBlock = blockStack[currStack--];
            currIndex = blockStack[currStack--];
        }
    }

    return;
}

RcppExport void findUnweightedPositions(ByteStringVector motifs, IntegerVector *unweightedPosStart,
                                        uint32_t **unweightedPos)
{
    int i, j, arraySize, currIndex, currPos;
    bool substGroupOpen;

    arraySize = motifs.length * 2; // start with two unweighted positions per motif on average
    pUnweightedPos = (uint32_t *) R_Calloc(arraySize, uint32_t);
    *unweightedPos = pUnweightedPos;
    currIndex = 0;

    for (i = 0; i < motifs.length; i++)
    {
        (*unweightedPosStart)[i] = currIndex;
        currPos = -1;
        substGroupOpen = FALSE;

        for (j = 0; j < motifs.nchar[i]; j++)
        {
            if (substGroupOpen)
            {
                if (motifs.ptr[i][j] == ']')
                {
                    currPos++;
                    substGroupOpen = FALSE;
                }
            }
            else
            {
                if (motifs.ptr[i][j] == '[')
                    substGroupOpen = TRUE;
                else
                {
                    currPos++;

                    if (motifs.ptr[i][j] == '.')
                    {
                        if (currIndex >=  arraySize)
                        {
                            arraySize = arraySize * 2;
                            pUnweightedPos = R_Realloc(pUnweightedPos, arraySize, uint32_t);
                            *unweightedPos = pUnweightedPos;
                        }

                        (*unweightedPos)[currIndex++] = currPos;
                    }
                }
            }
        }
    }

    (*unweightedPosStart)[motifs.length] = currIndex;

    return;
}

bool findAnnotatedMotif(uint32_t leafBlock, uint32_t motifStart, uint32_t currPos,
                        struct intfFindMotifs *intf)
{
    char *pKey, localKey[1024];
    int i, j, result, motifLength, motifIndex, motifIndexTemp, pos, temp, parts;
    uint32_t startIndex;
    double weight;
    struct fmData entryData;
    khiter_t iter;

    motifIndex = intf->pTree->node[leafBlock].ib.idx[MAX_ALPHA_SIZE - 2];
    motifIndexTemp = motifIndex;
    motifLength = (*intf->motifLengths)[motifIndex];
    pos = 0;

    for (i = 0; i < 10; i++)
    {
        localKey[9 - pos++] = '0' + motifIndexTemp % 10;
        motifIndexTemp = motifIndexTemp / 10;

    }

    localKey[pos++] = '_';
    memcpy(&(localKey[pos]), &(intf->annptr[motifStart]), motifLength);
    pos += motifLength;
    localKey[pos++] = '\0';

    // wildcard positions for annotation
    for (i = (*intf->unweightedPosStart)[motifIndex]; i < (*intf->unweightedPosStart)[motifIndex + 1]; i++)
        localKey[pos - motifLength + intf->unweightedPos[i] - 1] = '.';

    iter = kh_get(fim, intf->featMap, localKey);

    if (iter == kh_end(intf->featMap))
    {
        if (intf->pProf != NULL)
            return(TRUE);

        // new key
        if (*intf->keyPoolSize < *intf->poolNextFree + 10 + motifLength)
        {
            // realloc pool
            *intf->keyPoolSize *= 2;
            pKeypool = (char *) R_Realloc(intf->keyPool, *intf->keyPoolSize, char);
            intf->keyPool = pKeypool;
        }

        pKey = &(intf->keyPool[*intf->poolNextFree]);
        memcpy(&(intf->keyPool[*intf->poolNextFree]), localKey, pos);
        (*intf->poolNextFree) += pos;

        iter = kh_put(fim, intf->featMap, pKey, &result);

        if (result != -1)
        {
            entryData.counter = 0;
            entryData.featureIndex = intf->rowIndex; // temporary storage of sample in col index!
            entryData.sampleIndex = MAXUINT32;

            kh_value(intf->featMap, iter) = entryData;
            intf->featuresPerSample++;
            intf->numUsedMotifs += 1;
            intf->numNonzeroFeatures++;
        }
        else
        {
            Rprintf("Annotated motif could not be stored in hash map\n");
            return(FALSE);
        }
    }
    else
    {
        entryData = kh_value(intf->featMap, iter);

        if (!intf->markUsedOnly)
        {
            if (intf->pErd != NULL)
            {
                if (intf->presence)
                {
                    if ((*intf->pErd)(intf->rowIndex, entryData.featureIndex) == 0)
                    {
                        (*intf->pErd)(intf->rowIndex, entryData.featureIndex) = 1;

                        if (intf->getKernelValue)
                            intf->kernelValue++;
                    }
                }
                else
                {
                    temp = (*intf->pErd)(intf->rowIndex, entryData.featureIndex);
                    (*intf->pErd)(intf->rowIndex, entryData.featureIndex) = temp + 1;

                    if (intf->getKernelValue)
                        intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                }
            }
            else
            {
                if (intf->pProf != NULL)
                {
                    // $$$ TODO change weight lookup for position specific and quadratic kernel
                    weight = (*intf->featureWeights)(intf->svmIndex, entryData.featureIndex);

                    if ((*intf->unweightedPosStart)[motifIndex] < (*intf->unweightedPosStart)[motifIndex + 1])
                    {
                        startIndex = (*intf->unweightedPosStart)[motifIndex];
                        parts = motifLength - ((*intf->unweightedPosStart)[motifIndex + 1] - (*intf->unweightedPosStart)[motifIndex]);
                    }
                    else
                    {
                        startIndex = MAXUINT32;
                        parts = motifLength;
                    }

                    weight /= parts;

                    for (j = motifStart; j < (int)motifStart + motifLength; j++)
                    {
                        while ((startIndex < (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1]) &&
                               (j - motifStart) > intf->unweightedPos[startIndex])
                            startIndex++;

                        if (startIndex == (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1])
                            startIndex = MAXUINT32;

                        if (startIndex == MAXUINT32 || (intf->unweightedPos[startIndex] !=
                                                        (j - motifStart)))
                            (*intf->pProf)(intf->rowIndex, j) += weight;
                    }
                }
                else
                {
                    if (entryData.sampleIndex != (uint32_t) intf->rowIndex)
                    {
                        entryData.sampleIndex = intf->rowIndex;
                        entryData.counter = 1;
                        kh_value(intf->featMap, iter) = entryData;

                        if (intf->getKernelValue)
                            intf->kernelValue++;
                    }
                    else
                    {
                        if (!intf->presence)
                        {
                            temp = entryData.counter;
                            entryData.counter += 1;
                            kh_value(intf->featMap, iter) = entryData;

                            if (intf->getKernelValue)
                                intf->kernelValue += -temp * temp + (temp + 1) * (temp + 1);
                        }
                    }
                }
            }
        }
        else
        {
            if (entryData.featureIndex != (uint32_t) intf->rowIndex)
            {
                entryData.featureIndex = intf->rowIndex; // temporary storage of sample in col index!
                kh_value(intf->featMap, iter) = entryData;
                intf->featuresPerSample++;
                intf->numNonzeroFeatures++;
            }
        }
    }

    return(TRUE);
}

bool descendOnBranch(uint32_t startPos, uint32_t endPos, uint32_t startBlock, uint32_t motifStart,
                     struct intfFindMotifs *intf)
{
    uint32_t i, j, g, l, n, upper, curr, currNext, currGroup, motifIndex, startIndex;
    int index, temp, parts;
    double weight;

    for (i = startPos; i < endPos; i++)
    {
        if (startBlock == 0)
            motifStart = i;

        curr = startBlock;
        l = i;

        if ((l + intf->maxMotifLength) > endPos)
            upper = endPos - l;
        else
            upper = intf->maxMotifLength;

        for (j=0; j < upper; j++)
        {
            index = intf->alphaInf->seqIndexMap[(int)intf->seqptr[l++]];

            if (index > -1)
            {
                if (intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars + 1] != 0)
                {
                    // descend on wildcard
                    currNext = intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars + 1];

                    if (intf->pTree->node[currNext].flags & TF_LEAF)
                    {
                        if (intf->annptr != NULL)
                        {
                            if (!findAnnotatedMotif(currNext, motifStart, l - 1, intf))
                                return(FALSE);
                        }
                        else
                        {
                            if (intf->markUsedOnly)
                            {
                                if (!(intf->pTree->node[currNext].flags & TF_USED_MOTIF))
                                {
                                    intf->pTree->node[currNext].flags |= TF_USED_MOTIF;
                                    intf->numUsedMotifs++;
                                }

                                if (intf->markMotifsInSample)
                                {
                                    if (!(intf->pTree->node[currNext].flags & TF_MOTIF_IN_SAMPLE))
                                    {
                                        intf->pTree->node[currNext].flags |= TF_MOTIF_IN_SAMPLE;
                                        intf->numNonzeroFeatures++;
                                    }
                                }
                            }
                            else
                            {
                                if (!intf->presence)
                                {
                                    if (intf->pErd != NULL)
                                    {
                                        temp = (*intf->pErd)(intf->rowIndex, intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 1]);
                                        (*intf->pErd)(intf->rowIndex, intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 1]) = temp + 1;
                                        intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                                    }
                                    else
                                    {
                                        if (intf->pProf != NULL &&
                                            (intf->pTree->node[currNext].flags & TF_USED_MOTIF))
                                        {
                                            motifIndex = intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 2];
                                            weight = (*intf->featureWeights)(intf->svmIndex, intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 1]);

                                            if ((*intf->unweightedPosStart)[motifIndex] < (*intf->unweightedPosStart)[motifIndex + 1])
                                            {
                                                startIndex = (*intf->unweightedPosStart)[motifIndex];
                                                parts = l - motifStart - ((*intf->unweightedPosStart)[motifIndex + 1] - (*intf->unweightedPosStart)[motifIndex]);
                                            }
                                            else
                                            {
                                                startIndex = MAXUINT32;
                                                parts = l - motifStart;
                                            }

                                            weight /= parts;

                                            for (n = motifStart; n < l; n++)
                                            {
                                                while ((startIndex < (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1]) &&
                                                       (n - motifStart) > intf->unweightedPos[startIndex])
                                                    startIndex++;

                                                if (startIndex == (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1])
                                                    startIndex = MAXUINT32;

                                                if (startIndex == MAXUINT32 || (intf->unweightedPos[startIndex] != (n - motifStart)))
                                                    (*intf->pProf)(intf->rowIndex, n) += weight;
                                            }
                                        }
                                        else
                                        {
                                            temp = intf->pTree->node[currNext].value;
                                            intf->pTree->node[currNext].value = temp + 1;
                                            intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                                        }
                                    }
                                }
                                else
                                {
                                    if (intf->pErd != NULL)
                                        (*intf->pErd)(intf->rowIndex, intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 1]) = 1;
                                    else
                                    {
                                        if (intf->pTree->node[currNext].value == 0)
                                        {
                                            intf->pTree->node[currNext].value = 1;
                                            intf->kernelValue++;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (!descendOnBranch(l, i + upper, currNext, motifStart, intf))
                        return(FALSE);
                }

                if (intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars] != 0)
                {
                    // descend on substitution group
                    currNext = intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars];

                    while (currNext != 0)
                    {
                        for (g = 0; g < 2 * intf->pTree->node[currNext].value; g = g + 2)
                        {
                            if ((intf->pTree->node[currNext].ib.idx[g + 1] & (1UL << index)) > 0)
                            {
                                currGroup = intf->pTree->node[currNext].ib.idx[g];

                                if (intf->pTree->node[currGroup].flags & TF_LEAF)
                                {
                                    if (intf->annptr != NULL)
                                    {
                                        if (!findAnnotatedMotif(currGroup, motifStart, l - 1, intf))
                                            return(FALSE);
                                    }
                                    else
                                    {
                                        if (intf->markUsedOnly)
                                        {
                                            if (!(intf->pTree->node[currGroup].flags & TF_USED_MOTIF))
                                            {
                                                intf->pTree->node[currGroup].flags |= TF_USED_MOTIF;
                                                intf->numUsedMotifs++;
                                            }

                                            if (intf->markMotifsInSample)
                                            {
                                                if (!(intf->pTree->node[currGroup].flags & TF_MOTIF_IN_SAMPLE))
                                                {
                                                    intf->pTree->node[currGroup].flags |= TF_MOTIF_IN_SAMPLE;
                                                    intf->numNonzeroFeatures++;
                                                }
                                            }
                                        }
                                        else
                                        {
                                            if (!intf->presence)
                                            {
                                                if (intf->pErd != NULL)
                                                {
                                                    temp = (*intf->pErd)(intf->rowIndex, intf->pTree->node[currGroup].ib.idx[MAX_ALPHA_SIZE - 1]);
                                                    (*intf->pErd)(intf->rowIndex, intf->pTree->node[currGroup].ib.idx[MAX_ALPHA_SIZE - 1]) = temp + 1;
                                                    intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                                                }
                                                else
                                                {
                                                    if (intf->pProf != NULL &&
                                                        (intf->pTree->node[currGroup].flags & TF_USED_MOTIF))
                                                    {
                                                        motifIndex = intf->pTree->node[currGroup].ib.idx[MAX_ALPHA_SIZE - 2];
                                                        weight = (*intf->featureWeights)(intf->svmIndex, intf->pTree->node[currGroup].ib.idx[MAX_ALPHA_SIZE - 1]);

                                                        if ((*intf->unweightedPosStart)[motifIndex] < (*intf->unweightedPosStart)[motifIndex + 1])
                                                        {
                                                            startIndex = (*intf->unweightedPosStart)[motifIndex];
                                                            parts = l - motifStart - ((*intf->unweightedPosStart)[motifIndex + 1] - (*intf->unweightedPosStart)[motifIndex]);
                                                        }
                                                        else
                                                        {
                                                            startIndex = MAXUINT32;
                                                            parts = l - motifStart;
                                                        }

                                                        weight /= parts;

                                                        for (n = motifStart; n < l; n++)
                                                        {
                                                            while ((startIndex < (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1]) &&
                                                                   (n - motifStart) > intf->unweightedPos[startIndex])
                                                                startIndex++;

                                                            if (startIndex == (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1])
                                                                startIndex = MAXUINT32;

                                                            if (startIndex == MAXUINT32 || (intf->unweightedPos[startIndex] != (n - motifStart)))
                                                                (*intf->pProf)(intf->rowIndex, n) += weight;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        temp = intf->pTree->node[currGroup].value;
                                                        intf->pTree->node[currGroup].value = temp + 1;
                                                        intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (intf->pErd != NULL)
                                                    (*intf->pErd)(intf->rowIndex, intf->pTree->node[currGroup].ib.idx[MAX_ALPHA_SIZE - 1]) = 1;
                                                else
                                                {
                                                    if (intf->pTree->node[currGroup].value == 0)
                                                    {
                                                        intf->pTree->node[currGroup].value = 1;
                                                        intf->kernelValue++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                if (!descendOnBranch(l, i + upper, currGroup, motifStart, intf))
                                    return(FALSE);
                            }
                        }

                        currNext = intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 2];
                    }
                }

                if (intf->pTree->node[curr].ib.idx[index] == 0)
                {
                    if (startBlock != 0)
                        return(TRUE);
                    else
                        break;
                }

                // descend on exact character match
                curr = intf->pTree->node[curr].ib.idx[index];

                if (intf->pTree->node[curr].flags & TF_LEAF)
                {
                    if (intf->annptr != NULL)
                    {
                        if (!findAnnotatedMotif(curr, motifStart, l - 1, intf))
                            return(FALSE);
                    }
                    else
                    {
                        if (intf->markUsedOnly)
                        {
                            if (!(intf->pTree->node[curr].flags & TF_USED_MOTIF))
                            {
                                intf->pTree->node[curr].flags |= TF_USED_MOTIF;
                                intf->numUsedMotifs++;
                            }

                            if (intf->markMotifsInSample)
                            {
                                if (!(intf->pTree->node[curr].flags & TF_MOTIF_IN_SAMPLE))
                                {
                                    intf->pTree->node[curr].flags |= TF_MOTIF_IN_SAMPLE;
                                    intf->numNonzeroFeatures++;
                                }
                            }
                        }
                        else
                        {
                            if (!intf->presence)
                            {
                                if (intf->pErd != NULL)
                                {
                                    temp = (*intf->pErd)(intf->rowIndex, intf->pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1]);
                                    (*intf->pErd)(intf->rowIndex, intf->pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1]) = temp + 1;
                                    intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                                }
                                else
                                {
                                    if (intf->pProf != NULL &&
                                        (intf->pTree->node[curr].flags & TF_USED_MOTIF))
                                    {
                                        motifIndex = intf->pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 2];
                                        weight = (*intf->featureWeights)(intf->svmIndex, intf->pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1]);

                                        if ((*intf->unweightedPosStart)[motifIndex] < (*intf->unweightedPosStart)[motifIndex + 1])
                                        {
                                            startIndex = (*intf->unweightedPosStart)[motifIndex];
                                            parts = l - motifStart - ((*intf->unweightedPosStart)[motifIndex + 1] -
                                                                      (*intf->unweightedPosStart)[motifIndex]);
                                        }
                                        else
                                        {
                                            startIndex = MAXUINT32;
                                            parts = l - motifStart;
                                        }

                                        weight /= parts;

                                        for (n = motifStart; n < l; n++)
                                        {
                                            while ((startIndex < (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1]) &&
                                                   (n - motifStart) > intf->unweightedPos[startIndex])
                                                startIndex++;

                                            if (startIndex == (uint32_t) (*intf->unweightedPosStart)[motifIndex + 1])
                                                startIndex = MAXUINT32;

                                            if (startIndex == MAXUINT32 || (intf->unweightedPos[startIndex] != (n - motifStart)))
                                                (*intf->pProf)(intf->rowIndex, n) += weight;
                                        }
                                    }
                                    else
                                    {
                                        temp = intf->pTree->node[curr].value;
                                        intf->pTree->node[curr].value = temp + 1;
                                        intf->kernelValue = intf->kernelValue - (double) temp * temp + (double) (temp + 1) * (temp + 1);
                                    }
                                }
                            }
                            else
                            {
                                if (intf->pErd != NULL)
                                    (*intf->pErd)(intf->rowIndex, intf->pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1]) = 1;
                                else
                                {
                                    if (intf->pTree->node[curr].value == 0)
                                    {
                                        intf->pTree->node[curr].value = 1;
                                        intf->kernelValue++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
                break;
        }
    }

    return(TRUE);
}

void findMotifs(struct intfFindMotifs *intf)
{
    intf->kernelValue = 0;

    if (!descendOnBranch(0, intf->seqnchar, 0, 0, intf))
        intf->kernelValue = -1;

    return;
}

void descendOnBranchPos(uint32_t startPos, uint32_t endPos, uint32_t startBlock, uint32_t *motifBegin,
                        struct intfFindMotifs *intf)
{
    uint32_t i, j, g, l, upper, curr, currNext, currGroup;
    int index;

    for (i = startPos; i < endPos; i++)
    {
        if (startBlock == 0)
            *motifBegin = i;

        curr = startBlock;
        l = i;

        if ((l + intf->maxMotifLength) > intf->seqnchar)
            upper = intf->seqnchar - l;
        else
            upper = intf->maxMotifLength;

        for (j=0; j < upper; j++)
        {
            index = intf->alphaInf->seqIndexMap[(int)intf->seqptr[l++]];

            if (index > -1)
            {
                if (intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars + 1] != 0)
                {
                    currNext = intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars + 1];

                    if (intf->pTree->node[currNext].flags & TF_LEAF)
                    {
                        if (intf->markUsedOnly)
                        {
                            if (!(intf->pTree->node[currNext].flags & TF_USED_MOTIF))
                            {
                                intf->pTree->node[currNext].flags |= TF_USED_MOTIF;
                                intf->numUsedMotifs++;
                            }

                            if (intf->markMotifsInSample)
                            {
                                if (!(intf->pTree->node[currNext].flags & TF_MOTIF_IN_SAMPLE))
                                {
                                    intf->pTree->node[currNext].flags |= TF_MOTIF_IN_SAMPLE;
                                    intf->numNonzeroFeatures++;
                                }
                            }
                        }
                        else
                        {
                            if (intf->elemIndex >= intf->currFeatVecLength)
                            {
                                intf->currFeatVecLength = intf->currFeatVecLength * 1.4;
                                intf->pFeatVecIndex = (uint32_t *) R_Realloc(intf->pFeatVecIndex, intf->currFeatVecLength, uint32_t);
                                intf->pFeatVecValue = (int32_t *) R_Realloc(intf->pFeatVecValue, intf->currFeatVecLength, int32_t);
                            }

                            intf->pFeatVecValue[intf->elemIndex] = *motifBegin - intf->offset + 1;
                            intf->pFeatVecIndex[intf->elemIndex++] = intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 1];
                            intf->kernelValue++;
                        }
                    }

                    descendOnBranchPos(l, l + 1, currNext, motifBegin, intf);
                }

                if (intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars] != 0)
                {
                    currNext = intf->pTree->node[curr].ib.idx[intf->alphaInf->numAlphabetChars];

                    while (currNext != 0)
                    {
                        for (g = 0; g < 2 * intf->pTree->node[currNext].value; g = g + 2)
                        {
                            if ((intf->pTree->node[currNext].ib.idx[g + 1] & (1UL << index)) > 0)
                            {
                                currGroup = intf->pTree->node[currNext].ib.idx[g];

                                if (intf->pTree->node[currGroup].flags & TF_LEAF)
                                {
                                    if (intf->markUsedOnly)
                                    {
                                        if (!(intf->pTree->node[currGroup].flags & TF_USED_MOTIF))
                                        {
                                            intf->pTree->node[currGroup].flags |= TF_USED_MOTIF;
                                            intf->numUsedMotifs++;
                                        }

                                        if (intf->markMotifsInSample)
                                        {
                                            if (!(intf->pTree->node[currGroup].flags & TF_MOTIF_IN_SAMPLE))
                                            {
                                                intf->pTree->node[currGroup].flags |= TF_MOTIF_IN_SAMPLE;
                                                intf->numNonzeroFeatures++;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (intf->elemIndex >= intf->currFeatVecLength)
                                        {
                                            intf->currFeatVecLength = intf->currFeatVecLength * 1.4;
                                            intf->pFeatVecIndex = (uint32_t *) R_Realloc(intf->pFeatVecIndex, intf->currFeatVecLength,
                                                                                       uint32_t);
                                            intf->pFeatVecValue = (int32_t *) R_Realloc(intf->pFeatVecValue, intf->currFeatVecLength,
                                                                                      int32_t);
                                        }

                                        intf->pFeatVecValue[intf->elemIndex] = *motifBegin - intf->offset + 1;
                                        intf->pFeatVecIndex[intf->elemIndex++] = intf->pTree->node[currGroup].ib.idx[MAX_ALPHA_SIZE - 1];
                                        intf->kernelValue++;
                                    }
                                }

                                descendOnBranchPos(l, l + 1, currGroup, motifBegin, intf);
                            }
                        }

                        currNext = intf->pTree->node[currNext].ib.idx[MAX_ALPHA_SIZE - 2];
                    }
                }

                if (intf->pTree->node[curr].ib.idx[index] == 0)
                    break;

                // descend on character match
                curr = intf->pTree->node[curr].ib.idx[index];

                if (intf->pTree->node[curr].flags & TF_LEAF)
                {
                    if (intf->markUsedOnly)
                    {
                        if (!(intf->pTree->node[curr].flags & TF_USED_MOTIF))
                        {
                            intf->pTree->node[curr].flags |= TF_USED_MOTIF;
                            intf->numUsedMotifs++;
                        }

                        if (intf->markMotifsInSample)
                        {
                            if (!(intf->pTree->node[curr].flags & TF_MOTIF_IN_SAMPLE))
                            {
                                intf->pTree->node[curr].flags |= TF_MOTIF_IN_SAMPLE;
                                intf->numNonzeroFeatures++;
                            }
                        }
                    }
                    else
                    {
                        if (intf->elemIndex >= intf->currFeatVecLength)
                        {
                            intf->currFeatVecLength = intf->currFeatVecLength * 1.4;
                            intf->pFeatVecIndex = (uint32_t *) R_Realloc(intf->pFeatVecIndex, intf->currFeatVecLength, uint32_t);
                            intf->pFeatVecValue = (int32_t *) R_Realloc(intf->pFeatVecValue, intf->currFeatVecLength, int32_t);
                        }

                        intf->pFeatVecValue[intf->elemIndex] = *motifBegin - intf->offset + 1;
                        intf->pFeatVecIndex[intf->elemIndex++] = intf->pTree->node[curr].ib.idx[MAX_ALPHA_SIZE - 1];
                        intf->kernelValue++;
                    }
                }
            }
            else
                break;
        }
    }
}

void findMotifsForPos(struct intfFindMotifs *intf)
{
    uint32_t motifBegin;

    intf->kernelValue = 0;

    descendOnBranchPos(0, intf->seqnchar, 0, &motifBegin, intf);

    return;
}

void getNonzeroMotifs(bool annSpec, int32_t *featVectorValue, uint32_t *featVectorIndex,
                      struct prefTreeMotif *pTree, khash_t(fim) *featMap, uint32_t fDim,
                      uint32_t sampleIndex, struct alphaInfo *alphaInf, int maxMotifLength)
{
    uint32_t currBlock, elemIndex;
    uint32_t blockStack[4 * (2 * maxMotifLength + 1)];
    int currIndex, maxBlockIndex, currStack, numEntries;
    khiter_t iter;

    elemIndex = sampleIndex * fDim;

    if (!annSpec)
    {
        // walk tree to generate sparse feature vector
        maxBlockIndex = MAX_ALPHA_SIZE - 3; // highest valid index for motif
        currBlock = 0;
        currIndex = 0;
        currStack = -1;

        while (currStack >= 0 || currIndex <= maxBlockIndex)
        {
            if (currIndex == 0 && (pTree->node[currBlock].flags & TF_LEAF))
            {
                if (pTree->node[currBlock].value > 0)
                {
                    featVectorValue[elemIndex] = pTree->node[currBlock].value;
                    pTree->node[currBlock].value = 0;
                    featVectorIndex[elemIndex++] = pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 1];
                }
            }

            if (pTree->node[currBlock].ib.idx[currIndex] != 0)
            {
                if (pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
                {
                    if ((currIndex + 2 <= (int) (2 * pTree->node[currBlock].value - 2)) ||
                        ((currIndex + 2 == MAX_ALPHA_SIZE - 2) &&
                         (pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] != 0)))

                    {
                        blockStack[++currStack] = currIndex + 2;
                        blockStack[++currStack] = currBlock;
                    }
                }
                else
                {
                    if (currIndex + 1 <= alphaInf->numAlphabetChars + 1)
                    {
                        blockStack[++currStack] = currIndex + 1;
                        blockStack[++currStack] = currBlock;
                    }
                }

                currBlock = pTree->node[currBlock].ib.idx[currIndex];
                currIndex = 0;
            }
            else
            {
                if (pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
                {
                    currIndex = currIndex + 2;

                    if (currIndex < (int)(2 * pTree->node[currBlock].value))
                        continue;

                    if (pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] != 0)
                    {
                        // current block finished afterwards - nothing put on stack
                        currBlock = pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2];
                        currIndex = 0;
                        continue;
                    }
                }
                else
                {
                    currIndex++;

                    while ((pTree->node[currBlock].ib.idx[currIndex] == 0) &&
                           (currIndex <= alphaInf->numAlphabetChars + 1))
                        currIndex++;

                    if (currIndex <= alphaInf->numAlphabetChars + 1)
                        continue;
                }

                if (currStack == -1)
                    continue;

                currBlock = blockStack[currStack--];
                currIndex = blockStack[currStack--];
            }
        }

        featVectorIndex[elemIndex] = MAXUINT32;
        featVectorValue[elemIndex] = MAXINT32;
    }
    else
    {
        numEntries = kh_size(featMap);

        struct fmData entryData;

        if (numEntries > 0)
        {
            // sorting happens before creation of kernel matrix
            for (iter = kh_begin(featMap); iter != kh_end(featMap); iter++)
            {
                if (kh_exist(featMap, iter))
                {
                    entryData = kh_value(featMap, iter);

                    if (entryData.sampleIndex == sampleIndex && entryData.counter > 0)
                    {
                        featVectorValue[elemIndex] = entryData.counter;
                        featVectorIndex[elemIndex++] = entryData.featureIndex;
                    }
                }
            }

            featVectorIndex[elemIndex] = MAXUINT32;
            featVectorValue[elemIndex] = MAXINT32;
        }
    }

    return;
}

void getNonzeroMotifsERS(bool annSpec, struct prefTreeMotif *pTree, khash_t(fim) *featMap, struct alphaInfo *alphaInf,
                         int maxMotifLength, uint32_t sampleIndex, int numUsedFeatures, SEXP slot_j, SEXP slot_x, int *jIdx,
                         uint32_t *featVectorIndex, int32_t *featVectorValue, double normValue, bool normalized, bool zeroFeatures)
{
    uint32_t currBlock, blockStack[4 * (maxMotifLength + 1)];
    int i, currIndex, maxBlockIndex, currStack, numEntries, numFeatures, curr;
    khiter_t iter;


    if (!annSpec)
    {
        // walk tree to generate sparse feature vector
        maxBlockIndex = MAX_ALPHA_SIZE - 3; // highest valid index for motif
        currBlock = 0;
        currIndex = 0;
        currStack = -1;
        curr = 0;

        while (currStack >= 0 || currIndex <= maxBlockIndex)
        {
            if (currIndex == 0 && (pTree->node[currBlock].flags & TF_LEAF))
            {
                if (pTree->node[currBlock].value != 0)
                {
                    featVectorIndex[curr] = pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 1];
                    featVectorValue[curr++] = pTree->node[currBlock].value;
                    pTree->node[currBlock].value = 0;
                }
            }

            if (pTree->node[currBlock].ib.idx[currIndex] != 0)
            {
                if (pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
                {
                    if ((currIndex + 2 <= (int) (2 * pTree->node[currBlock].value - 2)) ||
                        ((currIndex + 2 == MAX_ALPHA_SIZE - 2) &&
                         (pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] != 0)))

                    {
                        blockStack[++currStack] = currIndex + 2;
                        blockStack[++currStack] = currBlock;
                    }
                }
                else
                {
                    if (currIndex + 1 <= alphaInf->numAlphabetChars + 1)
                    {
                        blockStack[++currStack] = currIndex + 1;
                        blockStack[++currStack] = currBlock;
                    }
                }

                currBlock = pTree->node[currBlock].ib.idx[currIndex];
                currIndex = 0;
            }
            else
            {
                if (pTree->node[currBlock].flags & TF_PATTERN_BLOCK)
                {
                    currIndex = currIndex + 2;

                    if (currIndex < (int) (2 * pTree->node[currBlock].value))
                        continue;

                    if (pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2] != 0)
                    {
                        // current block finished afterwards - nothing put on stack
                        currBlock = pTree->node[currBlock].ib.idx[MAX_ALPHA_SIZE - 2];
                        currIndex = 0;
                        continue;
                    }
                }
                else
                {
                    currIndex++;

                    while ((pTree->node[currBlock].ib.idx[currIndex] == 0) &&
                           (currIndex <= alphaInf->numAlphabetChars + 1))
                        currIndex++;

                    if (currIndex <= alphaInf->numAlphabetChars + 1)
                        continue;
                }

                if (currStack == -1)
                    continue;

                currBlock = blockStack[currStack--];
                currIndex = blockStack[currStack--];
            }
        }

        if (curr > 0)
        {
            featVectorIndex[curr] = MAXUINT32;
            featVectorValue[curr] = MAXINT32;
            sort2Arrays(MAXUINT32, featVectorIndex, featVectorValue, 1, numUsedFeatures, NULL);

            curr = 0;
            while (featVectorIndex[curr] != MAXUINT32)
            {
                INTEGER(slot_j)[*jIdx] = featVectorIndex[curr];

                if (normalized)
                    REAL(slot_x)[(*jIdx)++] = featVectorValue[curr++] / normValue;
                else
                    REAL(slot_x)[(*jIdx)++] = featVectorValue[curr++];
            }
        }
    }
    else
    {
        numEntries = kh_size(featMap);
        numFeatures = 0;
        curr = 0;

        struct fmData entryData;

        if (numEntries > 0)
        {
            // get key pointers
            char **keys = (char **) R_Calloc(numEntries, char *);
            pKeys = keys;
            curr = 0;

            for (iter = kh_begin(featMap); iter != kh_end(featMap); iter++)
            {
                if (kh_exist(featMap, iter))
                {
                    entryData = kh_value(featMap, iter);

                    if (entryData.sampleIndex == sampleIndex)
                    {
                        keys[curr++] = (char *)kh_key(featMap, iter);
                        numFeatures++;
                    }
                }
            }

            // sort keys
            ks_mergesort(str, numFeatures, (const char **) keys, NULL);

            // loop through keys
            for (i = 0; i < numFeatures; i++)
            {
                iter = kh_get(fim, featMap, keys[i]);

                if (iter != kh_end(featMap))
                {
                    entryData = kh_value(featMap, iter);

                    if (entryData.sampleIndex == sampleIndex && entryData.counter > 0)
                    {
                        INTEGER(slot_j)[*jIdx] = entryData.featureIndex;

                        if (normalized)
                            REAL(slot_x)[(*jIdx)++] = entryData.counter / normValue;
                        else
                            REAL(slot_x)[(*jIdx)++] = entryData.counter;
                    }
                }
            }

            R_Free(keys);
            pKeys = NULL;
        }
    }

    return;
}

void getKMStdAnnMotif(NumericMatrix km, ByteStringVector x, ByteStringVector y, int sizeX, int sizeY,
                      IntegerVector selX, IntegerVector selY, ByteStringVector annCharset,
                      ByteStringVector annX, ByteStringVector annY, ByteStringVector motifs,
                      IntegerVector motifLengths, int nodeLimit, int maxMotifLength, int maxPatternLength,
                      bool normalized, bool symmetric, bool presence, int maxSeqLength, struct alphaInfo *alphaInf)
{
    int i, freeNode, maxNoOfNodes, iX, iY, numSamples, fDim;
    int32_t *featVectorValue;
    uint32_t *featVectorIndex, *unweightedPositions, **unweightedPos;
    uint64_t keyPoolSize, poolNextFree;
    bool printWarning = TRUE;
    struct prefTreeMotif *pTree;
    struct indexBlock nullBlock;
    struct intfFindMotifs intf;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    IntegerVector unweightedPosStart(motifs.length + 1);

    // setup init mask
    for (i = 0; i < MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    numSamples = sizeX;

    if (!symmetric)
        numSamples += sizeY;

    double *normValues = (double *) R_alloc(numSamples, sizeof(double));

    // hashing only relevant for motif with annotation
    // for motif without annotation 32 bit feature index is sufficient

    // allocate arrays for sparse feature vectors
    if (annX.length > 0)
    {
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);
        fDim = 0;  // will be determined later

        // allocate key pool
        keyPoolSize = INIT_POOL_SIZE;
        intf.keyPoolSize = &keyPoolSize;
        poolNextFree = 0;
        intf.poolNextFree = &poolNextFree;
        intf.keyPool = (char *) R_Calloc(keyPoolSize, char);
        pKeypool = intf.keyPool;

        unweightedPos = &unweightedPositions;

        findUnweightedPositions(motifs, &unweightedPosStart, unweightedPos);

        intf.unweightedPosStart = &unweightedPosStart;
        intf.unweightedPos = *unweightedPos;
    }
    else
    {
        fDim = motifs.length + 1;
        intf.keyPoolSize = NULL;
        intf.poolNextFree = NULL;
        intf.keyPool = NULL;
        pKeypool = NULL;
        intf.annptr = NULL;
        intf.featMap = NULL;
    }

    // alloc mem for prefix tree
    if (nodeLimit < MAX_BLOCK)
        maxNoOfNodes = nodeLimit;
    else
        maxNoOfNodes = MAX_BLOCK;

    freeNode = 1;

    pTree = (struct prefTreeMotif *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    if (pTree == NULL)
    {
        Rprintf("Allocation of heap for tree failed\n");
        initMatrixWithNA(km, sizeX, sizeY);
        return;

    }

    //create tree for motifs
    if (!createMotifTree(motifs, maxMotifLength, pTree, maxNoOfNodes, &freeNode, &nullBlock,
                         &printWarning, alphaInf, TRUE))
    {
        Rprintf("Creation of tree failed\n");
        initMatrixWithNA(km, sizeX, sizeY);
        return;
    }

    intf.markUsedOnly = TRUE;
    intf.markMotifsInSample = FALSE;
    intf.zeroFeatures = FALSE;
    intf.presence = presence;
    intf.offset = 0;
    intf.pTree = pTree;
    intf.pErd = NULL;
    intf.pProf = NULL;
    intf.alphaInf = alphaInf;
    intf.maxMotifLength = maxMotifLength;
    intf.motifLengths = &motifLengths;

    if (annX.length > 0)
    {
        intf.annotationIndexMap = &annotationIndexMap;
        intf.maxNoOfNodes = maxNoOfNodes;
        intf.freeNode = &freeNode;
        intf.nullBlock = &nullBlock;
        intf.printWarning = &printWarning;
        intf.featMap = kh_init(fim);
        pFeatureHMap = intf.featMap;

        // calculate kernel matrix
        for (i = 0; i < numSamples; i++)
        {
            R_CheckUserInterrupt();

            if (i < sizeX)
            {
                iX = selX[i];
                intf.seqptr = x.ptr[iX];
                intf.seqnchar = x.nchar[iX];

                if (annX.length > 0)
                    intf.annptr = annX.ptr[iX];
            }
            else
            {
                iY = selY[i - sizeX];
                intf.seqptr = y.ptr[iY];
                intf.seqnchar = y.nchar[iY];

                if (annY.length > 0)
                    intf.annptr = annY.ptr[iY];
            }

            intf.featuresPerSample = 0;
            intf.rowIndex = i;

            findMotifs(&intf);

            if (intf.kernelValue == -1)
            {
                initMatrixWithNA(km, sizeX, sizeY);
                return;
            }

            if ((intf.featuresPerSample + 1) > fDim)
                fDim = intf.featuresPerSample + 1;
        }
    }

    featVectorValue = (int32_t *) R_alloc(numSamples * fDim, sizeof(int32_t));
    featVectorIndex = (uint32_t *) R_alloc(numSamples * fDim, sizeof(uint32_t));

    setFeatureIndex(pTree, maxMotifLength, maxPatternLength, alphaInf, FALSE, TRUE, motifs, &motifLengths,
                    annX.length > 0, NULL, intf.featMap, &intf.keyPool, &keyPoolSize, &poolNextFree, FALSE, FALSE);

    intf.markUsedOnly = FALSE;
    intf.markMotifsInSample = FALSE;
    intf.getKernelValue = TRUE;

    // calculate kernel matrix
    for (i = 0; i < numSamples; i++)
    {
        R_CheckUserInterrupt();

        if (i < sizeX)
        {
            iX = selX[i];
            intf.seqptr = x.ptr[iX];
            intf.seqnchar = x.nchar[iX];

            if (annX.length > 0)
                intf.annptr = annX.ptr[iX];
        }
        else
        {
            iY = selY[i - sizeX];
            intf.seqptr = y.ptr[iY];
            intf.seqnchar = y.nchar[iY];

            if (annY.length > 0)
                intf.annptr = annY.ptr[iY];
        }

        intf.rowIndex = i;

        findMotifs(&intf);

        if (intf.kernelValue == -1)
        {
            initMatrixWithNA(km, sizeX, sizeY);
            return;
        }

        if (normalized)
            normValues[i] = sqrt(intf.kernelValue);
        else
            normValues[i] = intf.kernelValue;

        R_CheckUserInterrupt();

        getNonzeroMotifs(annX.length > 0, featVectorValue, featVectorIndex, pTree, intf.featMap,
                         fDim, i, alphaInf, maxMotifLength);
    }

    if (pKeypool != NULL)
    {
        R_Free(pKeypool);
        pKeypool = NULL;
    }

    // sort feature vectors
    sort2Arrays(MAXUINT32, featVectorIndex, featVectorValue, numSamples, fDim, NULL);

    computeKernelMatrix(MAXUINT32, featVectorIndex, featVectorValue, km, normValues, fDim,
                        sizeX, sizeY, normalized, symmetric);

    return;
}

void getKMPosDistMotif(NumericMatrix km, ByteStringVector x, ByteStringVector y, int sizeX, int sizeY,
                       IntegerVector selX, IntegerVector selY, IntegerVector offsetX, IntegerVector offsetY,
                       ByteStringVector motifs, IntegerVector motifLengths, int nodeLimit, int maxMotifLength,
                       int maxPatternLength, bool normalized,  bool symmetric, bool posSpec,
                       NumericVector distWeight, int maxSeqLength, struct alphaInfo *alphaInf)
{
    uint32_t fDim, *featVectorIndex;
    int32_t *featVectorValue;
    uint64_t *featVectorsStart;
    int i, freeNode, maxNoOfNodes, iX, iY, numSamples, maxFeaturesPerSample;
    bool printWarning = TRUE;
    struct prefTreeMotif *pTree;
    struct indexBlock nullBlock;
    struct intfFindMotifs intf;

    // setup init mask
    for (i = 0; i < MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    numSamples = sizeX;

    if (!symmetric)
        numSamples += sizeY;

    double *normValues = (double *) R_Calloc(numSamples, double);

    // allocate arrays for sparse feature vectors with 32 or 64 bit index
    // store only unnormalized k-mer counts to avoid double space usage
    // add one for the sentinel
    fDim = 2 * maxSeqLength + 1;
    
    featVectorValue = (int32_t *) R_Calloc(numSamples * fDim, int32_t);
    featVectorIndex = (uint32_t *) R_Calloc(numSamples * fDim, uint32_t);
    featVectorsStart = (uint64_t *) R_Calloc(numSamples + 1, uint64_t);

    // alloc mem for prefix tree
    maxNoOfNodes = MAX_BLOCK;

    if (nodeLimit < maxNoOfNodes)
        maxNoOfNodes = nodeLimit;

    freeNode = 1;

    pTree = (struct prefTreeMotif *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    if (pTree == NULL)
    {
        Rprintf("Allocation of heap for tree failed\n");
        initMatrixWithNA(km, sizeX, sizeY);
        R_Free(featVectorIndex);
        R_Free(featVectorValue);
        R_Free(featVectorsStart);
        R_Free(normValues);
        return;
    }

    //create tree for motifs
    if (!createMotifTree(motifs, maxMotifLength, pTree, maxNoOfNodes, &freeNode, &nullBlock,
                         &printWarning, alphaInf, TRUE))
    {
        Rprintf("Creation of tree failed\n");
        initMatrixWithNA(km, sizeX, sizeY);
        R_Free(featVectorIndex);
        R_Free(featVectorValue);
        R_Free(featVectorsStart);
        R_Free(normValues);
        return;
    }

    intf.markUsedOnly = FALSE;
    intf.markMotifsInSample = FALSE;
    intf.zeroFeatures = TRUE;
    intf.presence = FALSE;
    intf.pTree = pTree;
    intf.pErd = NULL;
    intf.pProf = NULL;
    intf.alphaInf = alphaInf;
    intf.maxMotifLength = maxMotifLength;
    intf.currFeatVecLength = numSamples * fDim;
    intf.pFeatVecValue = featVectorValue;
    intf.pFeatVecIndex = featVectorIndex;
    intf.fDim = fDim;

    featVectorsStart[0] = 0;
    maxFeaturesPerSample = 0;
    intf.elemIndex = 0;

    // calculate kernel matrix
    for (i = 0; i < numSamples; i++)
    {
        R_CheckUserInterrupt();
        intf.rowIndex = i;
        intf.offset = 0;

        if (i < sizeX)
        {
            iX = selX[i];
            intf.seqptr = x.ptr[iX];
            intf.seqnchar = x.nchar[iX];

            if (offsetX.length() > 0)
                intf.offset = offsetX[iX];
        }
        else
        {
            iY = selY[i - sizeX];
            intf.seqptr = y.ptr[iY];
            intf.seqnchar = y.nchar[iY];

            if (offsetY.length() > 0)
                intf.offset = offsetY[iY];
        }

        findMotifsForPos(&intf);

        if (intf.kernelValue == -1)
        {
            initMatrixWithNA(km, sizeX, sizeY);
            return;
        }

        featVectorsStart[i + 1] = intf.elemIndex;

        if (maxFeaturesPerSample < (int) (featVectorsStart[i + 1] - featVectorsStart[i]))
            maxFeaturesPerSample = featVectorsStart[i + 1] - featVectorsStart[i];

        if (normalized)
            normValues[i] = sqrt(intf.kernelValue);
        else
            normValues[i] = intf.kernelValue;
    }

    featVectorIndex = intf.pFeatVecIndex;
    featVectorValue = intf.pFeatVecValue;

    computeKernelMatrixPos(MAXUINT32, featVectorIndex, featVectorValue, featVectorsStart, km,
                           normValues, maxFeaturesPerSample, motifs.length, sizeX, sizeY, normalized, symmetric,
                           FALSE, distWeight);
    
    R_Free(featVectorIndex);
    R_Free(featVectorValue);
    R_Free(featVectorsStart);
    R_Free(normValues);

    return;
}

void getERDMotif(NumericMatrix erd, ByteStringVector x, int sizeX, IntegerVector selX, ByteStringVector annCharset,
                 ByteStringVector annX, struct intfFindMotifs *intf, ByteStringVector motifs,
                 IntegerVector *motifLengths, int maxPatternLength, bool normalized, uint64_t *dimFeatureSpace,
                 bool useHash, bool useRowNames, bool useColNames)
{
    int i, iX, j, numProtect;
    double *normValues;
    const void *vmax;
    SEXP rownames, colnames, dimnames;

    intf->pErd = &erd;

    if (useColNames)
        PROTECT(colnames = Rf_allocVector(STRSXP, intf->numUsedMotifs));
    else
        PROTECT(colnames = Rf_allocVector(STRSXP, 0));

    PROTECT(rownames = Rf_allocVector(STRSXP, 0));
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, rownames);
    SET_VECTOR_ELT(dimnames, 1, colnames);
    Rf_setAttrib(erd, R_DimNamesSymbol, dimnames);
    numProtect = 3;

    vmax = vmaxget();

    setFeatureIndex(intf->pTree, intf->maxMotifLength, maxPatternLength, intf->alphaInf, useColNames && (intf->numUsedMotifs > 0),
                    intf->zeroFeatures, motifs, motifLengths, annX.length > 0, colnames, intf->featMap, NULL, NULL, NULL, FALSE,
                    FALSE);

    normValues = (double *) R_alloc(sizeX, sizeof(double));

    intf->markUsedOnly = FALSE;
    intf->markMotifsInSample = FALSE;
    intf->getKernelValue = normalized;

    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();

        iX = selX[i];
        intf->rowIndex = i;
        intf->seqptr = x.ptr[iX];
        intf->seqnchar = x.nchar[iX];

        if (annX.length > 0)
            intf->annptr = annX.ptr[iX];

        findMotifs(intf);

        if (normalized)
            normValues[i] = sqrt(intf->kernelValue);
    }

    if (normalized)
    {
        for (i = 0; i < sizeX; i++)
        {
            if (normValues[i] > 0)
            {
                for (j = 0; j < intf->numUsedMotifs; j++)
                    erd(i, j) = erd(i, j) / normValues[i];
            }
        }
    }

    vmaxset(vmax);

    if (numProtect > 0)
        UNPROTECT(numProtect);

    return;
}

void getERSMotif(SEXP *pErs, ByteStringVector x, int sizeX, IntegerVector selX, ByteStringVector annCharset,
                 ByteStringVector annX, struct intfFindMotifs *intf, ByteStringVector motifs,
                 IntegerVector *motifLengths, int maxPatternLength, bool normalized, uint64_t *dimFeatureSpace,
                 bool useHash, bool useRowNames, bool useColNames)
{
    int i, iX, numProtect, jIdx;
    int32_t *featVectorValue;
    uint32_t *featVectorIndex;
    double normValue = 1;
    const void *vmax;
    SEXP rownames, colnames, dimnames, ers, ers1, dims, dims1, slot_j, slot_x, slot_p;

    // allocate explicit representation
    numProtect = 0;
    ers = PROTECT(NEW_OBJECT(MAKE_CLASS("ExplicitRepresentationSparse")));
    *pErs = ers;
    dims = PROTECT(Rf_allocVector(INTSXP, 2));
    SET_SLOT(ers, Rf_mkChar("Dim"), dims);
    INTEGER(dims)[0] = sizeX;
    INTEGER(dims)[1] = intf->numUsedMotifs;
    slot_p = PROTECT(Rf_allocVector(INTSXP, sizeX + 1));
    SET_SLOT(ers, Rf_mkChar("p"), slot_p);
    numProtect = 3;

    if (useRowNames || useColNames)
    {
        PROTECT(rownames = Rf_allocVector(STRSXP, 0));
        PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
        numProtect +=2;

        if (useColNames && intf->numUsedMotifs > 0)
        {
            PROTECT(colnames = Rf_allocVector(STRSXP, intf->numUsedMotifs));
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

    slot_j = PROTECT(Rf_allocVector(INTSXP, intf->numNonzeroFeatures));
    SET_SLOT(ers, Rf_mkChar("j"), slot_j);

    slot_x = PROTECT(Rf_allocVector(REALSXP, intf->numNonzeroFeatures));
    SET_SLOT(ers, Rf_mkChar("x"), slot_x);
    numProtect = numProtect + 2;

    vmax = vmaxget();

    setFeatureIndex(intf->pTree, intf->maxMotifLength, maxPatternLength, intf->alphaInf, (useColNames && intf->numUsedMotifs > 0),
                    intf->zeroFeatures, motifs, motifLengths, annX.length > 0, colnames, intf->featMap, NULL, NULL, NULL, FALSE,
                    FALSE);

    // set feature values in ER
    intf->markUsedOnly = FALSE;
    intf->markMotifsInSample = FALSE;
    intf->getKernelValue = TRUE;

    featVectorValue = (int32_t *) R_alloc(intf->numUsedMotifs + 1, sizeof(int32_t));
    featVectorIndex = (uint32_t *) R_alloc(intf->numUsedMotifs + 1, sizeof(uint32_t));

    jIdx = 0;
    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();

        INTEGER(slot_p)[i] = jIdx;

        iX = selX[i];
        intf->rowIndex = i;
        intf->seqptr = x.ptr[iX];
        intf->seqnchar = x.nchar[iX];

        if (annX.length > 0)
            intf->annptr = annX.ptr[iX];

        findMotifs(intf);

        if (intf->kernelValue > 0)
        {
            if (normalized)
                normValue = sqrt(intf->kernelValue);

            R_CheckUserInterrupt();

            getNonzeroMotifsERS(annX.length > 0, intf->pTree, intf->featMap, intf->alphaInf, intf->maxMotifLength, i,
                                intf->numUsedMotifs + 1, slot_j, slot_x, &jIdx, featVectorIndex, featVectorValue, normValue,
                                normalized, intf->zeroFeatures);
        }
        else
        {
            if (intf->kernelValue == -1)
            {
                // return empty ers instead
                UNPROTECT(numProtect);
                vmaxset(vmax);

                Rprintf("Determination of feature values failed\n");
                ers1 = PROTECT(NEW_OBJECT(MAKE_CLASS("ExplicitRepresentationSparse")));
                *pErs = ers1;
                dims1 = PROTECT(Rf_allocVector(INTSXP, 2));
                SET_SLOT(ers1, Rf_mkChar("Dim"), dims1);
                INTEGER(dims1)[0] = 0;
                INTEGER(dims1)[1] = 0;
                UNPROTECT(2);
                return;
            }
        }
    }

    INTEGER(slot_p)[sizeX] = jIdx;

    vmaxset(vmax);
    UNPROTECT(numProtect);

    return;
}

RcppExport SEXP genExplRepMotif(ByteStringVector x, int sizeX, IntegerVector selX, ByteStringVector annCharset,
                                ByteStringVector annX, int maxSeqLength, int bioCharset, ByteStringVector motifs,
                                IntegerVector motifLengths, int maxMotifLength, int maxPatternLength,
                                int nodeLimit, bool presence, bool normalized, bool unmapped, bool lowercase,
                                bool useRowNames, bool useColNames, bool zeroFeatures, bool sparse)
{
    int i, iX, freeNode, maxNoOfNodes;
    uint32_t *unweightedPositions, **unweightedPos;
    uint64_t poolNextFree, keyPoolSize, dimFeatureSpace;
    bool useHash, printWarning = TRUE;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    struct prefTreeMotif *pTree;
    struct indexBlock nullBlock;
    struct intfFindMotifs intf;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    IntegerVector unweightedPosStart(motifs.length + 1);
    SEXP explicitRepSparse;

    pUnweightedPos = NULL;
    pKeys = NULL;

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    if (annX.length > 0)
    {
        useHash = TRUE;
        dimFeatureSpace = 0;

        for (i = 0; i < motifs.length; i++)
            dimFeatureSpace += pow(annCharset.nchar[0], motifLengths[i]);

        // allocate key pool
        keyPoolSize = INIT_POOL_SIZE;
        intf.keyPoolSize = &keyPoolSize;
        poolNextFree = 0;
        intf.poolNextFree = &poolNextFree;
        intf.keyPool = (char *) R_Calloc(keyPoolSize, char);
        pKeypool = intf.keyPool;
        unweightedPos = &unweightedPositions;
        
        findUnweightedPositions(motifs, &unweightedPosStart, unweightedPos);
        
        intf.unweightedPosStart = &unweightedPosStart;
        intf.unweightedPos = *unweightedPos;
    }
    else
    {
        useHash = FALSE;
        dimFeatureSpace = motifs.length;
        intf.keyPoolSize = NULL;
        intf.keyPool = NULL;
        intf.poolNextFree = NULL;
        pKeypool = NULL;
    }

    // check if features space to large
    if (dimFeatureSpace > MAXINDEX32 && zeroFeatures)
        return(generateEmptyExplicitRep(sizeX, sparse));

    if (dimFeatureSpace > FEATURE_SPACE_LIMIT)
    {
        Rprintf("feature space too large\n");
        return(generateEmptyExplicitRep(sizeX, sparse));
    }

    // setup init mask
    for (i = 0; i < MAX_ALPHA_SIZE; i++)
        nullBlock.idx[i] = 0;

    if (annX.length > 0)
        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);

    if (useHash)
    {
        intf.featMap = kh_init(fim);
        pFeatureHMap = intf.featMap;
    }
    else
    {
        intf.featMap = NULL;
        pFeatureHMap = NULL;
    }

    // alloc mem for prefix tree
    maxNoOfNodes = MAX_BLOCK;

    if (nodeLimit < maxNoOfNodes)
        maxNoOfNodes = nodeLimit;

    freeNode = 1;

    pTree = (struct prefTreeMotif *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    if (pTree == NULL)
    {
        Rprintf("Allocation of heap for tree failed\n");
        return(generateEmptyExplicitRep(sizeX, sparse));
    }

    //create tree for motifs
    if (!createMotifTree(motifs, maxMotifLength, pTree, maxNoOfNodes, &freeNode, &nullBlock,
                         &printWarning, &alphaInf, FALSE))
    {
        Rprintf("Creation of tree failed\n");
        return(generateEmptyExplicitRep(sizeX, sparse));
    }

    intf.presence = presence;
    intf.zeroFeatures = zeroFeatures;
    intf.pTree = pTree;
    intf.pErd = NULL;
    intf.pProf = NULL;
    intf.alphaInf = &alphaInf;
    intf.maxMotifLength = maxMotifLength;
    intf.numUsedMotifs = 0;
    intf.numNonzeroFeatures = 0;
    intf.allIndexMaps = &allIndexMaps;
    intf.motifLengths = &motifLengths;
    intf.markUsedOnly = FALSE;
    intf.markMotifsInSample = FALSE;

    if (annX.length > 0)
    {
        intf.annotationIndexMap = &annotationIndexMap;
        intf.maxNoOfNodes = maxNoOfNodes;
        intf.freeNode = &freeNode;
        intf.nullBlock = &nullBlock;
        intf.printWarning = &printWarning;
    }
    else
        intf.annptr = NULL;

    // get no of nonzero elements and used features
    if(!zeroFeatures || sparse)
    {
        intf.markUsedOnly = TRUE;

        if (sparse)
            intf.markMotifsInSample = TRUE;

        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            iX = selX[i];
            intf.rowIndex = i;
            intf.seqptr = x.ptr[iX];
            intf.seqnchar = x.nchar[iX];

            if (annX.length > 0)
                intf.annptr = annX.ptr[iX];

            findMotifs(&intf);

            if (intf.kernelValue == -1)
            {
                Rprintf("Determination of used features failed\n");
                return(generateEmptyExplicitRep(sizeX, sparse));
            }

            if (sparse && annX.length == 0)
                resetInfoInTree(TF_RESET_FLAG, TF_MOTIF_IN_SAMPLE, pTree, maxMotifLength, &alphaInf, zeroFeatures);
        }

        if (zeroFeatures)
            intf.numUsedMotifs = dimFeatureSpace;
    }
    else
        intf.numUsedMotifs = dimFeatureSpace;

    // check if no features found
    if (intf.numUsedMotifs < 1)
        return(generateEmptyExplicitRep(sizeX, sparse));

    if (!sparse)
    {
        NumericMatrix erd(sizeX, intf.numUsedMotifs);

        getERDMotif(erd, x, sizeX, selX, annCharset, annX, &intf, motifs, &motifLengths,
                    maxPatternLength, normalized, &dimFeatureSpace, useHash, useRowNames, useColNames);

        return(erd);
    }
    else
    {
        getERSMotif(&explicitRepSparse, x, sizeX, selX, annCharset, annX, &intf, motifs, &motifLengths,
                    maxPatternLength, normalized, &dimFeatureSpace, useHash, useRowNames, useColNames);

        return(explicitRepSparse);
    }
}

RcppExport SEXP validateMotifsC(SEXP motifsR, SEXP motifLengthsR)
{
    ByteStringVector motifs = charVector2ByteStringVec(motifsR);
    IntegerVector motifLengths(motifLengthsR);
    IntegerVector errorPositions(motifs.length);
    char substitionGroupChars[256];
    int state, i, j, k, patternLength, groupStart;
    bool errorFound, patternError;
    SEXP result;

    errorFound = FALSE;

    for (i = 0; i < motifs.length; i++)
    {
        patternLength = 0;
        state = TF_CHAR_MATCH;
        patternError = FALSE;
        errorPositions[i] = -1;
        groupStart = -1;

        for (j = 0; j < motifs.nchar[i]; j++)
        {
            switch (motifs.ptr[i][j])
            {
                case '[':
                {
                    if (state == TF_CHAR_MATCH)
                    {
                        patternLength++;
                        state = TF_GROUP_MATCH;
                        groupStart = j;

                        for (k = 0; k < MAX_CHAR; k++)
                            substitionGroupChars[k] = 0;
                    }
                    else
                        patternError = TRUE;

                    break;
                }

                case '^':
                {
                    if (state == TF_GROUP_MATCH)
                    {
                        if (j == (groupStart + 1))
                            state = TF_NEG_GROUP_MATCH;
                        else
                            patternError = TRUE;
                    }
                    else
                        patternError = TRUE;

                    break;
                }

                case ']':
                {
                    if (state > TF_CHAR_MATCH)
                        state = TF_CHAR_MATCH;
                    else
                        patternError = TRUE;

                    break;
                }

                case '.':
                {
                    if (state == TF_CHAR_MATCH)
                        patternLength++;
                    else
                        patternError = TRUE;

                    break;
                }

                default:
                {
                    if (state > TF_CHAR_MATCH)
                    {
                        if (substitionGroupChars[(int)motifs.ptr[i][j]] == 0)
                            substitionGroupChars[(int)motifs.ptr[i][j]] = 1;
                        else
                            patternError = TRUE;
                    }
                    else
                        patternLength++;

                    break;
                }
            }

            if (patternError)
            {
                errorFound = TRUE;
                errorPositions[i] = j + 1;   // R adapted position
                break;
            }
        }

        if (!errorFound && (state > TF_CHAR_MATCH))
        {
            errorFound = TRUE;
            errorPositions[i] = motifs.nchar[i];
        }
        else
            motifLengths[i] = patternLength;
    }

    result = PROTECT(Rf_allocVector(LGLSXP, 1));

    if (errorFound)
    {
        for (i = 0; i < motifs.length; i++)
            motifLengths[i] = errorPositions[i];

        LOGICAL(result)[0] = FALSE;
    }
    else
        LOGICAL(result)[0] = TRUE;

    UNPROTECT(1);

    return(result);
}

RcppExport SEXP motifKernelMatrixC(SEXP xR, SEXP yR, SEXP selXR, SEXP selYR, SEXP sizeXR, SEXP sizeYR,
                                   SEXP isXStringSetR, SEXP symmetricR, SEXP offsetXR, SEXP offsetYR,
                                   SEXP annCharsetR, SEXP annXR, SEXP annYR, SEXP motifsR, SEXP motifLengthsR,
                                   SEXP nodeLimitR, SEXP maxMotifLengthR, SEXP maxPatternLengthR, SEXP bioCharsetR,
                                   SEXP ignoreLowerR, SEXP unmappedR, SEXP maxSeqLengthR, SEXP posSpecR,
                                   SEXP distWeightR, SEXP normalizedR, SEXP presenceR)
{
    int sizeX = as<int>(sizeXR);
    int sizeY = as<int>(sizeYR);
    int i;
    uint64_t dimFeatureSpace;
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
    ByteStringVector motifs = charVector2ByteStringVec(motifsR);
    ByteStringVector x, y, annX, annY, annCharset;
    IntegerVector selX(selXR);
    IntegerVector selY(selYR);
    IntegerVector motifLengths(motifLengthsR);
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

    int bioCharset = as<int>(bioCharsetR);
    int maxMotifLength = as<int>(maxMotifLengthR);
    int maxPatternLength = as<int>(maxPatternLengthR);
    int nodeLimit = as<int>(nodeLimitR);
    int maxSeqLength = as<int>(maxSeqLengthR);
    bool lowercase = !as<bool>(ignoreLowerR);
    bool posSpec = as<bool>(posSpecR);
    bool unmapped = as<bool>(unmappedR);
    bool normalized = as<bool>(normalizedR);
    bool presence = as<bool>(presenceR);

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    if (annX.length > 0)
    {
        dimFeatureSpace = 0;

        for (i = 0; i < motifs.length; i++)
            dimFeatureSpace += pow(annCharset.nchar[0], motifLengths[i]);
    }
    else
        dimFeatureSpace = motifs.length;

    if (posSpec || distWeight.length() > 0)
    {
        getKMPosDistMotif(km, x, y, sizeX, sizeY, selX, selY, offsetX, offsetY, motifs,
                          motifLengths, nodeLimit, maxMotifLength, maxPatternLength, normalized,
                          symmetric, posSpec, distWeight, maxSeqLength, &alphaInf);
    }
    else
    {
        getKMStdAnnMotif(km, x, y, sizeX, sizeY, selX, selY, annCharset, annX, annY, motifs,
                         motifLengths, nodeLimit, maxMotifLength, maxPatternLength, normalized,
                         symmetric, presence, maxSeqLength, &alphaInf);
    }

    vmaxset(vmax);

    return(km);
}

// no support for annotation because annotation part cannot be represented as index
uint64_t * featureNamesToIndexMotif(SEXP featureNames, int numFeatures, void **pMotifTree, int *freeNode,
                                    ByteStringVector motifs, IntegerVector motifLengths, int maxMotifLength,
                                    int maxPatternLength, int nodeLimit, struct alphaInfo *alphaInf)
{
    int i, maxNoOfNodes;
    uint64_t *featIndex;
    bool printWarning = TRUE;
    const void *vmax;
    struct indexBlock nullBlock;
    struct prefTreeMotif *pTree;
    struct intfStorePattern intf;

    pTree = (struct prefTreeMotif *) *pMotifTree;
    memset(nullBlock.idx, 0, sizeof(indexBlock));

    vmax = vmaxget();

    featIndex = (uint64_t *) R_alloc(numFeatures, sizeof(uint64_t));

    if (nodeLimit < MAX_BLOCK)
        maxNoOfNodes = nodeLimit;
    else
        maxNoOfNodes = MAX_BLOCK;

    if (pTree == NULL)
    {
        *freeNode = 1;

        pTree = (struct prefTreeMotif *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));
        *pMotifTree = (void *) pTree;

        // create motif tree
        if (!createMotifTree(motifs, maxMotifLength, pTree, maxNoOfNodes, freeNode, &nullBlock,
                             &printWarning, alphaInf, TRUE))
        {
            Rprintf("Creation of motif tree failed\n");
            vmaxset(vmax);
            return(NULL);
        }
    }

    intf.pTree = pTree;
    intf.alphaInf = alphaInf;
    intf.annSpec = FALSE;
    intf.zeroFeatures = TRUE;

    for (i = 0; i < numFeatures; i++)
    {
        intf.pattern = (char *) CHAR(STRING_ELT(featureNames, i));
        intf.patternLength = strlen(intf.pattern);

        if (!getLeaf(&intf))
        {
            Rprintf("Traversal of motif tree failed\n");
            vmaxset(vmax);
            return(NULL);
        }

        // motif not found
        if (intf.leafBlock == 0)
            featIndex[i] = NA_INTEGER;
        else
            featIndex[i] = pTree->node[intf.leafBlock].ib.idx[MAX_ALPHA_SIZE - 2];
    }

    return(featIndex);
}

// delivers the position specific unnormalized feature vectors in compressed format together with kernel values
// no support for annotation in this function
void genFeatVectorsMotif(ByteStringVector x, int sizeX, IntegerVector selX, IntegerVector offsetX,
                         int maxSeqLength, void **pMotifTree, int *freeNode, ByteStringVector motifs,
                         IntegerVector motifLengths, int maxMotifLength, int maxPatternLength,
                         int nodeLimit, struct alphaInfo *alphaInf, bool presence, bool normalized,
                         bool posSpecific, NumericVector distWeight, int sortType, uint64_t **startIndex,
                         uint32_t **featVectorIndex, int32_t **featVectorValue, double **kernelValue)
{
    int i, maxNoOfNodes;
    uint32_t maxNoElements;
    bool printWarning = TRUE;
    const void *vmax;
    struct indexBlock nullBlock;
    struct prefTreeMotif *pTree;
    struct intfFindMotifs intf;
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    IntegerVector selCurr(1), selY(0), offsetY(0);
    NumericMatrix kmOne(1,1);
    ByteStringVector y;
    
    y.length = 0;
    pTree = (struct prefTreeMotif *) *pMotifTree;
    memset(nullBlock.idx, 0, sizeof(indexBlock));

    vmax = vmaxget();

    if (nodeLimit < MAX_BLOCK)
        maxNoOfNodes = nodeLimit;
    else
        maxNoOfNodes = MAX_BLOCK;

    if (pTree == NULL)
    {
        *freeNode = 1;

        pTree = (struct prefTreeMotif *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));
        *pMotifTree = (void *) pTree;

        // create motif tree
        if (!createMotifTree(motifs, maxMotifLength, pTree, maxNoOfNodes, freeNode, &nullBlock,
                             &printWarning, alphaInf, TRUE))
        {
            Rprintf("Creation of motif tree failed\n");
            vmaxset(vmax);
            return;
        }
    }

    intf.currFeatVecLength = sizeX * maxSeqLength * 2;

    if (intf.currFeatVecLength > MAX_FEAT_VEC_LENGTH)
        intf.currFeatVecLength = MAX_FEAT_VEC_LENGTH;

    *featVectorIndex = (uint32_t *) R_Calloc(intf.currFeatVecLength, uint32_t);
    *featVectorValue = (int32_t *) R_Calloc(intf.currFeatVecLength, int32_t);
    *startIndex = (uint64_t *) R_Calloc(sizeX + 1, uint64_t);

    if (normalized)
        *kernelValue = (double *) R_Calloc(sizeX, double);

    intf.markUsedOnly = FALSE;
    intf.markMotifsInSample = FALSE;
    intf.zeroFeatures = TRUE;
    intf.presence = FALSE;
    intf.pTree = pTree;
    intf.pErd = NULL;
    intf.pProf = NULL;
    intf.alphaInf = alphaInf;
    intf.maxMotifLength = maxMotifLength;
    intf.pFeatVecValue = *featVectorValue;
    intf.pFeatVecIndex = *featVectorIndex;
    intf.fDim = motifs.length;
    intf.elemIndex = 0;
    maxNoElements = 0;

    (*startIndex)[0] = 0;

    // calculate kernel matrix
    for (i = 0; i < sizeX; i++)
    {
        R_CheckUserInterrupt();
        intf.rowIndex = i;
        intf.offset = 0;
        (*startIndex)[i] = intf.elemIndex;

        if (intf.elemIndex > maxNoElements)
            maxNoElements = intf.elemIndex;

        intf.seqptr = x.ptr[selX[i]];
        intf.seqnchar = x.nchar[selX[i]];

        if (offsetX.length() > 0)
            intf.offset = offsetX[selX[i]];

        findMotifsForPos(&intf);

        if (intf.kernelValue == -1)
        {
            Rprintf("Error in generating sparse feature vectors");
            return;
        }

        (*startIndex)[i + 1] = intf.elemIndex;

        if (normalized)
        {
            if (distWeight.length() == 0)
                (*kernelValue)[i] = intf.kernelValue;
            else
            {
                selCurr[0] = selX[i];
                int currSeqLength = intf.seqnchar;
                
                getKMPosDistMotif(kmOne, x, y, 1, 1, selCurr, selY, offsetX, offsetY, motifs,
                                  motifLengths, nodeLimit, maxMotifLength, maxPatternLength, FALSE,
                                  TRUE, FALSE, distWeight, currSeqLength, alphaInf);
                (*kernelValue)[i] = kmOne(0,0);

            }
        }
    }

    (*startIndex)[sizeX] = intf.elemIndex;
    *featVectorIndex = intf.pFeatVecIndex;
    *featVectorValue = intf.pFeatVecValue;
    
    if (intf.elemIndex > maxNoElements)
        maxNoElements = intf.elemIndex;

    if (sortType == KBS_SORT_BY_POSITION)
    {
        sort2Arrays((int32_t) (MAXUINT32 >> 1), *featVectorValue, *featVectorIndex, 1,
                    maxNoElements, *startIndex);
    }

    return;
}

template<typename T>
bool getSVWeightsFeatMotif(T maxUnSignedIndex, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, ByteStringVector x,
                           int sizeX, IntegerVector selX, IntegerVector offsetX, int maxSeqLength, NumericVector coefs,
                           bool posSpecific, NumericVector distWeight, double weightLimit, ByteStringVector motifs,
                           IntegerVector motifLengths, int maxMotifLength, int maxPatternLength, int nodeLimit,
                           int minPos, int maxPos, struct alphaInfo *alphaInf, bool normalized, uint64_t *numKeys,
                           T **keys)
{
    int i, j, freeNode, result;
    uint32_t *featVectorIndex;
    int32_t *featVectorValue;
    uint64_t *startIndex, numEntries, key;
    double normFactor, limit, *kernelValue;
    void *pMotifTree;
    khiter_t iter;
    IntegerVector sel(1);

    normFactor = 1;
    pMotifTree = NULL;

    for (i = 0; i < sizeX; i++)
    {
        if (i % USER_INTERRUPT_LIMIT == 0)
            R_CheckUserInterrupt();

        sel[0] = selX[i]; // pass samples individually

        genFeatVectorsMotif(x, 1, sel, offsetX, maxSeqLength, &pMotifTree, &freeNode, motifs, motifLengths,
                            maxMotifLength, maxPatternLength, nodeLimit, alphaInf, FALSE, normalized,
                            posSpecific, distWeight, KBS_UNSORTED, &startIndex, &featVectorIndex, &featVectorValue,
                            &kernelValue);

        if (normalized)
            normFactor = 1.0 / sqrt(kernelValue[0]);

        for (j = 0; j < (int) startIndex[1]; j++)
        {
            iter = kh_put(pdfi, pdfimap, featVectorIndex[j], &result);
            key = (featVectorValue[j] - minPos) * motifs.length + featVectorIndex[j];
            iter = kh_put(pdfw, pdfwmap, key, &result);

            if (result)
                kh_value(pdfwmap, iter) = normFactor * coefs[sel[0]];
            else
                kh_value(pdfwmap, iter) = kh_value(pdfwmap, iter) + normFactor * coefs[sel[0]];
        }

        R_Free(featVectorIndex);
        R_Free(featVectorValue);
        R_Free(startIndex);

        if (normalized)
            R_Free(kernelValue);
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
void getWeightedFeatOfSVMotif(T maxUnSignedIndex, SEXP **pdFeatWeights, khash_t(pdfw) *pdfwmap,
                              khash_t(pdfi) *pdfimap, ByteStringVector x, int sizeX, IntegerVector selX,
                              IntegerVector offsetX, int maxSeqLength, NumericVector coefs, bool posSpecific,
                              NumericVector distWeight, double weightLimit, ByteStringVector motifs,
                              IntegerVector motifLengths, int maxMotifLength, int maxPatternLength,
                              int nodeLimit, int minPos, int maxPos, struct alphaInfo *alphaInf,
                              bool normalized, uint64_t *numKeys, T **keys)
{
    int i, j, row, numProtect;
    char kmer[maxPatternLength + 1], position[12];
    uint64_t featIndex;
    khiter_t iter;
    SEXP rownames, colnames, dimnames, slot_p, slot_i, slot_x, dims;

    if (!getSVWeightsFeatMotif(maxUnSignedIndex, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength, coefs,
                               posSpecific, distWeight, weightLimit, motifs, motifLengths, maxMotifLength,
                               maxPatternLength, nodeLimit, minPos, maxPos, alphaInf, normalized, numKeys, keys))
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

    for (i = minPos; i <= maxPos; i++)
    {
        snprintf(position, 12, "%d", i);
        SET_STRING_ELT(colnames, i - minPos, Rf_mkChar(position));
    }

    for (iter = kh_begin(pdfimap); iter != kh_end(pdfimap); iter++)
    {
        if (kh_exist(pdfimap, iter))
        {
            row = kh_value(pdfimap, iter);
            featIndex = kh_key(pdfimap, iter);
            kmer[motifs.nchar[featIndex]] = '\0';

            for (j = 0; j < motifs.nchar[featIndex]; j++)
                kmer[j] = motifs.ptr[featIndex][j];

            SET_STRING_ELT(rownames, row, Rf_mkChar(kmer));
        }
    }

    UNPROTECT(numProtect);

    return;
}

void getFeaturesOfSVMotif(SEXP **pdFeatWeights, khash_t(pdfw) *pdfwmap, khash_t(pdfi) *pdfimap, ByteStringVector x,
                          int sizeX, IntegerVector selX, IntegerVector offsetX, int maxSeqLength, NumericVector coefs,
                          bool posSpecific, NumericVector distWeight, double weightLimit, ByteStringVector motifs,
                          IntegerVector motifLengths, int maxMotifLength, int maxPatternLength, int nodeLimit,
                          int minPos, int maxPos, uint64_t dimFeatureSpace, struct alphaInfo *alphaInf,
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
            getWeightedFeatOfSVMotif(maxUIndex8, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                     coefs, posSpecific, distWeight, weightLimit, motifs, motifLengths, maxMotifLength,
                                     maxPatternLength, nodeLimit, minPos, maxPos, alphaInf, normalized, numKeys,
                                     (uint8_t **) keys);
            break;
        }

        case 2:
        {
            getWeightedFeatOfSVMotif(maxUIndex16, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                     coefs, posSpecific, distWeight, weightLimit, motifs, motifLengths, maxMotifLength,
                                     maxPatternLength, nodeLimit, minPos, maxPos, alphaInf, normalized, numKeys,
                                     (uint16_t **) keys);
            break;
        }

        case 3:
        case 4:
        {
            getWeightedFeatOfSVMotif(maxUIndex32, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                     coefs, posSpecific, distWeight, weightLimit, motifs, motifLengths, maxMotifLength,
                                     maxPatternLength, nodeLimit, minPos, maxPos, alphaInf, normalized, numKeys,
                                     (uint32_t **) keys);
            break;
        }

        default:
        {
            getWeightedFeatOfSVMotif(maxUIndex64, pdFeatWeights, pdfwmap, pdfimap, x, sizeX, selX, offsetX, maxSeqLength,
                                     coefs, posSpecific, distWeight, weightLimit, motifs, motifLengths, maxMotifLength,
                                     maxPatternLength, nodeLimit, minPos, maxPos, alphaInf, normalized, numKeys,
                                     (uint64_t **) keys);
            break;
        }
    }

    return;
}

RcppExport bool featuresToMotifTree(ByteStringVector motifs, int maxMotifLength, struct alphaInfo *alphaInf,
                                    bool annSpec, ByteStringVector annCharset, struct prefTreeMotif **pTree,
                                    int *freeNode, int nodeLimit, struct indexBlock *nullBlock, int maxNoOfNodes,
                                    bool *printWarning, IntegerVector *unweightedPosStart, uint32_t **unweightedPos,
                                    struct allIndMaps *allIndexMaps)
{
    *freeNode = 1;

    *pTree = (struct prefTreeMotif *) R_alloc(maxNoOfNodes, sizeof(struct treeNode));

    if (*pTree == NULL)
    {
        Rprintf("Allocation of heap for tree failed\n");
        return(FALSE);
    }

    // create motif tree
    if (!createMotifTree(motifs, maxMotifLength, *pTree, maxNoOfNodes, freeNode, nullBlock,
                         printWarning, alphaInf, FALSE))
    {
        Rprintf("Creation of tree failed\n");
        return(FALSE);
    }

    findUnweightedPositions(motifs, unweightedPosStart, unweightedPos);

    return(TRUE);
}

void genPredProfileMotif(NumericMatrix profiles, ByteStringVector x, IntegerVector selX, int numSamples,
                         ByteStringVector annCharset, ByteStringVector annX, int maxSeqLength, bool unmapped,
                         int kernelType, int k, int m, int bioCharset, NumericMatrix featureWeights,
                         int svmIndex, ByteStringVector motifs, IntegerVector *motifLengths, int maxMotifLength,
                         int maxPatternLength, ByteStringVector fwMotifs, IntegerVector *fwMotifLengths,
                         int fwMaxMotifLength, int fwMaxPatternLength, int nodeLimit, bool lowercase,
                         bool normalized, bool presence)
{
    int i, j, iX, freeNode, maxNoOfNodes, currMaxMotifLength, currSeqLength;
    uint32_t *unweightedPositions, **unweightedPos;
    uint64_t keyPoolSize, poolNextFree, dimFeatureSpace;
    double normValue;
    bool useHash, printWarning = TRUE;
    char *keyPool;
    struct prefTreeMotif *pTree;
    struct intfFindMotifs intf;
    struct indexBlock nullBlock;
    struct alphaInfo alphaInf;
    struct allIndMaps allIndexMaps;
    IntegerVector *currMotifLengths, selCurr(1), selY(0);
    IntegerVector annotationIndexMap(MAX_CHAR);
    IntegerVector reverseAnnotationMap(MAX_CHAR);
    IntegerVector unweightedPosStart(motifs.length + 1);
    NumericMatrix kmOne(1,1);
    ByteStringVector *currMotifs, y, annY;

    y.length = 0;
    annY.length = 0;
    memset(nullBlock.idx, 0, sizeof(indexBlock));

    pUnweightedPos = NULL;
    pKeys = NULL;
    unweightedPos = &unweightedPositions;

    // for ann spec kernel with original motifs
    // without annotation with motifs from feature names
    if (!(annX.length > 0) && !normalized)
    {
        currMotifs = &fwMotifs;
        currMotifLengths = fwMotifLengths;
        currMaxMotifLength = fwMaxMotifLength;
    }
    else
    {
        currMotifs = &motifs;
        currMotifLengths = motifLengths;
        currMaxMotifLength = maxMotifLength;
    }

    if (annX.length > 0)
    {
        useHash = TRUE;
        dimFeatureSpace = 0;

        for (i = 0; i < currMotifs->length; i++)
            dimFeatureSpace += pow(annCharset.nchar[0], (*currMotifLengths)[i]);

        initAnnotationMaps(annCharset, &annotationIndexMap, &reverseAnnotationMap);

        // allocate key pool
        keyPoolSize = INIT_POOL_SIZE;
        intf.keyPoolSize = &keyPoolSize;
        poolNextFree = 0;
        intf.poolNextFree = &poolNextFree;
        pKeypool = (char *) R_Calloc(keyPoolSize, char);
        keyPool = pKeypool;
        intf.keyPool = keyPool;
    }
    else
    {
        useHash = FALSE;
        dimFeatureSpace = currMotifs->length;
    }

    if (useHash)
    {
        intf.featMap = kh_init(fim);
        pFeatureHMap = intf.featMap;
    }
    else
    {
        intf.featMap = NULL;
        pFeatureHMap = NULL;
    }

    if (nodeLimit < MAX_BLOCK)
        maxNoOfNodes = nodeLimit;
    else
        maxNoOfNodes = MAX_BLOCK;

    getAlphabetInfo(bioCharset, lowercase, unmapped, &alphaInf, &allIndexMaps);

    if (!featuresToMotifTree(motifs, maxMotifLength, &alphaInf, annX.length > 0,
                             annCharset, &pTree, &freeNode, nodeLimit, &nullBlock, maxNoOfNodes,
                             &printWarning, &unweightedPosStart, unweightedPos, &allIndexMaps))
         return;

    setFeatureIndex(pTree, fwMaxMotifLength, fwMaxPatternLength, &alphaInf, FALSE, TRUE, fwMotifs,
                    fwMotifLengths, annX.length > 0, NULL, intf.featMap, &intf.keyPool, &keyPoolSize,
                    &poolNextFree, annX.length > 0, TRUE);

    intf.markUsedOnly = FALSE;
    intf.markMotifsInSample = FALSE;
    intf.presence = FALSE;
    intf.zeroFeatures = TRUE;
    intf.pTree = pTree;
    intf.pErd = NULL;
    intf.pProf = &profiles;
    intf.alphaInf = &alphaInf;
    intf.motifLengths = motifLengths;
    intf.maxMotifLength = currMaxMotifLength;
    intf.svmIndex = svmIndex;
    intf.featureWeights = &featureWeights;
    intf.unweightedPosStart = &unweightedPosStart;
    intf.unweightedPos = *unweightedPos;
    intf.allIndexMaps = &allIndexMaps;

    if (annX.length > 0)
    {
        intf.annotationIndexMap = &annotationIndexMap;
        intf.maxNoOfNodes = maxNoOfNodes;
        intf.freeNode = &freeNode;
        intf.nullBlock = &nullBlock;
        intf.printWarning = &printWarning;
    }
    else
    {
        useHash = FALSE;
        dimFeatureSpace = motifs.length;
        intf.keyPoolSize = NULL;
        intf.keyPool = NULL;
        intf.poolNextFree = NULL;
        pKeypool = NULL;
        intf.annptr = NULL;
    }

    for (i = 0 ; i < numSamples; i++)
    {
        R_CheckUserInterrupt();

        iX = selX[i];
        intf.rowIndex = i;
        intf.seqptr = x.ptr[iX];
        intf.seqnchar = x.nchar[iX];

        if (annX.length > 0)
            intf.annptr = annX.ptr[iX];

        findMotifs(&intf);
        
        if (intf.kernelValue == -1)
        {
            initMatrixWithNA(profiles, profiles.nrow(), profiles.ncol());
            return;
        }

        if (normalized)
        {
            selCurr[0] = iX;
            currSeqLength = intf.seqnchar;

            getKMStdAnnMotif(kmOne, x, y, 1, 1, selCurr, selY, annCharset, annX, annY, motifs,
                             *motifLengths, nodeLimit, maxMotifLength, maxPatternLength, FALSE,
                             TRUE, presence, currSeqLength, &alphaInf);

            normValue = sqrt(kmOne(0,0));

            for (j = 0; j < profiles.ncol(); j++)
            {
                if (profiles(i,j) != 0)
                    profiles(i,j) /= normValue;
            }
        }
    }

    if (annX.length > 0)
        pKeypool = keyPool;

    return;
}

void freeHeapMotif()
{
    if (pFeatureHMap != NULL)
    {
        kh_destroy(fim, pFeatureHMap);
        pFeatureHMap = NULL;
    }

    if (pKeypool != NULL)
    {
        R_Free(pKeypool);
        pKeypool = NULL;
    }

    if (pKeys != NULL)
    {
        R_Free(pKeys);
        pKeys = NULL;
    }

    if (pUnweightedPos != NULL)
    {
        R_Free(pUnweightedPos);
        pUnweightedPos = NULL;
    }
}
