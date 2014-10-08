/*
 *   Utilitiy Routines for SequenceKernels
 *
 *   Source: KernelUtils.cpp
 *   Author: Johannes Palme
 *
 *   Copyright (c) 2014, Johannes Palme
 *
 */

#include "KernelUtils.h"

#define BC_DNA_EXACT        1
#define BC_DNA_IUPAC        2
#define BC_RNA_EXACT        3
#define BC_RNA_IUPAC        4
#define BC_AA_EXACT         5
#define BC_AA_IUPAC         6
#define BC_MAX_BIO_CHARSET  6


// define charactersets for biological sequences

//DNA EXACT
static char * DNA_EXACT = (char *) "ACGT";

//DNA IUPAC
static char * DNA_IUPAC = (char *) "ACGTMRWSYKVHDBN-+";

//RNA EXACT
static char * RNA_EXACT = (char *) "ACGU";

//RNA IUPAC
static char * RNA_IUPAC = (char *) "ACGUMRWSYKVHDBN-+";

//AA EXACT including Selenocystein - Sec - U, without Pyrrolysin O
static char * AA_EXACT = (char *) "ACDEFGHIKLMNPQRSTUVWY";

//AA IUPAC
//static char * AA_IUPAC = (char *) "ABCDEFGHIJKLMNPQRSTUVWXYZ";


using namespace Rcpp;

// Single translation table
// as union of all allowed Biostrings characters
// string type not needed!
// *+-ABCDEFGHIKLMNPQRSTUVWXY
const int allIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 1,-1, 2,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 3, 4, 5, 6, 7, 8, 9,10,11,-1,12,13,14,15,-1,
    16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

// same including lower case characters
const int allLowerIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 1,-1, 2,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 3, 4, 5, 6, 7, 8, 9,10,11,-1,12,13,14,15,-1,
    16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,-1,
    -1, 3, 4, 5, 6, 7, 8, 9,10,11,-1,12,13,14,15,-1,
    16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

//pure DNA 4 Letter without lower case
const int dnaIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

//pure DNA 4 Letter with lower case
const int dnaLowerIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

//pure RNA 4 Letter without lower case
const int rnaIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

//pure RNA 4 Letter with lower case
const int rnaLowerIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

// pure AA 21 Letter without lower case
// including Selenocystein U (or Sec)
const int aaIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1, 2, 3, 4, 5, 6, 7,-1, 8, 9,10,11,-1,
    12,13,14,15,16,17,18,19,-1,20,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

// pure AA 21 Letter with lower case
// including Selenocystein U (or Sec)
const int aaLowerIndexMap[MAX_CHAR]
={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1, 2, 3, 4, 5, 6, 7,-1, 8, 9,10,11,-1,
    12,13,14,15,16,17,18,19,-1,20,-1,-1,-1,-1,-1,-1,
    -1, 0,-1, 1, 2, 3, 4, 5, 6, 7,-1, 8, 9,10,11,-1,
    12,13,14,15,16,17,18,19,-1,20,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

static int reverseIndexMap[MAX_CHAR];
static int indicesUnmapped[MAX_CHAR];
static int indicesReverseUnmapped[MAX_CHAR];



void getAllIndexMaps(struct allIndMaps *indexMaps)
{
    indexMaps->dna = (int *) &dnaIndexMap;
    indexMaps->dnaLower = (int *) &dnaLowerIndexMap;
    indexMaps->rna = (int *) &rnaIndexMap;
    indexMaps->rnaLower = (int *) &rnaLowerIndexMap;
    indexMaps->aa = (int *) &aaIndexMap;
    indexMaps->aaLower = (int *) &aaLowerIndexMap;
    indexMaps->all = (int *) &allIndexMap;
    indexMaps->allLower = (int *) &allLowerIndexMap;
    indexMaps->reverse = (int *) &reverseIndexMap;
    indexMaps->unmapped = (int *) &indicesUnmapped;
    indexMaps->reverseUnmapped = (int *) &indicesReverseUnmapped;
    return;
}

uint64_t getDimFeatureSpace(int kernelType, int k, int m, int numAlphabetChars, int numAnnChars,
                                   int numMotifs, int maxMotifLength)
{
    switch (kernelType)
    {
        case SPECTRUM:
        case MISMATCH:
            if (numAnnChars > 0)
                return(pow((uint64_t) numAlphabetChars, k) * pow((uint64_t) numAnnChars, k));
            else
                return(pow((uint64_t) numAlphabetChars, k));

        case GAPPY_PAIR:
            if (numAnnChars > 0)
            {
                return((m + 1) * pow((uint64_t) numAlphabetChars, 2 * k) *
                       pow((uint64_t) numAnnChars, 2 * k));
            }
            else
                return((m + 1) * pow((uint64_t) numAlphabetChars, 2 * k));

        case MOTIF:
            if (numAnnChars > 0)
                return(numMotifs * pow((uint64_t) numAnnChars, maxMotifLength));
            else
                return(numMotifs);

        default:
            return(0);
    }
}

void getAlphabetInfo(int seqType, bool lowercase, bool unmapped,
                      struct alphaInfo *info, struct allIndMaps *allIndexMaps)
{
    getAllIndexMaps(allIndexMaps);

    int i, code;

    if (info != NULL)
    {
        info->seqType = seqType;
        info->lowercase = lowercase;
        info->unmapped = unmapped;

        switch (seqType)
        {
            case BC_DNA_EXACT:
                if (lowercase)
                    info->indexMap = allIndexMaps->dnaLower;
                else
                    info->indexMap = allIndexMaps->dna;

                info->currAlphabet = DNA_EXACT;

                for (i = strlen(info->currAlphabet); i < MAX_CHAR; i++)
                    allIndexMaps->reverseUnmapped[i] = -1;

                for (i = 0; i < (int) strlen(info->currAlphabet); i++)
                {
                    if (unmapped)
                    {
                        code = DNAorRNAencode(DNA_EXACT[i], DT_DNA_BIOSTRING);

                        if (code >= 0)
                        {
                            allIndexMaps->unmapped[i] = code;
                            allIndexMaps->reverseUnmapped[code] = i;
                        }
                        else
                            Rprintf("Error in reverse mapping of code for char %c",
                                    DNA_EXACT[i]);
                    }
                    else
                        allIndexMaps->unmapped[i] = i;
                }

                if (unmapped)
                    info->seqIndexMap = allIndexMaps->reverseUnmapped;
                else
                {
                    if (lowercase)
                        info->seqIndexMap = allIndexMaps->dnaLower;
                    else
                        info->seqIndexMap = allIndexMaps->dna;
                }

                break;

            case BC_DNA_IUPAC:
                if (lowercase)
                    info->indexMap = allIndexMaps->dnaLower;
                else
                    info->indexMap = allIndexMaps->dna;

                info->currAlphabet = DNA_IUPAC;

                if (unmapped)
                    info->seqIndexMap = allIndexMaps->reverseUnmapped;
                else
                {
                    if (lowercase)
                        info->seqIndexMap = allIndexMaps->dnaLower;
                    else
                        info->seqIndexMap = allIndexMaps->dna;
                }

                break;

            case BC_RNA_EXACT:
                if (lowercase)
                    info->indexMap = allIndexMaps->rnaLower;
                else
                    info->indexMap = allIndexMaps->rna;

                info->currAlphabet = RNA_EXACT;

                for (i = strlen(info->currAlphabet); i < MAX_CHAR; i++)
                    allIndexMaps->reverseUnmapped[i] = -1;

                for (i = 0; i < (int) strlen(info->currAlphabet); i++)
                {
                    if (unmapped)
                    {
                        code = DNAorRNAencode(RNA_EXACT[i], DT_RNA_BIOSTRING);

                        if (code >= 0)
                        {
                            allIndexMaps->unmapped[i] = code;
                            allIndexMaps->reverseUnmapped[code] = i;
                        }
                        else
                            Rprintf("Error in reverse mapping of code for char %c",
                                    RNA_EXACT[i]);
                    }
                    else
                        allIndexMaps->unmapped[i] = i;
                }

                if (unmapped)
                    info->seqIndexMap = allIndexMaps->reverseUnmapped;
                else
                {
                    if (lowercase)
                        info->seqIndexMap = allIndexMaps->rnaLower;
                    else
                        info->seqIndexMap = allIndexMaps->rna;
                }

                break;

            case BC_RNA_IUPAC:
                if (lowercase)
                    info->indexMap = allIndexMaps->rnaLower;
                else
                    info->indexMap = allIndexMaps->rna;

                info->currAlphabet = RNA_IUPAC;

                if (unmapped)
                    info->seqIndexMap = allIndexMaps->reverseUnmapped;
                else
                {
                    if (lowercase)
                        info->seqIndexMap = allIndexMaps->rnaLower;
                    else
                        info->seqIndexMap = allIndexMaps->rna;
                }

                break;

            case BC_AA_EXACT:
                if (lowercase)
                    info->indexMap = allIndexMaps->aaLower;
                else
                    info->indexMap = allIndexMaps->aa;

                info->currAlphabet = AA_EXACT;

                for (i = 0; i < (int) strlen(info->currAlphabet); i++)
                    allIndexMaps->unmapped[i] = i;

                if (lowercase)
                    info->seqIndexMap = allIndexMaps->aaLower;
                else
                    info->seqIndexMap = allIndexMaps->aa;

                break;

                // BC_AA_IUPAC does not make sense because it only contains
                // ambiguities for 2 specific aa pairs => treat AA always exact

            default:
                if (lowercase)
                    info->indexMap = allIndexMaps->allLower;
                else
                    info->indexMap = allIndexMaps->all;

                if (lowercase)
                    info->seqIndexMap = allIndexMaps->allLower;
                else
                    info->seqIndexMap = allIndexMaps->all;

                break;
        }

        // create reverse mapping
        for (i = 0; i < MAX_CHAR; i++)
            allIndexMaps->reverse[i] = -1;

        info->maxAlphaIndex = -1;
        info->numAlphabetChars = 0;

        // run reverse to overwrite lowercase with uppercase
        for (i=MAX_CHAR-1; i>=0; i--)
        {
            if (info->indexMap[i] >= 0)
            {
                if (info->indexMap[i] > info->maxAlphaIndex)
                    info->maxAlphaIndex = info->indexMap[i];

                allIndexMaps->reverse[info->indexMap[i]] = i;
            }
        }

        for (i = 0; i < MAX_CHAR; i++)
        {
            if (allIndexMaps->reverse[i] >= 0)
                info->numAlphabetChars++;
        }

        info->reverseIndexMap = allIndexMaps->reverse;

        // $$$ TODO
        // include handling for indices unmapped in case of XStringset
    }
}

void initAnnotationMaps(ByteStringVector annCharset, Rcpp::IntegerVector *annIndexMap,
                        Rcpp::IntegerVector *revAnnMap)
{
    int i;

    for (i = 0; i < MAX_CHAR; i++)
    {
        (*annIndexMap)[i] = -1;
        (*revAnnMap)[i] = -1;
    }

    for (i = 0; i < annCharset.nchar[0]; i++)
    {
        (*revAnnMap)[i] = annCharset.ptr[0][i];
        (*annIndexMap)[annCharset.ptr[0][i]] = i;
    }
}

void initMatrixWithNA(Rcpp::NumericMatrix m, int sizeX, int sizeY)
{
    int i, j;

    for (i = 0; i < sizeX; i++)
    {
        for (j = 0; j < sizeY; j++)
            m(i,j) = NA_REAL;
     }

    return;
}

RcppExport SEXP generateEmptyExplicitRep(int sizeX, bool sparse)
{
    int i;

    if (sparse)
    {
        SEXP ers = PROTECT(NEW_OBJECT(MAKE_CLASS("ExplicitRepresentationSparse")));

        SEXP dims = PROTECT(Rf_allocVector(INTSXP, 2));
        SET_SLOT(ers, Rf_mkChar("Dim"), dims);
        INTEGER(dims)[0] = sizeX;
        INTEGER(dims)[1] = 0;
        SEXP slot_p = PROTECT(Rf_allocVector(INTSXP, sizeX + 1));
        SET_SLOT(ers, Rf_mkChar("p"), slot_p);

        for (i = 0; i < (sizeX+1); i++)
            INTEGER(slot_p)[i] = 0;

        UNPROTECT(3);

        return(ers);
    }
    else
    {
        return(createNAMatrix(sizeX, 0));
    }

}

