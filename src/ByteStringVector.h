#ifndef __ByteStringVector_C_H__

#define __ByteStringVector_C_H__

#include <Rinternals.h>

#define DT_STANDARD        1
#define DT_DNA_BIOSTRING   2
#define DT_RNA_BIOSTRING   3

typedef struct {
    int length;
    int *nchar;
    char **ptr;
} ByteStringVector;

ByteStringVector charVector2ByteStringVec(SEXP cvR);
ByteStringVector XStringSet2ByteStringVec(SEXP xssR);
char DNAorRNAdecode(int code, int decodeType);
char DNAorRNAencode(int c, int encodeType);
//void getBioStringCodes(const char * xStringsetClass, int * codes, int * noOfCodes);

#endif
