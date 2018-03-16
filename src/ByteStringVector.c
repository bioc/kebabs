#include "ByteStringVector.h"
#include "IRanges_interface.h"
#include "Biostrings_interface.h"

ByteStringVector charVector2ByteStringVec(SEXP cvR)
{
    ByteStringVector result;
    result.length = LENGTH(cvR);

    if (result.length > 0)
    {
        result.nchar = (int *)R_alloc(result.length, sizeof(int));
        result.ptr   = (char **)R_alloc(result.length, sizeof(char *));

        for (R_len_t i = 0; i < result.length; i++)
        {
            result.ptr[i] = (char *)CHAR(STRING_ELT(cvR, i));
            result.nchar[i] = strlen(result.ptr[i]);
        }
    }
    else
    {
        result.nchar = NULL;
        result.ptr = NULL;
    }

    return(result);
}

ByteStringVector XStringSet2ByteStringVec(SEXP xssR)
{
    ByteStringVector result;
    XStringSet_holder holder = hold_XStringSet(xssR);
    result.length = get_XStringSet_length(xssR);

    if (result.length > 0)
    {
        result.nchar = (int *)R_alloc(result.length, sizeof(int));
        result.ptr   = (char **)R_alloc(result.length, sizeof(char *));

        for (int i = 0; i < result.length; i++)
        {
            Chars_holder s = get_elt_from_XStringSet_holder(&holder, i);
            result.nchar[i] = s.length;
            result.ptr[i] = (char *)s.ptr;
        }
    }
    else
    {
        result.nchar = NULL;
        result.ptr = NULL;
    }

    return(result);
}

char DNAorRNAdecode(int c, int decodeType)
{
    if (decodeType == DT_RNA_BIOSTRING)
        return RNAdecode((char) c);
    else
        return DNAdecode((char) c);
}

char DNAorRNAencode(int c, int encodeType)
{
    if (encodeType == DT_RNA_BIOSTRING)
        return RNAencode((char) c);
    else
        return DNAencode((char) c);
}
