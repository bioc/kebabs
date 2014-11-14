//
// C Routines for Sparse Matrix Random Access
//
// Source : SparseMatrixHash.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014   J o h a n n e s  P a l m e
//

//#include "Kebabs.h"
#include "KernelUtils.h"

extern "C"
{
    #include "khash.h"
}

using namespace Rcpp;

#if __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

KHASH_MAP_INIT_INT(access32, uint32_t)
KHASH_MAP_INIT_INT64(access64, uint32_t)

static khash_t(access32) *pAccessHashMap32;
static khash_t(access64) *pAccessHashMap64;

#if __clang__
#pragma clang diagnostic pop
#endif

// only supports dgCMatrix from the Matrix package for bulk random access to elements
void * generateAccessHashmap(int nrows, int ncols, uint64_t *rowIndices, int dimFeatureSpace, int minPos,
                             int *i, int *p, double *x, bool bits64)
{
    int j, result;
    uint32_t l, key32;
    uint64_t key64;
    khiter_t iter;
    khash_t(access32) *map32;
    khash_t(access64) *map64;

    if (bits64)
    {
        map64 = kh_init(access64);
        pAccessHashMap64 = map64;

        for (j = 0; j < ncols; j++)
        {
            for (l = p[j]; l < (uint32_t) p[j+1]; l++)
            {
                key64 = (minPos + j) * dimFeatureSpace + rowIndices[i[l]];
                iter = kh_put(access64, map64, key64, &result);

                if (result != -1)
                    kh_value(map64, iter) = l;
            }
        }

        return((void *) map64);
    }
    else
    {
        map32 = kh_init(access32);
        pAccessHashMap32 = map32;


        for (j = 0; j < ncols; j++)
        {
            for (l = p[j]; l < (uint32_t) p[j+1]; l++)
            {
                key32 = (minPos + j) * dimFeatureSpace + rowIndices[i[l]];
                iter = kh_put(access32, map32, key32, &result);

                if (result != -1)
                    kh_value(map32, iter) = l;
            }
        }

        return((void *) map32);
    }
}


// $$$ TODO REMOVE
//RcppExport SEXP sumRandomAccessMatrixElements(SEXP rowIndR, SEXP colIndR, SEXP nrowsR, SEXP ncolsR,
//                                              SEXP iR, SEXP pR, SEXP xR)
//{
//    int j, nrows, ncols;
//    double sum;
//    void *hashmap;
//    uint32_t key32;
//    uint64_t key64;
//    bool bits64;
//    khiter_t iter;
//    khash_t(access32) *map32;
//    khash_t(access64) *map64;
//    nrows = as<int>(nrowsR);
//    ncols = as<int>(ncolsR);
//    IntegerVector rowInd(rowIndR);
//    IntegerVector colInd(colIndR);
//    IntegerVector i(iR);
//    IntegerVector p(pR);
//    NumericVector x(xR);
//    SEXP result;
//
//    key64 = nrows * ncols;
//    bits64 = (key64 > (uint64_t) MAXUINT32);
//
//    hashmap = generateAccessHashmap(nrows, ncols, i, p, x, bits64);
//
//    sum = 0;
//
//    if (bits64)
//    {
//        map64 = (khash_t(access64) *) hashmap;
//
//        for (j = 0; j < rowInd.size(); j++)
//        {
//            key64 = (uint64_t) ((colInd[j] - 1) * nrows + rowInd[j] - 1);
//            iter = kh_get(access64, map64, key64);
//
//            if (iter != kh_end(map64))
//                sum += x[kh_value(map64, iter)];
//        }
//    }
//    else
//    {
//        map32 = (khash_t(access32) *) hashmap;
//
//        for (j = 0; j < rowInd.size(); j++)
//        {
//            key32 = (uint32_t) ((colInd[j] - 1) * nrows + rowInd[j] - 1);
//            iter = kh_get(access32, map32, key32);
//
//            if (iter != kh_end(map32))
//                sum += x[kh_value(map32, iter)];
//        }
//    }
//
//    PROTECT(result = Rf_allocVector(REALSXP, 1));
//    REAL(result)[0] = sum;
//    UNPROTECT(1);
//    return(result);
//}

//(void *) convertSparseMatrixToHashMap()
//{
//
//}
//
//(void *) convertHashMapToSparseMatrix()
//{
//}
