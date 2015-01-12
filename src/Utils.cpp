/*
 *   Utilitiy Routines
 *
 *   Source: Utils.cpp
 *   Author: Johannes Palme
 *
 *   Copyright (c) 2014, Johannes Palme
 *
 */

//#include <Rcpp.h>
#include "Kebabs.h"
#include "KernelUtils.h"

using namespace Rcpp;

static int *ptrP;
static int *ptrI;
static double *ptrX;


RcppExport SEXP createNAMatrix(int sizeX, int sizeY)
{
    int i, j;
    NumericMatrix km(sizeX, sizeY);

    if (sizeX > 0 && sizeY > 0)
    {
        for (i = 0; i < sizeX; i++)
        {
            R_CheckUserInterrupt();

            km(i, i) = NA_REAL;

            for (j = i + 1; j < sizeY; j++)
            {
                km(i, j) = NA_REAL;
                km(j, i) = NA_REAL;
            }
        }
    }

    return(km);
}

// matrix product between matrix and dgRMatrix - not implemented in Matrix
RcppExport SEXP matrixdgRMatrixProductC(SEXP mR, SEXP noRowXR, SEXP noColXR,
                                        SEXP noRowYR, SEXP noColYR,
                                        SEXP pR, SEXP jR, SEXP xR)
{
    int i, k, l;
    const void *vmax;
    int nrowx = as<int>(noRowXR);
    int ncolx = as<int>(noColXR);
    int nrowy = as<int>(noRowYR);
    int ncoly = as<int>(noColYR);

    if (ncolx != nrowy)
    {
        NumericMatrix res(0,0);
        return(res);
    }

    NumericMatrix res(nrowx, ncoly);
    vmax = vmaxget();
    NumericMatrix m(mR);
    IntegerVector p(pR);
    IntegerVector j(jR);
    NumericVector x(xR);

    for (l = 0; l < nrowy; l++)
    {
        R_CheckUserInterrupt();

        for (k = p[l]; k < p[l+1]; k++)
            for (i = 0; i < nrowx; i++)
                res(i, j[k]) += x[k] * m(i, l);
    }

    vmaxset(vmax);

    return(res);
}

// matrix product between dgRMatrix and vector - not implemented in Matrix
RcppExport SEXP dgRMatrixNumericVectorProductC(SEXP pR, SEXP jR, SEXP xR,
                                               SEXP noRowXR, SEXP noColXR,
                                               SEXP yR, SEXP lengthYR)
{
    int i, k;
    const void *vmax;
    int nrowx = as<int>(noRowXR);
    int ncolx = as<int>(noColXR);
    int lengthy = as<int>(lengthYR);

    if (ncolx != lengthy)
    {
        NumericMatrix res(0,0);
        return(res);
    }

    NumericVector res(nrowx);
    vmax = vmaxget();
    NumericMatrix y(yR);
    IntegerVector p(pR);
    IntegerVector j(jR);
    NumericVector x(xR);

    for (i = 0; i < nrowx; i++)
    {
        for (k = p[i]; k < p[i+1]; k++)
            res(i) += x[k] * y[j[k]];
    }
    
    vmaxset(vmax);
    
    return(res);
}

// generate sparse kernel matrix with similarities larger than lower limit
RcppExport SEXP linearKernelSparseKMdgRMatrixC(SEXP sizeXR, SEXP pXR, SEXP jXR, SEXP xXR, SEXP selxR,
                                               SEXP sizeYR, SEXP pYR, SEXP jYR, SEXP xYR, SEXP selyR,
                                               SEXP rowNamesR, SEXP colNamesR, SEXP symmetricR,
                                               SEXP diagR, SEXP lowerLimitR)
{
    int i, j, i1, j1, ind1, ind2, nextFree, numProtect, sizeX, sizeY, offset;
    int *iptr, *pptr;
    double * xptr, growBy;
    SEXP spm, dims, dimnames, slot_i, slot_p, slot_x;
    double kv;
    uint64_t vectorSize;
    const void *vmax;

    vmax = vmaxget();

    IntegerVector pX(pXR);
    IntegerVector jX(jXR);
    NumericVector xX(xXR);
    IntegerVector selX(selxR);
    IntegerVector selY(selyR);
    CharacterVector rowNames(rowNamesR);
    CharacterVector colNames(colNamesR);
    bool symmetric = as<bool>(symmetricR);
    bool withDiagonal = as<bool>(diagR);
    double lowerLimit = as<double>(lowerLimitR);
    
    sizeX = selX.size();
    
    if (symmetric)
        sizeY = sizeX;
    else
        sizeY = selY.size();
    
    if (withDiagonal)
        offset = 0;
    else
        offset = 1;
    
    // set initial size and growth factor for i and x arrays
    vectorSize = sizeX + sizeY;
    growBy = 1.6;
    
    pptr = (int *) Calloc(sizeY + 1, int);
    iptr = (int *) Calloc(vectorSize, int);
    xptr = (double *) Calloc(vectorSize, double);
    ptrP = pptr;
    ptrI = iptr;
    ptrX = xptr;
    
    // do processing
    nextFree = 0;

    if (symmetric)
    {
        for (j=0; j < sizeY; j++)
        {
            R_CheckUserInterrupt();

            j1 = selX[j];
            pptr[j] = nextFree;

            for (i=j+offset; i < sizeX; i++)
            {
                i1 = selX[i];
                ind1 = pX[i1];
                ind2 = pX[j1];
                kv = 0;
                
                while (ind1 < pX[i1+1]  && ind2 < pX[j1+1])
                {
                    if (jX[ind1] < jX[ind2])
                        ind1++;
                    else if (jX[ind1] > jX[ind2])
                        ind2++;
                    else
                    {
                        kv += xX[ind1] * xX[ind2];
                        ind1++;
                        ind2++;
                    }
                }
                
                if (kv > lowerLimit)
                {
                    // new element - check if enough space
                    if (nextFree + 1 > vectorSize)
                    {
                        // realloc arrays i and x
                        vectorSize *= growBy;
                        iptr = (int *) Realloc(iptr, vectorSize, int);
                        xptr = (double *) Realloc(xptr, vectorSize, double);
                        ptrI = iptr;
                        ptrX = xptr;
                    }
                    
                    iptr[nextFree] = i;
                    xptr[nextFree++] = kv;
                }
            }
        }

        pptr[sizeY] = nextFree;
    }
    else
    {
        IntegerVector pY(pYR);
        IntegerVector jY(jYR);
        NumericVector xY(xYR);

        for (j=0; j < sizeY; j++)
        {
            R_CheckUserInterrupt();

            j1 = selY[j];
            pptr[j] = nextFree;

            for (i=0; i < sizeX; i++)
            {
                i1 = selX[i];
                ind1 = pX[i1];
                ind2 = pY[j1];
                kv = 0;
                
                while (ind1 < pX[i1+1]  && ind2 < pY[j1+1])
                {
                    if (jX[ind1] < jY[ind2])
                        ind1++;
                    else if (jX[ind1] > jY[ind2])
                        ind2++;
                    else
                    {
                        kv += xX[ind1] * xY[ind2];
                        ind1++;
                        ind2++;
                    }
                }
                
                if (kv > lowerLimit)
                {
                    // new element - check if enough space
                    if (nextFree + 1 > vectorSize)
                    {
                        // realloc arrays i and x
                        vectorSize = (uint64_t) (vectorSize * growBy);
                        iptr = (int *) Realloc(iptr, vectorSize, int);
                        xptr = (double *) Realloc(xptr, vectorSize, double);
                        ptrI = iptr;
                        ptrX = xptr;
                    }
                    
                    iptr[nextFree] = i;
                    xptr[nextFree++] = kv;
                }
            }
        }

        pptr[sizeY] = nextFree;
    }

    vmaxset(vmax);

    // shrink to actual size
    iptr = (int *) Realloc(iptr, nextFree, int);
    xptr = (double *) Realloc(xptr, nextFree, double);
    
    // allocate sparse km as dgCMatrix
    numProtect = 0;
    spm = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    dims = PROTECT(Rf_allocVector(INTSXP, 2));
    SET_SLOT(spm, Rf_mkChar("Dim"), dims);
    INTEGER(dims)[0] = sizeX;
    INTEGER(dims)[1] = sizeY;
    numProtect = 2;
    
    dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_SLOT(spm, Rf_mkChar("Dimnames"), dimnames);
    numProtect++;

    if (rowNames.size() > 0)
        SET_VECTOR_ELT(dimnames, 0, rowNames);
    else
        SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    
    if (colNames.size() > 0)
        SET_VECTOR_ELT(dimnames, 1, colNames);
    else
        SET_VECTOR_ELT(dimnames, 1, R_NilValue);
    
    slot_p = PROTECT(Rf_allocVector(INTSXP, sizeY+1));
    SET_SLOT(spm, Rf_mkChar("p"), slot_p);
    
    for (i=0; i < sizeY + 1; i++)
        INTEGER(slot_p)[i] = pptr[i];
    
    free(pptr);
    ptrP = NULL;
    
    slot_i = PROTECT(Rf_allocVector(INTSXP, nextFree));
    SET_SLOT(spm, Rf_mkChar("i"), slot_i);
    
    for (i=0; i < nextFree; i++)
        INTEGER(slot_i)[i] = iptr[i];
    
    free(iptr);
    ptrI = NULL;
    
    slot_x = PROTECT(Rf_allocVector(REALSXP, nextFree));
    SET_SLOT(spm, Rf_mkChar("x"), slot_x);
    
    for (i=0; i < nextFree; i++)
        REAL(slot_x)[i] = xptr[i];
    
    free(xptr);
    ptrX = NULL;
    
    numProtect += 3;
    
    UNPROTECT(numProtect);
    
    return(spm);
}

RcppExport SEXP linearKerneldgRMatrixC(SEXP sizeXR, SEXP pXR, SEXP jXR, SEXP xXR, SEXP selxR,
                                       SEXP sizeYR, SEXP pYR, SEXP jYR, SEXP xYR, SEXP selyR, SEXP symmetricR)
{
    int i1, i2, ind1, ind2, ix, iy, sizeX, sizeY;
    double kv;
    IntegerVector pX(pXR);
    IntegerVector jX(jXR);
    NumericVector xX(xXR);
    IntegerVector selX(selxR);
    IntegerVector selY(selyR);
    bool symmetric = as<bool>(symmetricR);

    sizeX = selX.size();
    
    if (symmetric)
        sizeY = sizeX;
    else
        sizeY = selY.size();

    NumericMatrix km(sizeX, sizeY);

    if(symmetric)
    {
        for (ind1 = 0; ind1 < sizeX; ind1++)
        {
            R_CheckUserInterrupt();

            i1 = selX[ind1];

            for (ind2 = ind1; ind2 < sizeX; ind2++)
            {
                i2 = selX[ind2];

                kv = 0;
                ix = pX[i1];
                iy = pX[i2];

                // sparse dot product between two rows
                while (ix < pX[i1+1] && iy < pX[i2+1])
                {
                    if (jX[ix] < jX[iy])
                        ix++;
                    else if (jX[ix] > jX[iy])
                        iy++;
                    else
                    {
                        kv += xX[ix] * xX[iy];
                        ix++;
                        iy++;
                    }
                }

                km(ind1, ind2) = kv;
                km(ind2, ind1) = kv;
            }
        }
    }
    else
    {
        IntegerVector pY(pYR);
        IntegerVector jY(jYR);
        NumericVector xY(xYR);

        for (ind1 = 0; ind1 < sizeX; ind1++)
        {
            R_CheckUserInterrupt();

            i1 = selX[ind1];

            for (ind2 = 0; ind2 < sizeY; ind2++)
            {
                i2 = selY[ind2];

                kv = 0;
                ix = pX[i1];
                iy = pY[i2];

                // sparse dot product between two rows
                while (ix < pX[i1+1] && iy < pY[i2+1])
                {
                    if (jX[ix] < jY[iy])
                        ix++;
                    else if (jX[ix] > jY[iy])
                        iy++;
                    else
                    {
                        kv += xX[ix] * xY[iy];
                        ix++;
                        iy++;
                    }
                }

                km(ind1, ind2) = kv;
            }
        }
    }

    return(km);
}

extern "C" {

void freeHeapLinearKernelC()
{
    if (ptrP != NULL)
    {
        Free(ptrP);
        ptrP = NULL;
    }

    if (ptrI != NULL)
    {
        Free(ptrI);
        ptrI = NULL;
    }

    if (ptrX != NULL)
    {
        Free(ptrX);
        ptrX = NULL;
    }
}

}
