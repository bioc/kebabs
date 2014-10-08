/*
 *   Utilitiy Routines
 *
 *   Source: Utils.cpp
 *   Author: Johannes Palme
 *
 *   Copyright (c) 2014, Johannes Palme
 *
 */

#include <Rcpp.h>

using namespace Rcpp;

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

RcppExport SEXP linearKerneldgRMatrixC(SEXP sizeXR, SEXP pXR, SEXP jXR, SEXP xXR, SEXP selxR,
                                       SEXP sizeYR, SEXP pYR, SEXP jYR, SEXP xYR, SEXP selyR, SEXP symmetricR)
{
    int i1, i2, ind1, ind2, ix, iy;
    double kv;
    bool subsetX = FALSE, subsetY =FALSE;
    IntegerVector pX(pXR);
    IntegerVector jX(jXR);
    NumericVector xX(xXR);
    IntegerVector *selX = NULL;
    IntegerVector *selY = NULL;

    int sizeX = as<int>(sizeXR);
    int sizeY = as<int>(sizeYR);
    bool symmetric = as<bool>(symmetricR);


    if (!Rf_isNull(selxR))
    {
        IntegerVector selX1(selxR);
        sizeX = selX1.size();
        selX = &selX1;
        subsetX = TRUE;
    }

    if (!Rf_isNull(selyR))
    {
        IntegerVector selY1(selyR);
        sizeY = selY1.size();
        selY = &selY1;
        subsetY = TRUE;
    }

    NumericMatrix km(sizeX, sizeY);

    if(symmetric)
    {
        for (ind1 = 0; ind1 < sizeX; ind1++)
        {
            R_CheckUserInterrupt();

            if (subsetX)
                i1 = (*selX)[ind1];
            else
                i1 = ind1;

            for (ind2 = ind1; ind2 < sizeX; ind2++)
            {
                if (subsetX)
                    i2 = (*selX)[ind2];
                else
                    i2 = ind2;

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

            if (subsetX)
                i1 = (*selX)[ind1];
            else
                i1 = ind1;

            for (ind2 = 0; ind2 < sizeY; ind2++)
            {
                if (subsetY)
                    i2 = (*selY)[ind2];
                else
                    i2 = ind2;

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
                km(ind2, ind1) = kv;
            }
        }
    }


    return(km);
}
