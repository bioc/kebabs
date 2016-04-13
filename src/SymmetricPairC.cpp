//
// C Routines for Symmetric Pair Kernel
//
// Source : SymmetricPairC.cpp
// Package: kebabs
// Author : J. P.
//
// Copyright (C) 2014-2016   J o h a n n e s  P a l m e
//

#include "Kebabs.h"
#include "KernelUtils.h"

using namespace Rcpp;

RcppExport SEXP symmetricPairKernelC(SEXP kmSIR, SEXP selXR, SEXP selYR,
                                     SEXP sizeXR, SEXP sizeYR, SEXP kernelTypeR,
                                     SEXP symmetricR)
{
    const void *vmax;
    int i, j, start;
    int sizeX = as<int>(sizeXR);
    int sizeY = as<int>(sizeYR);
    int kernelType = as<int>(kernelTypeR);
    bool symmetric = as<bool>(symmetricR);

    NumericMatrix km(sizeX, sizeY);

    vmax = vmaxget();

    // single instance kernel matrix
    NumericMatrix kmSI(kmSIR);
    NumericMatrix selx(selXR);
    NumericMatrix sely(selYR);

    start = 0;

    for (i = 0; i < sizeX; i++)
    {
        if (symmetric)
            start = i;

        switch (kernelType)
        {
            case KBS_ARITHMETIC_MEAN:
            {
                for (j = start; j < sizeY; j++)
                {
                    km(i, j) = (kmSI(selx(i,0), sely(j,0)) +
                                kmSI(selx(i,0), sely(j,1)) +
                                kmSI(selx(i,1), sely(j,0)) +
                                kmSI(selx(i,1), sely(j,1))) / 4;

                    if (symmetric)
                        km(j,i) = km(i,j);
                }

                break;
            }

            case KBS_TPPK: // tensor product pairwise kernel
            {
                for (j = start; j < sizeY; j++)
                {
                    km(i, j) = kmSI(selx(i,0), sely(j,0)) *
                               kmSI(selx(i,1), sely(j,1)) +
                               kmSI(selx(i,0), sely(j,1)) *
                               kmSI(selx(i,1), sely(j,0));

                    if (symmetric)
                        km(j,i) = km(i,j);
                }

                break;
            }

            default:
                break;
        }
    }

    vmaxset(vmax);

    return(km);
}
