#ifndef __Utils_CPP_H__

#define __Utils_CPP_H__


# ifdef __R_Interfaces_Only__

RcppExport SEXP matrixdgRMatrixProductC(SEXP mR, SEXP noRowXR, SEXP noColXR,
                                        SEXP noRowYR, SEXP noColYR,
                                        SEXP pR, SEXP jR, SEXP xR);

RcppExport SEXP dgRMatrixNumericVectorProductC(SEXP pR, SEXP jR, SEXP xR,
                                               SEXP noRowXR, SEXP noColXR,
                                               SEXP yR, SEXP lengthYR);

RcppExport SEXP linearKerneldgRMatrixC(SEXP sizeXR, SEXP pXR, SEXP jXR, 
                                       SEXP xXR, SEXP selxR, SEXP sizeYR, 
                                       SEXP pYR, SEXP jYR, SEXP xYR, 
                                       SEXP selyR, SEXP symmetricR);

#else

#include <Rdefines.h>
#include <stdint.h>


RcppExport SEXP createNAMatrix(int sizeX, int sizeY);

#endif

#endif

