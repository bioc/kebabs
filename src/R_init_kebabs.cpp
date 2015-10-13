#include <R_ext/Rdynload.h>

// with this define only the R interfaces are included
#define __R_Interfaces_Only__

#include "ExplicitRepC.h"
#include "GappyPairC.h"
#include "MismatchC.h"
#include "MotifC.h"
#include "SpectrumC.h"
#include "SymmetricPairC.h"
#include "PredictionC.h"
#include "FeatureWeightsPosDepC.h"
#include "PredictionProfileC.h"
#include "Utils.h"


static const R_CallMethodDef callMethods[] = 
{
    /* ExplicitRepC.cpp */
    {".genExplRepC", (DL_FUNC) & genExplRepC, 25},
    {".freeHeapCallocsC", (DL_FUNC) & freeHeapCallocsC, 1},
    /* GappyPairC.cpp */
    {".gappyPairKernelMatrixC", (DL_FUNC) & gappyPairKernelMatrixC, 23},
    /* MismatchC.cpp */
    {".mismatchKernelMatrixC", (DL_FUNC) & mismatchKernelMatrixC, 16},
    /* MotifC.cpp */
    {".motifKernelMatrixC", (DL_FUNC) & motifKernelMatrixC, 26},
    {".validateMotifsC", (DL_FUNC) & validateMotifsC, 2},
    /* SpectrumC.cpp */
    {".spectrumKernelMatrixC", (DL_FUNC) & spectrumKernelMatrixC, 23},
    /* SymmetricPairC.cpp */
    {".symmetricPairKernelC", (DL_FUNC) & symmetricPairKernelC, 7},
    /* PredictionC.cpp */
    {".getPosDepPredOrProfC", (DL_FUNC) & getPosDepPredOrProfC, 29},
    /* FeatureWeightsPosDepC.cpp */
    {".getFeatureWeightsPosDepC", (DL_FUNC) & getFeatureWeightsPosDepC, 25},
    {".freeHeapFeatureWeightsC", (DL_FUNC) & freeHeapFeatureWeightsC, 0},
    /* PredictionProfileC.cpp */
    {".generatePredictionProfilesC", (DL_FUNC) & 
                                     generatePredictionProfilesC, 27},
    /* Utils.cpp */
    {".linearKerneldgRMatrixC", (DL_FUNC) & linearKerneldgRMatrixC, 11},
    {".linearKernelSparseKMdgRMatrixC", (DL_FUNC) & linearKernelSparseKMdgRMatrixC, 15},
    {".freeHeapLinearKernelC", (DL_FUNC) & freeHeapLinearKernelC, 0},
    {".matrixdgRMatrixProductC", (DL_FUNC) & matrixdgRMatrixProductC, 8},
    {".dgRMatrixNumericVectorProductC", (DL_FUNC) & 
                                        dgRMatrixNumericVectorProductC, 7},
    {NULL, NULL, 0}
};

extern "C" 
{

    void R_init_kebabs(DllInfo *info)
    {
        /* Register routines,allocate resources. */
        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    }

    void R_unload_kebabs(DllInfo *info)
    {
        /* Release resources. */
    }
}
