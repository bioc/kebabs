#ifndef __SparseMatrixHash_C_H__

#define __SparseMatrixHash_C_H__

void * generateAccessHashmap(int nrows, int ncols, uint64_t *rowIndices, int dimFeatureSpace,
                             int minPos, int *i, int *p, double *x, bool bits64);

#endif
