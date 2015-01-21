#ifndef __Kebabs_C_H__

#define __Kebabs_C_H__

#include <stdint.h>

#if defined(_WIN32) || defined(_WIN64)
#define MAX_BLOCK            8388608         // 2^23 nodes maximum
#define drand48()            (rand() * (1. / RAND_MAX))
#else
#define MAX_BLOCK            33554432         // 2^25 nodes maximum
#endif

#define FEATURE_SPACE_LIMIT  0x100000000000   // 100 billion
#define MAX_FEATURES         1073741824       // 2^30 = 1 billion
#define MAX_ALPHA_SIZE       32
#define MAX_CHAR             256
#define CHAR_UNMAPPED        FALSE
#define MAXUINT8             0xFF
#define MAXUINT16            0xFFFF
#define MAXINT32             0x7FFFFFFF
#define MAXUINT32            0xFFFFFFFF
#define MAXUINT64            0xFFFFFFFFFFFFFFFF
#define MAXINDEX32           MAXUINT32 - 1
#define MAXINDEX64           MAXUINT64 - 1
#define ANNOT_BLOCK_SIZE     32
#define START_ANNOTATION     97
#define VALID_FEATURE        0xFFFFFFFE
#define HASH_MAP_LIMIT       16777216     // not larger than 2^24*4 = 64MB


#define DENSE_ER             0
#define SPARSE_ER            1

// Supported kernels
#define SPECTRUM             1
#define MIXED_SPECTRUM       2
#define MISMATCH             3
#define MOTIF                4
#define WEIGHTED_DEGREE      5
#define GAPPY_PAIR           6

// Pairwise kernel types
#define KBS_ARITHMETIC_MEAN  1
#define KBS_TPPK             2  // tensor product pairwise kernel

#define KBS_UNSORTED         0
#define KBS_SORT_BY_FEATURE  1
#define KBS_SORT_BY_POSITION 2

// definitions for sparse matrix from Matrix package
#define MATRIX_I_SLOT        "i"
#define MATRIX_J_SLOT        "j"
#define MATRIX_P_SLOT        "p"
#define MATRIX_X_SLOT        "x"
#define MATRIX_DIM_SLOT      "Dim"
#define MATRIX_DIMNAMES_SLOT "Dimnames"



struct indexBlock {
    uint32_t idx[MAX_ALPHA_SIZE];
};

struct treeNode {
    struct indexBlock ib;
    uint32_t          value;
    uint8_t           leaf;
};

struct prefTree {
    struct treeNode node[MAX_BLOCK];
};

#endif
