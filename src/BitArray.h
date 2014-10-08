#ifndef __BitArray_C_H__

#define __BitArray_C_H__

#define BITS_PER_UINT64    64
#define POS_TO_SHIFT_64    6
#define BITS_PER_UINT8     8
#define POS_TO_SHIFT_8     3

#define NO_BITS_TO_UINT_64(noBits)                 \
((noBits + BITS_PER_UINT64) >> POS_TO_SHIFT_64)

#define NO_BITS_TO_UINT_8(noBits)                 \
((noBits + BITS_PER_UINT8) >> POS_TO_SHIFT_8)

#define BITS_IN_LAST_WORD(noBits)                  \
(noBits % BITS_PER_UINT64)

#define DEFINE_BITARRAY(name, noBits)                \
uint64_t *name = (uint64_t *) R_alloc(NO_BITS_TO_UINT_64(noBits), sizeof(uint64_t))

#define SET_BITARRAY(name, noBits, bit)                  \
memset(name, bit, NO_BITS_TO_UINT_64(noBits) * sizeof(uint64_t))

#define COPY_BITARRAY(dest, source, noBits)                  \
memcpy(dest, source, NO_BITS_TO_UINT_8(noBits))

void static setBit(uint64_t *bitarray, uint64_t pos)
{
    bitarray[NO_BITS_TO_UINT_64(pos) - 1] = (bitarray[NO_BITS_TO_UINT_64(pos) - 1] |
                                             (0x1ull << (pos % BITS_PER_UINT64)));
}

//void static resetBit(uint64_t *bitarray, uint64_t pos)
//{
//    bitarray[NO_BITS_TO_UINT_64(pos) - 1] = (bitarray[NO_BITS_TO_UINT_64(pos) - 1] &
//                                             ~(0x1ull << (pos % BITS_PER_UINT64)));
//}

int static getBit(uint64_t *bitarray, uint64_t pos)
{
    return((bitarray[NO_BITS_TO_UINT_64(pos) - 1] &
            (0x1ull << (pos % BITS_PER_UINT64))) >> (pos % BITS_PER_UINT64));
}

#endif
