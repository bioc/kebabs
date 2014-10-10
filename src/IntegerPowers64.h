#ifndef __IntegerPowers64_C_H__

#define __IntegerPowers64_C_H__


// integer powers of uint64_t

static uint64_t ipow64(uint64_t base, uint8_t exponent)
{
    uint64_t result = 1;

    while (exponent)
    {
        if (exponent & 1)
            result *= base;
        exponent >>= 1;
        base *= base;
    }

    return result;
}

#endif
