#include "urand.h"

#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <omp.h>

#define SEED_SIZE   128


static uint64_t seed[SEED_SIZE];


void urand_init()
{
    unsigned k;
    uint64_t start;

    start = time(NULL);

    for(k = 0; k < SEED_SIZE; ++k)
        seed[k] = start + k * UINT32_MAX;
}

float urand()
{
    return erand48((unsigned short*)&seed[omp_get_thread_num()]);
}