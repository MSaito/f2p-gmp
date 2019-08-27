#include "tinymt32.h"
#include "f2p-gmp.h"
#include <stdio.h>
#include <inttypes.h>

tinymt32_t tiny32;

unsigned int dummy()
{
    return tinymt32_generate_uint32(&tiny32);
}

int f2p_annihilate(char * poly, f2rng gen, int mexp)
{
    int bit = 0;
    mpz_t pol;
    mpz_init(pol);
    f2p_set_hexstr(pol, poly);
    for (int i = 0; i <= mexp; i++) {
        if (mpz_tstbit(pol, i) == 1) {
            bit ^= gen() & 1;
        }
    }
    mpz_clear(pol);
    if (bit == 0) {
        return 1;
    } else {
        return 0;
    }
}

int test_minpoly()
{
    char poly[200];
    f2p_minpoly(poly, dummy, 127);
    printf("minpoly end\n");
    //int deg = f2p_degree(poly);
    //printf("deg(poly) = %d\n", deg);
    printf("poly = %s\n", poly);
    mpz_t p;
    mpz_init(p);
    f2p_set_hexstr(p, poly);
    int deg = f2p_degree(p);
    printf("deg(poly) = %d\n", deg);
    mpz_clear(p);
    int ok = 1;
    for (int i = 0; i < 100; i++) {
        int r = f2p_annihilate(poly, dummy, 127);
        ok &= r;
        dummy();
        if (r == 1) {
            printf("o");
        } else {
            printf("x");
        }
    }
    printf("\n");
    if (ok) {
        printf("annihilate OK\n");
    } else {
        printf("annihilate NG\n");
    }
    return 0;
}

int main(int argc, char * argv[])
{
    int verbose = 0;
    int r = 0;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = 1;
    }
    tiny32.mat1 = 0x8f7011ee;
    tiny32.mat2 = 0xfc78ff1f;
    tiny32.tmat = 0x3793fdff;

    tinymt32_init(&tiny32, 1);
    r += test_minpoly();
    if (r == 0) {
        return 0;
    } else {
        return -1;
    }
}
