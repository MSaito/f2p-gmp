#include "tinymt32.h"
#include "f2p-gmp.h"
#include <stdio.h>
#include <inttypes.h>

tinymt32_t tiny32;

unsigned int dummy()
{
    return tinymt32_generate_uint32(&tiny32);
}

int test_minpoly()
{
    char poly[200];
    minpoly(poly, dummy, 127);
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
