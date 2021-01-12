#include "test_ntl.hpp"
#include <NTL/ZZ.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#define LINEARITY_CHECK
#include "tinymt32.h"
#include "f2p_gmp.h"
#include <string>
#include <stdio.h>
#include <inttypes.h>
#include <ctype.h>

using namespace NTL;
using namespace std;

// parameter no check
void tinymt32_add(tinymt32_t *a, tinymt32_t *b)
{
    for (int i = 0; i < 4; i++) {
        a->status[i] ^= b->status[i];
    }
}

void ok_print(bool ok, bool verbose)
{
    if (verbose) {
        if (ok) {
            printf("o");
        } else {
            printf("x");
        }
    }
}

void ntl_calc_jump(GF2X& jump, const GF2X& minpoly, ZZ& step)
{
    PowerXMod(jump, step, minpoly);
}

void gmp_calc_jump(GF2X& jump, mpz_t minpoly, mpz_t step)
{
    mpz_t jpoly;
    static char buff[2000];
    mpz_init(jpoly);
    f2p_calc_jump(jpoly, minpoly, step);
    f2p_get_hexstr(buff, jpoly);
    hexto_poly(jump, buff);
    mpz_clear(jpoly);
}

bool test_calc_jump(long degree, bool verbose)
{
    bool ok = true;
    GF2X poly;
    GF2X ntl_minpoly;
    BuildIrred(poly, degree);
    long step = degree * 10 + 7;
    ZZ ntl_step(step);
    GF2X ntl_jump;
    GF2X f2p_jump;
    mpz_t f2p_minpoly;
    mpz_t f2p_step;
    mpz_inits(f2p_minpoly, f2p_step, NULL);
    mpz_set_ui(f2p_step, step);
    string work;
    for (int i = 0; i < 10; i++) {
        BuildRandomIrred(ntl_minpoly, poly);
        to_hexstring(work, ntl_minpoly);
        f2p_set_hexstr(f2p_minpoly, work.c_str());
        ntl_calc_jump(ntl_jump, ntl_minpoly, ntl_step);
        gmp_calc_jump(f2p_jump, f2p_minpoly, f2p_step);
        if (ntl_jump != f2p_jump) {
            to_hexstring(work, ntl_jump);
            printf("ntl_jump:%s\n", work.c_str());
            to_hexstring(work, f2p_jump);
            printf("f2p_jump:%s\n", work.c_str());
            ok = false;
        }
        ok_print(ok, verbose);
    }
    if (verbose) {
        printf("\n");
    }
    mpz_clears(f2p_minpoly, f2p_step, NULL);
    return ok;
}


int main(int argc, char * argv[])
{
    bool verbose = false;
    bool ok = true;
    long degree = 127;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = true;
    }
    if (argc > 2) {
        degree = strtol(argv[2], NULL, 10);
    }
#if 0
    tiny32.mat1 = 0x8f7011ee;
    tiny32.mat2 = 0xfc78ff1f;
    tiny32.tmat = 0x3793fdff;

    tinymt32_init(&tiny32, 1);
#endif

    ok = ok && test_calc_jump(degree, verbose);
    if (ok) {
        return 0;
    } else {
        return -1;
    }
}
