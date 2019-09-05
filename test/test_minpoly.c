#include "debug_f2p.h"
#define LINEARITY_CHECK
#include "tinymt32.h"
#include "f2p_gmp.h"
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>

int f2p_annihilate(mpz_t minpoly, tinymt32_t * tiny32, int mexp)
{
    int bit = 0;
    for (int i = 0; i <= mexp; i++) {
        if (mpz_tstbit(minpoly, i) == 1) {
            //if (mpz_tstbit(pol, mexp - i) == 1) {
            bit ^= tinymt32_generate_uint32(tiny32) & 1;
        }
    }
    if (bit == 0) {
        return 1;
    } else {
        return 0;
    }
}

int test_minpoly(int seed, int verbose)
{
    if (verbose) {
        printf("start test_minpoly\n");
    }
    char buff[2000];
    mpz_t poly;
    mpz_t seq;
    mpz_inits(poly, seq, NULL);
    tinymt32_t tiny32;
    int mexp = 127;
    tiny32.mat1 = 0x8f7011ee;
    tiny32.mat2 = 0xfc78ff1f;
    tiny32.tmat = 0x3793fdff;
    tinymt32_init(&tiny32, seed);
    for (int i = 0; i < mexp * 2; i++) {
        uint32_t x = tinymt32_generate_uint32(&tiny32);
        if (x & 1) {
            mpz_setbit(seq, i);
        } else {
            mpz_clrbit(seq, i);
        }
    }
    f2p_get_hexstr(buff, seq);
    printf("seq = %s\n", buff);
    f2p_minpoly(poly, seq, mexp);
    int deg = f2p_degree(poly);
    //int deg = f2p_degree(poly);
    //printf("deg(poly) = %d\n", deg);
    f2p_get_hexstr(buff, poly);
    printf("poly = %d,%s\n", deg, buff);
    int ok = 1;
#if 0
    for (int i = 0; i < 100; i++) {
        int r = f2p_annihilate(poly, &tiny32, 127);
        ok &= r;
        tinymt32_generate_uint32(&tiny32);
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
#endif
    mpz_clears(poly, seq, NULL);
    if (verbose) {
        printf("end test_minpoly\n");
    }
    return ok;
}


int chk_irre_bin(const char * pol, int expected)
{
    mpz_t p;
    //char buff[2000];
    mpz_init(p);
    f2p_set_binstr(p, pol);
    int r = f2p_is_irreducible(p);
    mpz_clear(p);
    if (r != expected) {
        printf("chk_irre_bin %s expected %d but %d\n", pol, expected, r);
        return 0;
    } else {
        return 1;
    }
}

int test_irreducible(int verbose)
{
    if (verbose) {
        printf("start test_irreducible\n");
    }
    int ok = 1;
    ok &= chk_irre_bin("11", 1);
    ok &= chk_irre_bin("10", 1);
    ok &= chk_irre_bin("111", 1);
    ok &= chk_irre_bin("110", 0);
    ok &= chk_irre_bin("100", 0);
    ok &= chk_irre_bin("101", 0);
    ok &= chk_irre_bin("1111", 0);
    mpz_t pol;
    mpz_init(pol);
    f2p_set_hexstr(pol, "d8524022ed8dff4a8dcc50c798faba43");
    int r = f2p_is_irreducible(pol);
    if (!r) {
        printf("d8524022ed8dff4a8dcc50c798faba43 is not irreducible\n");
        ok = 0;
    }
    f2p_set_hexstr(pol, "d8524022ed8dff4a8dcc50c798faba47");
    r = f2p_is_irreducible(pol);
    if (r) {
        printf("d8524022ed8dff4a8dcc50c798faba47 is irreducible\n");
        ok = 0;
    }
    mpz_clear(pol);
    if (verbose) {
        printf("start test_irreducible\n");
    }
    return ok;
}

#if 0
void f2p_exeuclid3(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y, int m,
                   int wp, f2p_wm_t *wm)
{
    PUTS("f2p_exeuclid3 start\n");
    //mp_bitcnt_t max_deg = m;
    int zcmp = mpz_cmp_ui(x, 0);
    assert(zcmp != 0);
    zcmp = mpz_cmp_ui(y, 0);
    assert(zcmp != 0);
    mpz_t *q1 = &(wm->ar[wp++]);
    mpz_t *r0 = &(wm->ar[wp++]);
    mpz_t *r1 = &(wm->ar[wp++]);
    mpz_t *r2 = &(wm->ar[wp++]);
    mpz_t *a0 = &(wm->ar[wp++]);
    mpz_t *a1 = &(wm->ar[wp++]);
    mpz_t *a2 = &(wm->ar[wp++]);
    mpz_t *b0 = &(wm->ar[wp++]);
    mpz_t *b1 = &(wm->ar[wp++]);
    mpz_t *b2 = &(wm->ar[wp++]);
    mpz_t *tmp = &(wm->ar[wp++]);
    //mpz_t *t; 後で使う
    assert(wp <= wm->max_size);
    PRT("x= ", x);
    PRT("y= ", y);
    mpz_set(*r0, x);
    mpz_set(*r1, y);
    mpz_set_ui(*a0, 1);
    mpz_set_ui(*a1, 0);
    mpz_set_ui(*b0, 0);
    mpz_set_ui(*b1, 1);
    PRT("r0 = ", *r0);
    PRT("r1 = ", *r1);
    PRT("a0 = ", *a0);
    PRT("a1 = ", *a1);
    PRT("b0 = ", *b0);
    PRT("b1 = ", *b1);
    //int dr = f2p_degree(*r0);
    //while (dr >= m) {
    for (;;) {
        //while (mpz_cmp_ui(*r1, 0) > 0) {
        PUTS("f2p_divrem\n");
        f2p_divrem(*q1, *r2, *r0, *r1, wp, wm); // 2 wm
        PUTS("f2p_mul\n");
        f2p_mul(*tmp, *q1, *a1, wp, wm); // 2 wm
        PUTS("f2p_add\n");
        f2p_add(*a2, *a0, *tmp);
        PUTS("f2p_mul\n");
        f2p_mul(*tmp, *q1, *b1, wp, wm);// 2 wm
        PUTS("f2p_add\n");
        f2p_add(*b2, *b0, *tmp);
        PRT("q1 = ", *q1);
        PRT("tmp = ", *tmp);
        PRT("r0 = ", *r0);
        PRT("r1 = ", *r1);
        PRT("r2 = ", *r2);
        PRT("a0 = ", *a0);
        PRT("a1 = ", *a1);
        PRT("a2 = ", *a2);
        PRT("b0 = ", *b0);
        PRT("b1 = ", *b1);
        PRT("b2 = ", *b2);
        int dr = f2p_degree(*r2);
        if (dr < m) {
            break;
        }
        mpz_set(*r0, *r1);
        mpz_set(*r1, *r2);
        mpz_set(*a0, *a1);
        mpz_set(*a1, *a2);
        mpz_set(*b0, *b1);
        mpz_set(*b1, *b2);
        //dr = f2p_degree(*r0);
        //int da = f2p_degree(*a0);
        //int db = f2p_degree(*b0);
        //if (dr < m) {
        //if (da <= m && db <= m && da >= 1 && db >= 1) {
        //break;
            //}
        PUTS("loop last\n");
    }
    PUTS("loop end\n");
    PRT("a2 = ", *a2);
    mpz_set(a, *a2);
    PRT("b2 = ", *b2);
    mpz_set(b, *b2);
    PRT("r2 = ", *r2);
    mpz_set(c, *r2);
    PUTS("f2p_exeuclid3 end\n");
}

void test_f2p_minpoly2(f2rng gen, int mexp)
{
    PUTS("minpoly start\n");
    f2p_wm_t wm;
    int wp = 0;
    f2p_wm_init(&wm, 25);
    mpz_t *poly = &(wm.ar[wp++]);
    mpz_t *seq = &(wm.ar[wp++]);
    mpz_t *b = &(wm.ar[wp++]);
    mpz_t *c = &(wm.ar[wp++]);
    mpz_t *x2t = &(wm.ar[wp++]);
    mpz_t *r1 = &(wm.ar[wp++]);
    mpz_t *r2 = &(wm.ar[wp++]);
    //mpz_t *r3 = &(wm.ar[wp++]);
    mpz_setbit(*x2t, 2 * mexp);
    //mpz_setbit(*x2t, 0);
//    for (int i = 0; i < 2 * mexp; i++) {
    for (int i = 2 * mexp - 1; i >= 0; i--) {
        unsigned int x = gen();
        if (x & 1) {
            mpz_setbit(*seq, i);
        } else {
            mpz_clrbit(*seq, i);
        }
    }
    PRT("seq = ", *seq);
    //f2p_exeuclid2(*poly, *b, *c, *seq, *x2t, mexp, wp, &wm); // 13 wm
    f2p_exeuclid3(*poly, *b, *c, *seq, *x2t, mexp, wp, &wm); // 13 wm
    PRT("poly = ", *poly);
    //PRT("b = ", *b);
    PRT("c = ", *c);
    PRT("seq = ", *seq);
    PRT("x2t = ", *x2t);
    f2p_mul(*r1, *poly, *seq, wp, &wm);
    f2p_mul(*r2, *b, *x2t, wp, &wm);
    f2p_add(*r1, *r1, *r2);
    int cmp = mpz_cmp(*r1, *c);
    if (cmp == 0) {
        printf("exeuclid3 OK\n");
    } else {
        printf("exeuclid3 NG cmp = %d\n", cmp);
    }
    char minpoly[20000];
    f2p_get_hexstr(minpoly, *poly);
    f2p_wm_clear(&wm);
    printf("minpoly = %s\n", minpoly);
    PUTS("minpoly end\n");
}

int test_minpoly2()
{
    test_f2p_minpoly2(dummy, 127);
    return 1;
}
#endif

int main(int argc, char * argv[])
{
    int verbose = 0;
    int r = 1;
    uint32_t seed = 1;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = 1;
    }
    if (argc > 2) {
        seed = strtoul(argv[2], NULL, 10);
    }
    r *= test_minpoly(seed, verbose);
    //r *= test_annihilate();
    //r *= test_minpoly2();
    r *= test_irreducible(verbose);
    if (r == 1) {
        return 0;
    } else {
        return -1;
    }
}
