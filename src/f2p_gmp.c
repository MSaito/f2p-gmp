/**
 * @file f2p_gmp.c
 *
 * @brief Simple F2 Polynomial Library for GMP (GNU Multi-Precision Library).
 *
 * @author Mutsuo Saito
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2019 Mutsuo Saito, Makoto Matsumoto
 * and Hiroshima University.
 * All rights reserved.
 *
 * The MIT License is applied to this software, see
 * LICENSE.txt
 */

#include "f2p_gmp.h"
#include "test_f2p.h"
#include <stdlib.h> // malloc

static void f2p_mul_aux(mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm);
static void f2p_mulmod_aux(mpz_t r, mpz_t a, mpz_t b, mpz_t mod,
                           int wp, f2p_wm_t *wm);

/**
 * Initialization of (Shared) Working Memory
 *
 * @param wm working memory
 */
void f2p_wm_init(f2p_wm_t *wm, int size)
{
    wm->max_size = size;
    wm->ar = malloc(size * sizeof(mpz_t));
    assert(wm->ar != NULL);
    for (int i = 0; i < size; i++) {
        mpz_init(wm->ar[i]);
    }
}

/**
 * Cearing (Shared) Working Memory
 *
 * @param wm working memory
 */
void f2p_wm_clear(f2p_wm_t *wm)
{
    for (int i = 0; i < wm->max_size; i++) {
        mpz_clear(wm->ar[i]);
    }
    free(wm->ar);
}

/**
 * Calculate residue of polynomial a divided by b polynomial.
 * a %= b
 *
 * use 1 wm
 * @param a dividend and result
 * @param b divisor
 * @param wp pointer of wm
 * @param wm shared working memory
 */
//void f2p_mod(mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm)
void f2p_mod(mpz_t a, mpz_t b, int wp, f2p_wm_t *wm)
{
    int zcmp = mpz_cmp_ui(b, 0);
    assert(zcmp != 0); // zero divide
    int deg = f2p_degree(b);
    int diff = f2p_degree(a) - deg;
    //mpz_t *z = &wm[wp++];
    mpz_t *x = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    if (mpz_cmp_ui(b, 1) == 0) {
        mpz_set_ui(a, 0);
        return;
    } else if (diff < 0) {
        //mpz_set(r, a);
        return;
    } else if (diff == 0) {
        f2p_add(a, a, b);
        return;
    }
    //mpz_set(*wz, a);
    mpz_set(*x, b);
    f2p_lshift(*x, diff);
    f2p_add(a, a, *x);
    int zdeg = f2p_degree(a);
    PRT("a = ", a);
    while (zdeg >= deg) {
        f2p_rshift(*x, f2p_degree(*x) - zdeg);
        f2p_add(a, a, *x);
        PRT("a = ", a);
        PRT("x = ", *x);
        zdeg = f2p_degree(a);
    }
    //mpz_set(r, *wz);
}

/**
 * Calculate r and q such that
 * a = q * b + r (degree(r) < degree(b))
 *
 * use 2 wm
 *
 * @param q quotient
 * @param r remainder
 * @param a dividend
 * @param b divisor
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_divrem(mpz_t q, mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm)
{
    int zcmp = mpz_cmp_ui(b, 0);
    assert(zcmp != 0); // zero divide
    int deg = f2p_degree(b);
    int diff = f2p_degree(a) - deg;
    PUTS("in f2p_divrem");
    PRT("a = ", a);
    PRT("b = ", b);
    if (diff < 0) {
        mpz_set(r, a);
        mpz_set_ui(q, 0);
        PUTS("out 1 f2p_divrem");
        return;
    } else if (diff == 0) {
        f2p_add(r, a, b);
        mpz_set_ui(q, 1);
        PUTS("out 2 f2p_divrem");
        return;
    } else if (mpz_cmp_ui(b, 1) == 0) {
        mpz_set_ui(r, 0);
        mpz_set(q, a);
        PUTS("out 2.5 f2p_divrem");
        return;
    }
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *qw = &(wm->ar[wp++]);
    //mpz_t *rw = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);

    mpz_set_ui(*qw, 0);
    mpz_set(*x, b);
    mpz_set(r, a);
    int pos = diff;
    f2p_lshift(*x, diff);
    mpz_setbit(*qw, pos);
    f2p_add(r, r, *x);
    int zdeg = f2p_degree(r);
    PRT("r = ", r);
    PRT("x = ", *x);
    PRT("qw = ", *qw);
    while (zdeg >= deg) {
        int d = f2p_degree(*x) - zdeg;
        pos = pos - d;
        f2p_rshift(*x, d);
        mpz_setbit(*qw, pos);
        f2p_add(r, r, *x);
        zdeg = f2p_degree(r);
        PRT("r = ", r);
        PRT("x = ", *x);
        PRT("qw = ", *qw);
    }
    mpz_set(q, *qw);
    PUTS("out 3 f2p_divrem");
    //mpz_set(r, *rw);
}

/**
 * Calculate r = a * b
 *
 * use 2 wm
 *
 * @param r result
 * @param a polynomial
 * @param b polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_mul(mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm)
{
    if (f2p_degree(a) < f2p_degree(b)) {
        f2p_mul_aux(r, b, a, wp, wm);
    } else {
        f2p_mul_aux(r, a, b, wp, wm);
    }
}

/**
 * Calculate r = a * b
 *
 * use 2 wm
 * @param r result
 * @param a polynomial
 * @param b polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
static void f2p_mul_aux(mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm)
{
    mpz_t *v = &(wm->ar[wp++]);
    mpz_t *w = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    mpz_set(*v, a);
    mpz_set(*w, b);
    mpz_set_ui(r, 0); // r should be set 0 after a is copy to v
    mp_bitcnt_t bposmax = f2p_degree(*w);
    mp_bitcnt_t bpos = 0;
    while(bpos <= bposmax) {
        if (f2p_coefficient(*w, bpos) == 1) {
            f2p_add(r, r, *v);
        }
        bpos++;
        f2p_lshift(*v, 1); // bottle neck
    }
}

/**
 * Calculate r = (a * b) % mod
 *
 * use 3 wm
 *
 * @param r result
 * @param a polynomial
 * @param b polynomial
 * @param mod modulus polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_mulmod(mpz_t r, mpz_t a, mpz_t b, mpz_t mod, int wp, f2p_wm_t *wm)
{
    if (f2p_degree(a) < f2p_degree(b)) {
        f2p_mulmod_aux(r, b, a, mod, wp, wm);
    } else {
        f2p_mulmod_aux(r, a, b, mod, wp, wm);
    }
}

/**
 * Calculate r = (a * b) % mod
 *
 * use 3 wm
 *
 * @param r result
 * @param a polynomial
 * @param b polynomial
 * @param mod modulus polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
static void f2p_mulmod_aux(mpz_t r, mpz_t a, mpz_t b, mpz_t mod,
                           int wp, f2p_wm_t *wm)
{
    //mp_bitcnt_t maxbits = mpz_sizeinbase(mod, 2);
    int zcmp = mpz_cmp_ui(mod, 0);
    assert(zcmp != 0); // zero divide
    mpz_t *v = &(wm->ar[wp++]);
    mpz_t *w = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    mpz_set(*v, a);
    mpz_set(*w, b);
    mpz_set_ui(r, 0); // r should be set 0 after a is copy to v
    f2p_mod(*v, mod, wp, wm); // 1 wm
    f2p_mod(*w, mod, wp, wm);
    mp_bitcnt_t bposmax = mpz_sizeinbase(*w, 2) - 1;
    mp_bitcnt_t bpos = 0;
    while(bpos <= bposmax) {
        if (f2p_coefficient(*w, bpos) == 1) {
            f2p_add(r, r, *v);
        }
        bpos++;
        f2p_lshift(*v, 1); // bottle neck
        if (f2p_degree(*v) == f2p_degree(mod)) {
            f2p_add(*v, *v, mod);
        }
    }
}

/**
 * Calculate r = (a * a) % mod
 *
 * use 3 wm
 *
 * @param r result
 * @param a polynomial
 * @param mod modulus polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_pow2mod(mpz_t r, mpz_t a, mpz_t mod, int wp, f2p_wm_t *wm)
{
    int zcmp = mpz_cmp_ui(mod, 0);
    assert(zcmp != 0); // zero divide
    mp_bitcnt_t bpos = 0;
    mp_bitcnt_t bposmax = f2p_degree(a);
    mpz_t *b = &(wm->ar[wp++]);
    mpz_t *v = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    mpz_set(*b, a);
    mpz_set(*v, a);
    mpz_set_ui(r, 0); // r should be set 0 after a is copied.
    f2p_mod(*v, mod, wp, wm); // 1 wm
    while(bpos <= bposmax) {
        if (f2p_coefficient(*b, bpos) == 1) {
            f2p_add(r, r, *v);
        }
        bpos++;
        f2p_lshift(*v, 1);
        if (f2p_degree(*v) == f2p_degree(mod)) {
            f2p_add(*v, *v, mod);
        }
    }
}

/**
 * calculate r = x^e % mod.
 *
 * use 4 wm
 *
 * @param r residue polynomial whose degree is less than mod polynomial
 * @param x polynomial
 * @param e exponent (big integer)
 * @param mod polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_powermod(mpz_t r, mpz_t x, mpz_t e, mpz_t mod, int wp, f2p_wm_t *wm)
{
    int zcmp = mpz_cmp_ui(mod, 0);
    assert(zcmp != 0); // zero divide
    mpz_t *s = &(wm->ar[wp++]);
    //mpz_t *z = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    mpz_set(*s, x);
    //mpz_set_ui(*z, 1);
    mpz_set_ui(r, 1);
    mp_bitcnt_t bposmax = mpz_sizeinbase(e, 2) - 1;
    mp_bitcnt_t bpos = 0;
    while (bpos <= bposmax) {
        if (mpz_tstbit(e, bpos) == 1) {
            f2p_mulmod(r, r, *s, mod, wp, wm); // 3 wm
        }
        f2p_pow2mod(*s, *s, mod, wp, wm); // 3 wm
        bpos++;
    }
    //mpz_set(r, *z);
}

/**
 * extended euclid
 * calculate a, b, c for given x, y
 * where a*x + b*y = c (c = GCD(x, y))
 *
 * use 13 wm
 *
 *@param a result polynomial
 *@param b result polynomial
 *@param c result polynomial GCD(x, y)
 *@param x input polynomial
 *@param y input polynomial
 *@param wp
 *@param wm
 */
void f2p_exeuclid(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y,
                  int wp, f2p_wm_t *wm)
{
    PUTS("f2p_exeuclid start\n");
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
    while (mpz_cmp_ui(*r1, 0) > 0) {
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
        PRT("r0 = ", *r0);
        PRT("r1 = ", *r1);
        PRT("r2 = ", *r2);
        PRT("a0 = ", *a0);
        PRT("a1 = ", *a1);
        PRT("a2 = ", *a2);
        PRT("b0 = ", *b0);
        PRT("b1 = ", *b1);
        PRT("b2 = ", *b2);
        mpz_set(*r0, *r1);
        mpz_set(*r1, *r2);
        mpz_set(*a0, *a1);
        mpz_set(*a1, *a2);
        mpz_set(*b0, *b1);
        mpz_set(*b1, *b2);
        PUTS("loop last\n");
    }
    PUTS("loop end\n");
    PRT("a0 = ", *a0);
    mpz_set(a, *a0);
    PRT("b0 = ", *b0);
    mpz_set(b, *b0);
    PRT("r0 = ", *r0);
    mpz_set(c, *r0);
    PUTS("f2p_exeuclid end\n");
}

/**
 * extended euclid
 * calculate a, b, c for given x, y
 * where a*x + b*y = c (c = GCD(x, y))
 *
 * use 10 wm
 *
 *@see https://planetmath.org/berlekampmasseyalgorithm
 *
 *@param a result polynomial
 *@param c result polynomial GCD(x, y)
 *@param x input polynomial
 *@param y input polynomial
 *@param m max degree
 *@param wp
 *@param wm
 */
//void f2p_exeuclid2(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y, int m,
//                  int wp, f2p_wm_t *wm)
void f2p_exeuclid2(mpz_t a, mpz_t c, mpz_t x, mpz_t y, int m,
                   int wp, f2p_wm_t *wm)
{
    PUTS("f2p_exeuclid2 start\n");
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
    //mpz_t *b0 = &(wm->ar[wp++]);
    //mpz_t *b1 = &(wm->ar[wp++]);
    //mpz_t *b2 = &(wm->ar[wp++]);
    mpz_t *tmp = &(wm->ar[wp++]);
    //mpz_t *t; 後で使う
    assert(wp <= wm->max_size);
    PRT("x= ", x);
    PRT("y= ", y);
    mpz_set(*r0, x);
    mpz_set(*r1, y);
    mpz_set_ui(*a0, 1);
    mpz_set_ui(*a1, 0);
    //mpz_set_ui(*b0, 0);
    //mpz_set_ui(*b1, 1);
    PRT("r0 = ", *r0);
    PRT("r1 = ", *r1);
    PRT("a0 = ", *a0);
    PRT("a1 = ", *a1);
    //PRT("b0 = ", *b0);
    //PRT("b1 = ", *b1);
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
        //f2p_mul(*tmp, *q1, *b1, wp, wm);// 2 wm
        //PUTS("f2p_add\n");
        //f2p_add(*b2, *b0, *tmp);
        PRT("q1 = ", *q1);
        PRT("tmp = ", *tmp);
        PRT("r0 = ", *r0);
        PRT("r1 = ", *r1);
        PRT("r2 = ", *r2);
        PRT("a0 = ", *a0);
        PRT("a1 = ", *a1);
        PRT("a2 = ", *a2);
        //PRT("b0 = ", *b0);
        //PRT("b1 = ", *b1);
        //PRT("b2 = ", *b2);
        int dr = f2p_degree(*r2);
        if (dr < m) {
            break;
        }
        mpz_set(*r0, *r1);
        mpz_set(*r1, *r2);
        mpz_set(*a0, *a1);
        mpz_set(*a1, *a2);
        //mpz_set(*b0, *b1);
        //mpz_set(*b1, *b2);
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
    //PRT("b2 = ", *b2);
    //mpz_set(b, *b2);
    PRT("r2 = ", *r2);
    mpz_set(c, *r2);
    PUTS("f2p_exeuclid2 end\n");
}

#if 0
void f2p_minpoly(char * minpoly, f2rng gen, int mexp)
{
    PUTS("minpoly start\n");
    f2p_wm_t wm;
    int wp = 0;
    f2p_wm_init(&wm, 20);
    mpz_t *poly = &(wm.ar[wp++]);
    mpz_t *seq = &(wm.ar[wp++]);
    //mpz_t *b = &(wm.ar[wp++]);
    mpz_t *c = &(wm.ar[wp++]);
    mpz_t *x2t = &(wm.ar[wp++]);
    mpz_setbit(*x2t, 2 * mexp);
    //mpz_realloc2(*seq, 2 * mexp);
    //mpz_setbit(*x2t, 0);
//    for (int i = 0; i < 2 * mexp; i++) {
    for (int i = 2 * mexp - 1; i >= 0; i--) {
        uint32_t x = gen();
        if (x & 1) {
            mpz_setbit(*seq, i);
        } else {
            mpz_clrbit(*seq, i);
        }
    }
    PRT("seq = ", *seq);
    //f2p_exeuclid2(*poly, *b, *c, *seq, *x2t, mexp, wp, &wm); // 13 wm
    f2p_exeuclid2(*poly, *c, *seq, *x2t, mexp, wp, &wm); // 13 wm
    PRT("poly = ", *poly);
    //PRT("b = ", *b);
    PRT("c = ", *c);
    PRT("seq = ", *seq);
    PRT("x2t = ", *x2t);
    f2p_get_hexstr(minpoly, *poly);
    f2p_wm_clear(&wm);
    PUTS("minpoly end\n");
}
#endif

/**
 * gcd
 * calculate gcd(x, y)
 *
 * use 5 wm
 *
 *@param gcd result polynomial
 *@param x input polynomial
 *@param y input polynomial
 *@param wp
 *@param wm
 */
void f2p_gcd(mpz_t gcd, mpz_t x, mpz_t y, int wp, f2p_wm_t *wm)
{
    PUTS("f2p_gcd start\n");
    int zcmp = mpz_cmp_ui(x, 0);
    assert(zcmp != 0);
    zcmp = mpz_cmp_ui(y, 0);
    assert(zcmp != 0);
    mpz_t *r0 = &(wm->ar[wp++]);
    mpz_t *r1 = &(wm->ar[wp++]);
    mpz_t *r2 = &(wm->ar[wp++]);
    //mpz_t *t; 後で使う
    assert(wp <= wm->max_size);
    PRT("x= ", x);
    PRT("y= ", y);
    mpz_set(*r0, x);
    mpz_set(*r1, y);
    PRT("r0 = ", *r0);
    PRT("r1 = ", *r1);
    while (mpz_cmp_ui(*r1, 0) > 0) {
        PUTS("f2p_mod\n");
        mpz_set(*r2, *r0);
        //f2p_mod(*r2, *r0, *r1, wp, wm); // 1 wm
        f2p_mod(*r2, *r1, wp, wm); // 1 wm
        PRT("r0 = ", *r0);
        PRT("r1 = ", *r1);
        //mpz_swap(*r0, *r1);
        mpz_set(*r0, *r1);
        mpz_set(*r1, *r2);
        PUTS("loop last\n");
    }
    PUTS("loop end\n");
    PRT("r0 = ", *r0);
    mpz_set(gcd, *r0);
    PUTS("f2p_gcd end\n");

}

/**
 * is irreducible aux
 *
 * use 8 wm
 *
 *@param poly input polynomial
 *@param wp
 *@param wm
 *@reaturns 1 irreducible, 0 reducible
 */
int f2p_is_irreducible_aux(mpz_t poly, int wp, f2p_wm_t * wm)
{
    PUTS("f2p_is_irreducible_aux start\n");
    if (mpz_cmp_ui(poly, 0) == 0) {
        return 0;
    }
    if (mpz_cmp_ui(poly, 2) == 0 || mpz_cmp_ui(poly, 3) == 0) {
        return 1;
    }
    mpz_t *t2m = &(wm->ar[wp++]);
    mpz_t *t1 = &(wm->ar[wp++]);
    mpz_t *t = &(wm->ar[wp++]);
    mpz_t *work = &(wm->ar[wp++]);
    mpz_setbit(*t2m, 2); // t^2m
    mpz_setbit(*t1, 1);  // t^1
    int degpol = f2p_degree(poly);
    int result = 1;
    f2p_add(*t, *t2m, *t1);
    for (int m = 1; m <= degpol / 2; m++) {
        f2p_gcd(*work, poly, *t, wp, wm); // 4 wm
        if (mpz_cmp_ui(*work, 1) != 0) {
            result = 0;
            break;
        }
        f2p_mul(*t2m, *t2m, *t2m, wp, wm); // 2 wm
        f2p_mod(*t2m, poly, wp, wm); // 1 wm
        f2p_add(*t, *t2m, *t1);
    }
    PUTS("f2p_is_irreducible_aux end\n");
    return result;
}

#if 0
int f2p_is_irreducible(const char * poly)
{
    PUTS("f2p_is_irreducible start\n");
    f2p_wm_t wm;
    int wp = 0;
    f2p_wm_init(&wm, 10);
    mpz_t *pol = &(wm.ar[wp++]);
    f2p_set_hexstr(*pol, poly);
    int result = f2p_is_irreducible_aux(pol, wp, &wm);
    f2p_wm_clear(&wm);
    PUTS("f2p_is_irreducible end\n");
    return result;
}
#endif

/**
 * minimal polynomial of a linear recurrence sequence.
 *
 * sequence length len shoud be len >= 2 * maxdeg.
 *
 *@param minpoly
 *@param seq linear recurrence sequence
 *@param maxdeg supporsed max degree of minpoly
 */
void f2p_minpoly(mpz_t minpoly, mpz_t seq, int maxdeg)
{
    PUTS("minpoly start\n");
    f2p_wm_t wm;
    int wp = 0;
    f2p_wm_init(&wm, 16);
    //mpz_t *poly = &(wm.ar[wp++]);
    mpz_t *rseq = &(wm.ar[wp++]);
    //mpz_t *b = &(wm.ar[wp++]);
    mpz_t *c = &(wm.ar[wp++]);
    mpz_t *x2t = &(wm.ar[wp++]);
    mpz_setbit(*x2t, 2 * maxdeg); // x^{2*maxdeg}
    //mpz_realloc2(*seq, 2 * maxdeg);
    //mpz_setbit(*x2t, 0);
//    for (int i = 0; i < 2 * maxdeg; i++) {
    for (int i = 0; i < 2 * maxdeg; i++) {
        if (mpz_tstbit(seq, i)) {
            mpz_setbit(*rseq, 2 * maxdeg - 1 - i);
        } else {
            mpz_clrbit(*rseq, 2 * maxdeg - 1 - i);
        }
    }
    PRT("rseq = ", *rseq);
    //f2p_exeuclid2(*poly, *b, *c, *seq, *x2t, maxdeg, wp, &wm); // 13 wm
    f2p_exeuclid2(minpoly, *c, *rseq, *x2t, maxdeg, wp, &wm); // 13 wm
    PRT("poly = ", minpoly);
    //PRT("b = ", *b);
    PRT("c = ", *c);
    PRT("seq = ", seq);
    PRT("x2t = ", *x2t);
    //f2p_get_hexstr(minpoly, *poly);
    f2p_wm_clear(&wm);
    PUTS("minpoly end\n");
}

/**
 * is irreducible
 *
 * If a polynomial is irreducible and its degree is mersenne expornent
 * then the polynomial is primitive.
 *
 *@param poly input polynomial
 *@reaturns 1 if poly is irreducible, 0 otherwise
 */
int f2p_is_irreducible(mpz_t poly)
{
    PUTS("f2p_is_irreducible start\n");
    f2p_wm_t wm;
    int wp = 0;
    f2p_wm_init(&wm, 10);
    int result = f2p_is_irreducible_aux(poly, wp, &wm);
    f2p_wm_clear(&wm);
    PUTS("f2p_is_irreducible end\n");
    return result;
}
