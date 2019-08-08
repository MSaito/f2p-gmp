/*
 */

#include "f2p-gmp.h"

/**
 * Initialization of (Shared) Working Memory
 *
 * @param wm working memory
 */
void f2p_wm_init(mpz_t *wm)
{
    for (int i = 0; i < F2P_WMSIZE; i++) {
        mpz_init(wm[i]);
    }
}

/**
 * Cearing (Shared) Working Memory
 *
 * @param wm working memory
 */
void f2p_wm_clear(mpz_t *wm)
{
    for (int i = 0; i < F2P_WMSIZE; i++) {
        mpz_clear(wm[i]);
    }
}

/**
 * Calculate residue of polynomial a divided by b polynomial.
 * a %= b
 *
 * @param a dividend and result
 * @param b divisor
 * @param wp pointer of wm
 * @param wm shared working memory
 */
//void f2p_mod(mpz_t r, mpz_t a, mpz_t b, int wp, mpz_t *wm)
void f2p_mod(mpz_t a, mpz_t b, int wp, mpz_t *wm)
{
    int zcmp = mpz_cmp_ui(b, 0);
    assert(zcmp != 0); // zero divide
    int deg = f2p_degree(b);
    int diff = f2p_degree(a) - deg;
    //mpz_t *z = &wm[wp++];
    mpz_t *x = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    if (diff < 0) {
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
    while (zdeg >= deg) {
        f2p_rshift(*x, f2p_degree(*x) - zdeg);
        f2p_add(a, a, *x);
        zdeg = f2p_degree(a);
    }
    //mpz_set(r, *wz);
}

/**
 * Calculate r and q such that
 * a = q * b + r (degree(r) < degree(b))
 *
 * @param q quotient
 * @param r remainder
 * @param a dividend
 * @param b divisor
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_divrem(mpz_t q, mpz_t r, mpz_t a, mpz_t b, int wp, mpz_t *wm)
{
    int zcmp = mpz_cmp_ui(b, 0);
    assert(zcmp != 0); // zero divide
    int deg = f2p_degree(b);
    int diff = f2p_degree(a) - deg;
    if (diff < 0) {
        mpz_set(r, a);
        mpz_set_ui(q, 0);
        return;
    } else if (diff == 0) {
        f2p_add(r, a, b);
        mpz_set_ui(q, 1);
        return;
    }
    mpz_t *x = &wm[wp++];
    mpz_t *qw = &wm[wp++];
    //mpz_t *rw = &wm[wp++];
    assert(wp <= F2P_WMSIZE);

    mpz_set_ui(*qw, 0);
    mpz_set(*x, b);
    mpz_set(r, a);
    f2p_lshift(*x, diff);
    mpz_setbit(*qw, diff);
    f2p_add(r, r, *x);
    int zdeg = f2p_degree(r);
    while (zdeg >= deg) {
        int d = f2p_degree(*x) - zdeg;
        f2p_rshift(*x, d);
        mpz_setbit(*qw, diff - d);
        f2p_add(r, r, *x);
        zdeg = f2p_degree(r);
    }
    mpz_set(q, *qw);
    //mpz_set(r, *rw);
}

/**
 * Calculate r = (a * b) % mod
 *
 * @param r result
 * @param a polynomial
 * @param b polynomial
 * @param mod modulus polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_mulmod(mpz_t r, mpz_t a, mpz_t b, mpz_t mod, int wp, mpz_t *wm)
{
    //mp_bitcnt_t maxbits = mpz_sizeinbase(mod, 2);
    int zcmp = mpz_cmp_ui(mod, 0);
    assert(zcmp != 0); // zero divide
    mpz_t *v = &wm[wp++];
    mpz_t *w = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    mpz_set(*v, a);
    mpz_set(*w, b);
    mpz_set_ui(r, 0); // r should be set 0 after a is copy to v
    f2p_mod(*v, mod, wp, wm);
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
 * @param r result
 * @param a polynomial
 * @param mod modulus polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_pow2mod(mpz_t r, mpz_t a, mpz_t mod, int wp, mpz_t *wm)
{
    int zcmp = mpz_cmp_ui(mod, 0);
    assert(zcmp != 0); // zero divide
    mp_bitcnt_t bpos = 0;
    mp_bitcnt_t bposmax = f2p_degree(a);
    mpz_t *b = &wm[wp++];
    mpz_t *v = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    mpz_set(*b, a);
    mpz_set(*v, a);
    mpz_set_ui(r, 0); // r should be set 0 after a is copied.
    f2p_mod(*v, mod, wp, wm);
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
 * @param r residue polynomial whose degree is less than mod polynomial
 * @param x polynomial
 * @param e exponent (big integer)
 * @param mod polynomial
 * @param wp pointer of wm
 * @param wm shared working memory
 */
void f2p_powermod(mpz_t r, mpz_t x, mpz_t e, mpz_t mod, int wp, mpz_t *wm)
{
    int zcmp = mpz_cmp_ui(mod, 0);
    assert(zcmp != 0); // zero divide
    mpz_t *s = &wm[wp++];
    //mpz_t *z = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    mpz_set(*s, x);
    //mpz_set_ui(*z, 1);
    mpz_set_ui(r, 1);
    mp_bitcnt_t bposmax = mpz_sizeinbase(e, 2) - 1;
    mp_bitcnt_t bpos = 0;
    while (bpos <= bposmax) {
        if (mpz_tstbit(e, bpos) == 1) {
            f2p_mulmod(r, r, *s, mod, wp, wm);
        }
        f2p_pow2mod(*s, *s, mod, wp, wm);
        bpos++;
    }
    //mpz_set(r, *z);
}
