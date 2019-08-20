#pragma once
#ifndef F2P_GMP_H
#define F2P_GMP_H
/*
 */

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <assert.h>

//#define F2P_WMSIZE 20
struct F2P_WM_T {
    int max_size;
    mpz_t * ar;
};
typedef struct F2P_WM_T f2p_wm_t;

typedef unsigned int (*f2rng)(void);

static inline void f2p_set_str(mpz_t poly, const char * str)
{
    mpz_set_str(poly, str, 2);
}

static inline void f2p_set_hexstr(mpz_t poly, const char * str)
{
    mpz_set_str(poly, str, 16);
}

static inline char * f2p_get_str(char * str, mpz_t poly)
{
    return mpz_get_str(str, 2, poly);
}

static inline char * f2p_get_hexstr(char * str, mpz_t poly)
{
    return mpz_get_str(str, 16, poly);
}

/**
 * degree of polynomial
 *
 * this function returns degree 0 for the polynomial 0.
 *@param poly polynomial
 *@return degree of polynomial
 */
static inline mp_bitcnt_t f2p_degree(mpz_t poly)
{
    return mpz_sizeinbase(poly, 2) - 1;
}

static inline void f2p_add(mpz_t result, mpz_t a, mpz_t b)
{
    mpz_xor(result, a, b);
}

static inline void f2p_lshift(mpz_t result, unsigned long int n)
{
    mpz_mul_2exp(result, result, n);
}

static inline void f2p_rshift(mpz_t r, unsigned long int n)
{
    mpz_fdiv_q_2exp(r, r, n);
}

static inline int f2p_coefficient(mpz_t poly, size_t index)
{
    return mpz_tstbit(poly, index);
}

void f2p_wm_init(f2p_wm_t *wm, int size);

void f2p_wm_clear(f2p_wm_t *wm);

void f2p_mod(mpz_t a, mpz_t b, int wp, f2p_wm_t *wm);

void f2p_divrem(mpz_t q, mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm);

void f2p_mul(mpz_t r, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm);

void f2p_mulmod(mpz_t r, mpz_t a, mpz_t b, mpz_t mod, int wp, f2p_wm_t *wm);

void f2p_powermod(mpz_t r, mpz_t x, mpz_t e, mpz_t mod, int wp, f2p_wm_t *wm);

void f2p_exeuclid(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y,
                  int wp, f2p_wm_t *wm);

void f2p_exeuclid2(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y, int m,
                   int wp, f2p_wm_t *wm);

void minpoly(char * minpoly, f2rng gen, int mexp);

#endif // F2P_GMP_H
