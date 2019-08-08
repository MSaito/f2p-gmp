#pragma once
#ifndef F2P_GMP_H
#define F2P_GMP_H
/*
 */

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <assert.h>

#define F2P_WMSIZE 10

static inline void f2p_set_str(mpz_t poly, const char * str)
{
    mpz_set_str(poly, str, 2);
}

static inline char * f2p_get_str(char * str, mpz_t poly)
{
    return mpz_get_str(str, 2, poly);
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

void f2p_wm_init(mpz_t *wm);

void f2p_wm_clear(mpz_t *wm);

void f2p_mod(mpz_t a, mpz_t b, int wp, mpz_t *wm);

void f2p_divrem(mpz_t q, mpz_t r, mpz_t a, mpz_t b, int wp, mpz_t *wm);

void f2p_mulmod(mpz_t r, mpz_t a, mpz_t b, mpz_t mod, int wp, mpz_t *wm);

void f2p_powermod(mpz_t r, mpz_t x, mpz_t e, mpz_t mod, int wp, mpz_t *wm);

#endif // F2P_GMP_H
