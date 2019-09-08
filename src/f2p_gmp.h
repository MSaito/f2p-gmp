#pragma once
#ifndef F2P_GMP_H
#define F2P_GMP_H
/**
 * @file f2p_gmp.h
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

#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <assert.h>

#if defined(__cplusplus)
extern "C" {
#endif

    struct F2P_WM_T {
        int max_size;
        mpz_t * ar;
    };

    typedef struct F2P_WM_T f2p_wm_t;

//typedef unsigned int (*f2rng)(void);

    static inline void f2p_set_binstr(mpz_t poly, const char * str)
    {
        mpz_set_str(poly, str, 2);
    }

    static inline void f2p_set_hexstr(mpz_t poly, const char * str)
    {
        mpz_set_str(poly, str, 16);
    }

    static inline char * f2p_get_binstr(char * str, mpz_t poly)
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
    static inline size_t f2p_degree(mpz_t poly)
    {
#if 0
        if (mpz_cmp_ui(poly, 0) == 0) {
            return -1;
        }
#endif
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

    void f2p_pow2mod(mpz_t r, mpz_t a, mpz_t mod, int wp, f2p_wm_t *wm);

    void f2p_powermod(mpz_t r, mpz_t x, mpz_t e, mpz_t mod,
                      int wp, f2p_wm_t *wm);

    void f2p_square(mpz_t r, mpz_t a, int wp, f2p_wm_t *wm);

    void f2p_exeuclid(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y,
                      int wp, f2p_wm_t *wm);

//void f2p_exeuclid2(mpz_t a, mpz_t b, mpz_t c, mpz_t x, mpz_t y, int m,
//                   int wp, f2p_wm_t *wm);
    void f2p_exeuclid2(mpz_t a, mpz_t c, mpz_t x, mpz_t y, int m,
                       int wp, f2p_wm_t *wm);

    void f2p_gcd(mpz_t gcd, mpz_t a, mpz_t b, int wp, f2p_wm_t *wm);

//void f2p_minpoly(char * minpoly, f2rng gen, int mexp);

    void f2p_minpoly(mpz_t minpoly, mpz_t seq, int mexp);

//    int f2p_is_irreducible(const char * poly);
    int f2p_is_irreducible(mpz_t poly);

    void f2p_calc_jump(mpz_t jump, mpz_t minpoly, mpz_t step);

#if defined(__cplusplus)
}
#endif

#endif // F2P_GMP_H
