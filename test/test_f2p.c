/**
 * @file test_f2p.c
 *
 * @brief test program for f2p_gmp.c
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
#include <string.h>

int test_string(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_string\n");
    }
    char * str1 = "101";
    char buffer[100];
    mpz_t *x = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    f2p_set_binstr(*x, str1);
    f2p_get_binstr(buffer, *x);
    int r;
    r = strcmp(str1, buffer);
    if (r != 0) {
        printf("test_string failure\n");
        printf("org str = |%s|\n", str1);
        printf("result str = |%s|\n", buffer);
    }
    if (verbose) {
        printf("end test_string\n");
    }
    if (r) {
        printf("x");
        return 0;
    } else {
        printf("o");
        return 1;
    }
}

int test_degree(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_degree\n");
    }
    //int r = 0;
    mpz_t *x = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);

    f2p_set_binstr(*x, "101");
    int deg = f2p_degree(*x);
    int ok = 1;
    if (deg != 2) {
        printf("test_degree failure expected 2 returns %d\n", deg);
        ok = 0;
    }
    f2p_set_binstr(*x, "1");
    deg = f2p_degree(*x);
    if (deg != 0) {
        printf("test_degree failure expected 0 returns %d\n", deg);
        ok = 0;
    }
    if (verbose) {
        printf("end test_degree\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_add(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_add\n");
    }
    //int r = 0;
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *y = &(wm->ar[wp++]);
    mpz_t *z = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);

    f2p_set_binstr(*x, "101");
    f2p_set_binstr(*y, "011");
    f2p_add(*z, *x, *y);
    char buffer[200];
    f2p_get_binstr(buffer, *z);
    if (strcmp(buffer, "110") != 0) {
        printf("test_add failure expected 110 returns %s\n", buffer);
        ok = 0;
    }
    f2p_set_binstr(*y, "101");
    f2p_add(*z, *x, *y);
    f2p_get_binstr(buffer, *z);
    if (strcmp(buffer, "0") != 0) {
        printf("test_add failure expected 0 returns %s\n", buffer);
        ok = 0;
    }
    if (verbose) {
        printf("end test_add\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }

    return ok;
}

int test_lshift(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_lshift\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);

    f2p_set_binstr(*x, "101");
    f2p_lshift(*x, 2);
    char buffer[200];
    f2p_get_binstr(buffer, *x);
    if (strcmp(buffer, "10100") != 0) {
        printf("test_lshift failure expected 1010 returns %s\n", buffer);
        ok = 0;
    }
    f2p_set_binstr(*x, "101");
    f2p_lshift(*x, 0);
    f2p_get_binstr(buffer, *x);
    if (strcmp(buffer, "101") != 0) {
        printf("test_lshift failure expected 0 returns %s\n", buffer);
        ok = 0;
    }
    if (verbose) {
        printf("end test_lshift\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }

    return ok;
}

int test_rshift(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_rshift\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);

    f2p_set_binstr(*x, "1011");
    f2p_rshift(*x, 2);
    char buffer[200];
    f2p_get_binstr(buffer, *x);
    if (strcmp(buffer, "10") != 0) {
        printf("test_rshift failure expected 1010 returns %s\n", buffer);
        ok = 0;
    }
    f2p_set_binstr(*x, "0");
    f2p_rshift(*x, 2);
    f2p_get_binstr(buffer, *x);
    if (strcmp(buffer, "0") != 0) {
        printf("test_rshift failure expected 0 returns %s\n", buffer);
        ok = 0;
    }
    f2p_set_binstr(*x, "11");
    f2p_rshift(*x, 2);
    f2p_get_binstr(buffer, *x);
    if (strcmp(buffer, "0") != 0) {
        printf("test_rshift failure expected 0 returns %s\n", buffer);
        ok = 0;
    }
    if (verbose) {
        printf("end test_rshift\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_coefficient(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_coefficient\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);

    f2p_set_binstr(*x, "1011");
    int c = f2p_coefficient(*x, 0);
    if (c != 1) {
        printf("test_cofficient failure 1 expected 1 returns %d\n", c);
        ok = 0;
    }
    c = f2p_coefficient(*x, 2);
    if (c != 0) {
        printf("test_cofficient failure 2 expected 0 returns %d\n", c);
        ok = 0;
    }
    c = f2p_coefficient(*x, 3);
    if (c != 1) {
        printf("test_cofficient failure 3 expected 1 returns %d\n", c);
        ok = 0;
    }
    c = f2p_coefficient(*x, 4);
    if (c != 0) {
        printf("test_cofficient failure 4 expected 0 returns %d\n", c);
        ok = 0;
    }
    if (verbose) {
        printf("end test_coefficient\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_mod(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_mod\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *y = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    char buff[200];

    f2p_set_binstr(*x, "1011");
    f2p_set_binstr(*y, "11");
    f2p_mod(*x, *y, wp, wm);
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "1") != 0) {
        printf("test_mod failure 1 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_set_binstr(*x, "11");
    f2p_set_binstr(*y, "1011");
    f2p_mod(*x, *y, wp, wm);
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_mod failure 2 expected 11 returns %s\n", buff);
        ok = 0;
    }
    if (verbose) {
        printf("end test_mod\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_divrem(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_divrem\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *y = &(wm->ar[wp++]);
    mpz_t *q = &(wm->ar[wp++]);
    mpz_t *rem = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    char buff[200];

    f2p_set_binstr(*x, "1011");
    f2p_set_binstr(*y, "11");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_binstr(buff, *q);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 1 expected 110 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 2 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "1011") != 0) {
        printf("test_divrem failure 2 expected 1011 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *y);
    if (strcmp(buff, "11") != 0) {
        printf("test_divrem failure 2 expected 11 returns %s\n", buff);
        ok = 0;
    }
    f2p_set_binstr(*x, "11");
    f2p_set_binstr(*y, "1011");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_binstr(buff, *q);
    if (strcmp(buff, "0") != 0) {
        printf("test_divrem failure 3 expected 0 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "11") != 0) {
        printf("test_divrem failure 4 expected 11 returns %s\n", buff);
        ok = 0;
    }
    f2p_set_binstr(*x, "1101");
    f2p_set_binstr(*y, "1011");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_binstr(buff, *q);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 5 expected 0 returns %s\n", buff);
        ok = 0;
        ok = 0;
    }
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 6 expected 110 returns %s\n", buff);
        ok = 0;
    }
    f2p_set_binstr(*x, "1011");
    f2p_set_binstr(*y, "11");
    f2p_divrem(*x, *y, *x, *y, wp, wm);
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 7 expected 0 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *y);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 8 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_set_binstr(*x, "1011");
    f2p_set_binstr(*y, "11");
    f2p_divrem(*y, *x, *x, *y, wp, wm);
    f2p_get_binstr(buff, *y);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 9 expected 0 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 10 expected 1 returns %s\n", buff);
        ok = 0;
    }

    f2p_set_binstr(*x, "1011");
    f2p_set_binstr(*y, "1");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_binstr(buff, *q);
    if (strcmp(buff, "1011") != 0) {
        printf("test_divrem failure 11 expected 1011 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "0") != 0) {
        printf("test_divrem failure 12 expected 0 returns %s\n", buff);
        ok = 0;
    }

    f2p_set_binstr(*x, "1110110");
    f2p_set_binstr(*y, "1001");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_binstr(buff, *q);
    if (strcmp(buff, "1111") != 0) {
        printf("test_divrem failure 13 expected 1111 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 14 expected 1 returns %s\n", buff);
        ok = 0;
    }

    if (verbose) {
        printf("end test_divrem\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_mulmod(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_mulmod\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *y = &(wm->ar[wp++]);
    mpz_t *mod = &(wm->ar[wp++]);
    mpz_t *rem = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    char buff[200];

    f2p_set_binstr(*x, "101");
    f2p_set_binstr(*y, "11");
    f2p_set_binstr(*mod, "10");
    f2p_mulmod(*rem, *x, *y, *mod, wp, wm);
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_mulmod failure 1 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *mod);
    if (strcmp(buff, "10") != 0) {
        printf("test_mulmod failure 2 expected 10 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "101") != 0) {
        printf("test_mulmod failure 3 expected 101 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *y);
    if (strcmp(buff, "11") != 0) {
        printf("test_mulmod failure 4 expected 11 returns %s\n", buff);
        ok = 0;
    }

    f2p_set_binstr(*x, "11");
    f2p_set_binstr(*y, "1011");
    f2p_set_binstr(*mod, "101");
    f2p_mulmod(*x, *x, *y, *mod, wp, wm);
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_mulmod failure 5 expected 11 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *y);
    if (strcmp(buff, "1011") != 0) {
        printf("test_mulmod failure 6 expected 1011 returns %s\n", buff);
        ok = 0;
    }
    if (verbose) {
        printf("end test_mulmod\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_powermod(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_powermod\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *e = &(wm->ar[wp++]);
    mpz_t *mod = &(wm->ar[wp++]);
    mpz_t *rem = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    char buff[200];

    f2p_set_binstr(*x, "11");
    mpz_set_ui(*e, 2);
    f2p_set_binstr(*mod, "111");
    f2p_powermod(*rem, *x, *e, *mod, wp, wm);
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "10") != 0) {
        printf("test_powermod failure 1 expected 10 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *mod);
    if (strcmp(buff, "111") != 0) {
        printf("test_powermod failure 2 expected 111 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_powermod failure 3 expected 11 returns %s\n", buff);
        ok = 0;
    }
    int d = mpz_cmp_ui(*e, 2);
    if (d != 0) {
        printf("test_powermod failure 4 expected 0 returns %d\n", d);
        ok = 0;
    }
    f2p_set_binstr(*x, "11");
    mpz_set_ui(*e, 3);
    f2p_set_binstr(*mod, "111");
    f2p_powermod(*rem, *x, *e, *mod, wp, wm);
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_powermod failure 5 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *mod);
    if (strcmp(buff, "111") != 0) {
        printf("test_powermod failure 6 expected 111 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_powermod failure 7 expected 11 returns %s\n", buff);
        ok = 0;
    }
    d = mpz_cmp_ui(*e, 3);
    if (d != 0) {
        printf("test_powermod failure 8 expected 0 returns %d\n", d);
        ok = 0;
    }

    f2p_set_binstr(*x, "11");
    mpz_set_ui(*e, 5);
    f2p_set_binstr(*mod, "111");
    f2p_powermod(*rem, *x, *e, *mod, wp, wm);
    f2p_get_binstr(buff, *rem);
    if (strcmp(buff, "10") != 0) {
        printf("test_powermod failure 10 expected 10 returns %s\n", buff);
        ok = 0;
    }
    if (verbose) {
        printf("end test_powermod\n");
    }
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int test_gcd_aux(const char * binx, const char * biny, const char * bingcd,
                 int verbose, int wp, f2p_wm_t *wm)
{
    static char buff[2000];
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *y = &(wm->ar[wp++]);
    mpz_t *expected = &(wm->ar[wp++]);
    mpz_t *gcd = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    f2p_set_binstr(*x, binx);
    f2p_set_binstr(*y, biny);
    f2p_set_binstr(*expected, bingcd);
    f2p_gcd(*gcd, *x, *y, wp, wm);
    int r = mpz_cmp(*gcd, *expected);
    if (r != 0 && verbose) {
        f2p_get_binstr(buff, *gcd);
        printf("gcd(%s, %s) expected = %s, returns %s\n", binx, biny, bingcd,
                buff);
        fflush(stdout);
        return 0;
    } else {
        return 1;
    }
}

int test_gcd(int verbose, int wp, f2p_wm_t *wm)
{
    int ok = 1;
    ok *= test_gcd_aux("10", "1", "1", verbose, wp, wm);
    if (ok) {
        printf("o");
    } else {
        printf("x");
    }
    return ok;
}

int main(int argc, char * argv[])
{
    int verbose = 0;
    int ok = 1;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = 1;
    }
    int wp = 0;
    f2p_wm_t wm;

    f2p_wm_init(&wm, 20);

    ok *= test_string(verbose, wp, &wm);
    ok *= test_degree(verbose, wp, &wm);
    ok *= test_add(verbose, wp, &wm);
    ok *= test_lshift(verbose, wp, &wm);
    ok *= test_rshift(verbose, wp, &wm);
    ok *= test_coefficient(verbose, wp, &wm);
    ok *= test_mod(verbose, wp, &wm);
    ok *= test_divrem(verbose, wp, &wm);
    ok *= test_mulmod(verbose, wp, &wm);
    ok *= test_powermod(verbose, wp, &wm);
    ok *= test_gcd(verbose, wp, &wm);
    printf("\n");
    f2p_wm_clear(&wm);

    if (ok == 1) {
        return 0;
    } else {
        return -1;
    }
}
