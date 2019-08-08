#include "f2p-gmp.h"
#include <string.h>

int test_string(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_string\n");
    }
    char * str1 = "101";
    char buffer[100];
    mpz_t *x = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    f2p_set_str(*x, str1);
    f2p_get_str(buffer, *x);
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
    return r;
}

int test_degree(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_degree\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    assert(wp <= F2P_WMSIZE);

    f2p_set_str(*x, "101");
    int deg = f2p_degree(*x);
    if (deg != 2) {
        printf("test_degree failure expected 2 returns %d\n", deg);
        r = 1;
    }
    f2p_set_str(*x, "1");
    deg = f2p_degree(*x);
    if (deg != 0) {
        printf("test_degree failure expected 0 returns %d\n", deg);
        r += 1;
    }
    if (verbose) {
        printf("end test_degree\n");
    }
    return r;
}

int test_add(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_add\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    mpz_t *y = &wm[wp++];
    mpz_t *z = &wm[wp++];
    assert(wp <= F2P_WMSIZE);

    f2p_set_str(*x, "101");
    f2p_set_str(*y, "011");
    f2p_add(*z, *x, *y);
    char buffer[200];
    f2p_get_str(buffer, *z);
    if (strcmp(buffer, "110") != 0) {
        printf("test_add failure expected 110 returns %s\n", buffer);
        r = 1;
    }
    f2p_set_str(*y, "101");
    f2p_add(*z, *x, *y);
    f2p_get_str(buffer, *z);
    if (strcmp(buffer, "0") != 0) {
        printf("test_add failure expected 0 returns %s\n", buffer);
        r += 1;
    }
    if (verbose) {
        printf("end test_add\n");
    }
    return r;
}

int test_lshift(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_lshift\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    assert(wp <= F2P_WMSIZE);

    f2p_set_str(*x, "101");
    f2p_lshift(*x, 2);
    char buffer[200];
    f2p_get_str(buffer, *x);
    if (strcmp(buffer, "10100") != 0) {
        printf("test_lshift failure expected 1010 returns %s\n", buffer);
        r = 1;
    }
    f2p_set_str(*x, "101");
    f2p_lshift(*x, 0);
    f2p_get_str(buffer, *x);
    if (strcmp(buffer, "101") != 0) {
        printf("test_lshift failure expected 0 returns %s\n", buffer);
        r += 1;
    }
    if (verbose) {
        printf("end test_lshift\n");
    }
    return r;
}

int test_rshift(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_rshift\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    assert(wp <= F2P_WMSIZE);

    f2p_set_str(*x, "1011");
    f2p_rshift(*x, 2);
    char buffer[200];
    f2p_get_str(buffer, *x);
    if (strcmp(buffer, "10") != 0) {
        printf("test_rshift failure expected 1010 returns %s\n", buffer);
        r = 1;
    }
    f2p_set_str(*x, "0");
    f2p_rshift(*x, 2);
    f2p_get_str(buffer, *x);
    if (strcmp(buffer, "0") != 0) {
        printf("test_rshift failure expected 0 returns %s\n", buffer);
        r += 1;
    }
    f2p_set_str(*x, "11");
    f2p_rshift(*x, 2);
    f2p_get_str(buffer, *x);
    if (strcmp(buffer, "0") != 0) {
        printf("test_rshift failure expected 0 returns %s\n", buffer);
        r += 1;
    }
    if (verbose) {
        printf("end test_rshift\n");
    }
    return r;
}

int test_coefficient(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_coefficient\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    assert(wp <= F2P_WMSIZE);

    f2p_set_str(*x, "1011");
    int c = f2p_coefficient(*x, 0);
    if (c != 1) {
        printf("test_cofficient failure 1 expected 1 returns %d\n", c);
        r = 1;
    }
    c = f2p_coefficient(*x, 2);
    if (c != 0) {
        printf("test_cofficient failure 2 expected 0 returns %d\n", c);
        r += 1;
    }
    c = f2p_coefficient(*x, 3);
    if (c != 1) {
        printf("test_cofficient failure 3 expected 1 returns %d\n", c);
        r += 1;
    }
    c = f2p_coefficient(*x, 4);
    if (c != 0) {
        printf("test_cofficient failure 4 expected 0 returns %d\n", c);
        r += 1;
    }
    if (verbose) {
        printf("end test_coefficient\n");
    }
    return r;
}

int test_mod(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_mod\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    mpz_t *y = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    char buff[200];

    f2p_set_str(*x, "1011");
    f2p_set_str(*y, "11");
    f2p_mod(*x, *y, wp, wm);
    f2p_get_str(buff, *x);
    if (strcmp(buff, "1") != 0) {
        printf("test_mod failure 1 expected 1 returns %s\n", buff);
        r = 1;
    }
    f2p_set_str(*x, "11");
    f2p_set_str(*y, "1011");
    f2p_mod(*x, *y, wp, wm);
    f2p_get_str(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_mod failure 2 expected 11 returns %s\n", buff);
        r += 1;
    }
    if (verbose) {
        printf("end test_mod\n");
    }
    return r;
}

int test_divrem(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_divrem\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    mpz_t *y = &wm[wp++];
    mpz_t *q = &wm[wp++];
    mpz_t *rem = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    char buff[200];

    f2p_set_str(*x, "1011");
    f2p_set_str(*y, "11");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_str(buff, *q);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 1 expected 110 returns %s\n", buff);
        r = 1;
    }
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 2 expected 1 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *x);
    if (strcmp(buff, "1011") != 0) {
        printf("test_divrem failure 2 expected 1011 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *y);
    if (strcmp(buff, "11") != 0) {
        printf("test_divrem failure 2 expected 11 returns %s\n", buff);
        r += 1;
    }
    f2p_set_str(*x, "11");
    f2p_set_str(*y, "1011");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_str(buff, *q);
    if (strcmp(buff, "0") != 0) {
        printf("test_divrem failure 3 expected 0 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "11") != 0) {
        printf("test_divrem failure 4 expected 11 returns %s\n", buff);
        r += 1;
    }
    f2p_set_str(*x, "1101");
    f2p_set_str(*y, "1011");
    f2p_divrem(*q, *rem, *x, *y, wp, wm);
    f2p_get_str(buff, *q);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 5 expected 0 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 6 expected 110 returns %s\n", buff);
        r += 1;
    }
    f2p_set_str(*x, "1011");
    f2p_set_str(*y, "11");
    f2p_divrem(*x, *y, *x, *y, wp, wm);
    f2p_get_str(buff, *x);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 7 expected 0 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *y);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 8 expected 1 returns %s\n", buff);
        r += 1;
    }
    f2p_set_str(*x, "1011");
    f2p_set_str(*y, "11");
    f2p_divrem(*y, *x, *x, *y, wp, wm);
    f2p_get_str(buff, *y);
    if (strcmp(buff, "110") != 0) {
        printf("test_divrem failure 9 expected 0 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *x);
    if (strcmp(buff, "1") != 0) {
        printf("test_divrem failure 10 expected 1 returns %s\n", buff);
        r += 1;
    }

    if (verbose) {
        printf("end test_divrem\n");
    }
    return r;
}

int test_mulmod(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_mulmod\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    mpz_t *y = &wm[wp++];
    mpz_t *mod = &wm[wp++];
    mpz_t *rem = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    char buff[200];

    f2p_set_str(*x, "101");
    f2p_set_str(*y, "11");
    f2p_set_str(*mod, "10");
    f2p_mulmod(*rem, *x, *y, *mod, wp, wm);
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_mulmod failure 1 expected 1 returns %s\n", buff);
        r = 1;
    }
    f2p_get_str(buff, *mod);
    if (strcmp(buff, "10") != 0) {
        printf("test_mulmod failure 2 expected 10 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *x);
    if (strcmp(buff, "101") != 0) {
        printf("test_mulmod failure 3 expected 101 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *y);
    if (strcmp(buff, "11") != 0) {
        printf("test_mulmod failure 4 expected 11 returns %s\n", buff);
        r += 1;
    }

    f2p_set_str(*x, "11");
    f2p_set_str(*y, "1011");
    f2p_set_str(*mod, "101");
    f2p_mulmod(*x, *x, *y, *mod, wp, wm);
    f2p_get_str(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_mulmod failure 5 expected 11 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *y);
    if (strcmp(buff, "1011") != 0) {
        printf("test_mulmod failure 6 expected 1011 returns %s\n", buff);
        r += 1;
    }
    if (verbose) {
        printf("end test_mulmod\n");
    }
    return r;
}

int test_powermod(int verbose, int wp, mpz_t *wm)
{
    if (verbose) {
        printf("start test_powermod\n");
    }
    int r = 0;
    mpz_t *x = &wm[wp++];
    mpz_t *e = &wm[wp++];
    mpz_t *mod = &wm[wp++];
    mpz_t *rem = &wm[wp++];
    assert(wp <= F2P_WMSIZE);
    char buff[200];

    f2p_set_str(*x, "11");
    mpz_set_ui(*e, 2);
    f2p_set_str(*mod, "111");
    f2p_powermod(*rem, *x, *e, *mod, wp, wm);
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "10") != 0) {
        printf("test_powermod failure 1 expected 10 returns %s\n", buff);
        r = 1;
    }
    f2p_get_str(buff, *mod);
    if (strcmp(buff, "111") != 0) {
        printf("test_powermod failure 2 expected 111 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_powermod failure 3 expected 11 returns %s\n", buff);
        r += 1;
    }
    int d = mpz_cmp_ui(*e, 2);
    if (d != 0) {
        printf("test_powermod failure 4 expected 0 returns %d\n", d);
        r += 1;
    }
    f2p_set_str(*x, "11");
    mpz_set_ui(*e, 3);
    f2p_set_str(*mod, "111");
    f2p_powermod(*rem, *x, *e, *mod, wp, wm);
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "1") != 0) {
        printf("test_powermod failure 5 expected 1 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *mod);
    if (strcmp(buff, "111") != 0) {
        printf("test_powermod failure 6 expected 111 returns %s\n", buff);
        r += 1;
    }
    f2p_get_str(buff, *x);
    if (strcmp(buff, "11") != 0) {
        printf("test_powermod failure 7 expected 11 returns %s\n", buff);
        r += 1;
    }
    d = mpz_cmp_ui(*e, 3);
    if (d != 0) {
        printf("test_powermod failure 8 expected 0 returns %d\n", d);
        r += 1;
    }

    f2p_set_str(*x, "11");
    mpz_set_ui(*e, 5);
    f2p_set_str(*mod, "111");
    f2p_powermod(*rem, *x, *e, *mod, wp, wm);
    f2p_get_str(buff, *rem);
    if (strcmp(buff, "10") != 0) {
        printf("test_powermod failure 10 expected 10 returns %s\n", buff);
        r += 1;
    }
    if (verbose) {
        printf("end test_powermod\n");
    }
    return r;
}

int main(int argc, char * argv[])
{
    int verbose = 0;
    int r = 0;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = 1;
    }
    int wp = 0;
    mpz_t wm[F2P_WMSIZE];

    f2p_wm_init(wm);

    r += test_string(verbose, wp, wm);
    r += test_degree(verbose, wp, wm);
    r += test_add(verbose, wp, wm);
    r += test_lshift(verbose, wp, wm);
    r += test_rshift(verbose, wp, wm);
    r += test_coefficient(verbose, wp, wm);
    r += test_mod(verbose, wp, wm);
    r += test_divrem(verbose, wp, wm);
    r += test_mulmod(verbose, wp, wm);
    r += test_powermod(verbose, wp, wm);

    f2p_wm_clear(wm);

    if (r == 0) {
        return 0;
    } else {
        return -1;
    }
}
