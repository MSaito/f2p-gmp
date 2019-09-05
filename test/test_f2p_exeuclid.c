#include "f2p_gmp.h"
#include <string.h>

int test_exeuclid(int verbose, int wp, f2p_wm_t *wm)
{
    if (verbose) {
        printf("start test_euclid\n");
    }
    int ok = 1;
    mpz_t *x = &(wm->ar[wp++]);
    mpz_t *y = &(wm->ar[wp++]);
    mpz_t *a = &(wm->ar[wp++]);
    mpz_t *b = &(wm->ar[wp++]);
    mpz_t *c = &(wm->ar[wp++]);
    assert(wp <= wm->max_size);
    char buff[200];

    f2p_set_binstr(*x, "110");
    f2p_set_binstr(*y, "1");
    f2p_exeuclid(*a, *b, *c, *x, *y, wp, wm);
    f2p_get_binstr(buff, *c);
    if (strcmp(buff, "1") != 0) {
        printf("failure 0.1 expected 1 returns %s\n", buff);
        ok = 0;
    }

    f2p_set_binstr(*x, "110");
    f2p_set_binstr(*y, "111");
    f2p_exeuclid(*a, *b, *c, *x, *y, wp, wm);
    f2p_get_binstr(buff, *x);
    if (strcmp(buff, "110") != 0) {
        printf("failure 1 expected 110 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *y);
    if (strcmp(buff, "111") != 0) {
        printf("failure 2 expected 111 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *a);
    if (strcmp(buff, "1") != 0) {
        printf("failure 3 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *b);
    if (strcmp(buff, "1") != 0) {
        printf("failure 4 expected 1 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *c);
    if (strcmp(buff, "1") != 0) {
        printf("failure 5 expected 1 returns %s\n", buff);
        ok = 0;
    }

    f2p_set_binstr(*x, "1101");
    f2p_set_binstr(*y, "11111");
    f2p_exeuclid(*a, *b, *c, *x, *y, wp, wm);
    f2p_get_binstr(buff, *a);
    if (strcmp(buff, "1100") != 0) {
        printf("failure 6 expected 1100 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *b);
    if (strcmp(buff, "111") != 0) {
        printf("failure 7 expected 111 returns %s\n", buff);
        ok = 0;
    }
    f2p_get_binstr(buff, *c);
    if (strcmp(buff, "1") != 0) {
        printf("failure 8 expected 1 returns %s\n", buff);
        ok = 0;
    }
    if (verbose) {
        printf("end test_euclid\n");
    }
    if (ok) {
        printf("o");
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

    f2p_wm_init(&wm, 30);

    ok *= test_exeuclid(verbose, wp, &wm);
    printf("\n");
    f2p_wm_clear(&wm);

    if (ok == 1) {
        return 0;
    } else {
        return -1;
    }
}
