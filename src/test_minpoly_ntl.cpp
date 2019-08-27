#include <NTL/GF2X.h>
#include "tinymt32.h"
#include <string>
#include <stdio.h>
#include <inttypes.h>
#include <ctype.h>

using namespace NTL;
using namespace std;

typedef unsigned int (*f2rng)(void);

tinymt32_t tiny32;

unsigned int dummy()
{
    return tinymt32_generate_uint32(&tiny32);
}

void to_string(string& str, GF2X& poly)
{
    int d = deg(poly);
    int m = d / 64;
    char buff[200];
    str = "";
    for (int c = m; c >= 0; c--) {
        uint64_t mask = 1;
        uint64_t p = 0;
        for(int i = 0; i < 64; i++) {
            if (IsOne(coeff(poly, 64 * c + i))) {
                p |= mask;
            }
            mask <<= 1;
        }
        sprintf(buff, "%016" PRIx64, p);
        str = str + buff;
    }
}

void to_poly(GF2X& poly, string& str)
{
    GF2X x1;
    SetCoeff(x1, 1, 1);
    int len = str.length();
    poly = 0;
    for (int i = 0; i < len; i++) {
        char c = tolower(str[i]);
        int d;
        if (isdigit(c)) {
            d = c - '0';
        } else if (isxdigit(c)) {
            d = c - 'a' + 10;
        }
        int mask = 8;
        for (int j = 0; j < 4; j++) {
            poly *= x1;
            if (d & mask) {
                SetCoeff(poly, 0, 1);
            } else {
                SetCoeff(poly, 0, 0);
            }
            mask = mask >> 1;
        }
    }
}

void ntl_minpoly(GF2X& poly, f2rng gen, int mexp)
{
    Vec<GF2> v;
    v.SetLength(2 * mexp);
    for (int i = 0; i < 2 * mexp; i++) {
        v[i] = gen() & 1;
    }
    MinPolySeq(poly, v, mexp);
}

void ntl_minpoly(char * minpoly, f2rng gen, int mexp)
{
    GF2X poly;
    ntl_minpoly(poly, gen, mexp);
    string str;
    to_string(str, poly);
    strcpy(minpoly, str.c_str());
}

int ntl_annihilate(const GF2X& poly, f2rng gen, int mexp)
{
    int bit = 0;
    for (int i = 0; i <= mexp; i++) {
        if (IsOne(poly[i])) {
            bit ^= gen() & 1;
        }
    }
    if (bit == 0) {
        return 1;
    } else {
        return 0;
    }
}

int test_minpoly()
{
    GF2X poly;
    ntl_minpoly(poly, dummy, 127);
    int ok = 1;
    for (int i = 0; i < 100; i++) {
        int r = ntl_annihilate(poly, dummy, 127);
        ok &= r;
        dummy();
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
    //char polystr[200];
    string str;
    to_string(str, poly);
    printf("poly = %ld,%s\n", deg(poly), str.c_str());
    to_poly(poly, str);
    to_string(str, poly);
    printf("poly = %s\n", str.c_str());
    printf("minpoly end\n");
    //int deg = f2p_degree(poly);
    //printf("deg(poly) = %d\n", deg);
    return 0;
}

int main(int argc, char * argv[])
{
    int verbose = 0;
    int r = 0;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = 1;
    }
    tiny32.mat1 = 0x8f7011ee;
    tiny32.mat2 = 0xfc78ff1f;
    tiny32.tmat = 0x3793fdff;

    tinymt32_init(&tiny32, 1);
    r += test_minpoly();
    if (r == 0) {
        return 0;
    } else {
        return -1;
    }
}
