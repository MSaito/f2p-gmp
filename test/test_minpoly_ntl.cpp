#include "test_ntl.hpp"
#include <NTL/GF2X.h>
#define LINEARITY_CHECK
#include "tinymt32.h"
#include <string>
#include <stdio.h>
#include <inttypes.h>
#include <ctype.h>

using namespace NTL;
using namespace std;

tinymt32_t tiny32;

unsigned int dummy()
{
    return tinymt32_generate_uint32(&tiny32);
}

void ntl_exeuclid2(GF2X& a, GF2X& c, GF2X& x, GF2X& y, int m)
{
    PUTS("f2p_exeuclid2 start\n");
    assert(!IsZero(x));
    assert(!IsZero(y));
    GF2X q1;
    GF2X r0(x);
    GF2X r1(y);
    GF2X r2;
    GF2X a0;
    GF2X a1;
    GF2X a2;
    GF2X tmp;
    PRT("x= ", x);
    PRT("y= ", y);
    a0 = 1;
    a1 = 0;
    PRT("r0 = ", r0);
    PRT("r1 = ", r1);
    PRT("a0 = ", a0);
    PRT("a1 = ", a1);
    int dr = deg(r0);
    while (dr >= m) {
        DivRem(q1, r2, r0, r1);
        tmp = q1 * a1;
        a2 = a0 + tmp;
        PRT("q1 = ", q1);
        PRT("tmp = ", tmp);
        PRT("r0 = ", r0);
        PRT("r1 = ", r1);
        PRT("r2 = ", r2);
        PRT("a0 = ", a0);
        PRT("a1 = ", a1);
        PRT("a2 = ", a2);
        r0 = r1;
        r1 = r2;
        a0 = a1;
        a1 = a2;
        dr = deg(r0);
        PUTS("loop last\n");
    }
    PUTS("loop end\n");
    PRT("a0 = ", a0);
    a = a0;
    PRT("r0 = ", r0);
    c = r0;
    PUTS("f2p_exeuclid2 end\n");
}

#if 0
void ntl_minpoly(GF2X& poly, f2rng gen, int mexp)
{
    Vec<GF2> v;
    v.SetLength(2 * mexp);
    for (int i = 0; i < 2 * mexp; i++) {
        v[i] = gen() & 1;
    }
    MinPolySeq(poly, v, mexp);
}
#endif

void ntl_minpoly(GF2X& poly, tinymt32_t * tiny32, int mexp)
{
    GF2X v;
    v.SetLength(2 * mexp);
    //for (int i = 0; i < 2 * mexp; i++) {
    for (int i = 2 * mexp - 1; i >= 0; i--) {
        v[i] = tinymt32_generate_uint32(tiny32) & 1;
    }
    PRT("seq = ", v);
    //MinPolySeq(poly, v, mexp);
    GF2X x2m;
    GF2X c;
    SetCoeff(x2m, 2 * mexp, 1);
    ntl_exeuclid2(poly, c, v, x2m, mexp);
    printf("deg(c) = %ld\n", deg(c));
}

void ntl_minpoly(char * minpoly, tinymt32_t * tiny32, int mexp)
{
    GF2X poly;
    ntl_minpoly(poly, tiny32, mexp);
    string str;
    to_hexstring(str, poly);
    strcpy(minpoly, str.c_str());
}

int ntl_annihilate(const GF2X& poly, tinymt32_t * tiny32, int mexp)
{
    int bit = 0;
    for (int i = 0; i <= mexp; i++) {
        if (IsOne(poly[i])) {
            bit ^= tinymt32_generate_uint32(tiny32) & 1;
        }
    }
    if (bit == 0) {
        return 1;
    } else {
        return 0;
    }
}

int test_minpoly(tinymt32_t * tiny32)
{
    GF2X poly;
    ntl_minpoly(poly, tiny32, 127);
    int ok = 1;
    for (int i = 0; i < 100; i++) {
        int r = ntl_annihilate(poly, tiny32, 127);
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
    to_hexstring(str, poly);
    printf("poly = %ld,%s\n", deg(poly), str.c_str());
    hexto_poly(poly, str);
    to_hexstring(str, poly);
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
    r += test_minpoly(&tiny32);
    if (r == 0) {
        return 0;
    } else {
        return -1;
    }
}
