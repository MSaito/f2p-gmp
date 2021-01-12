#pragma once
#ifndef TEST_NTL_HPP
#define TEST_NTL_HPP
#include <NTL/GF2X.h>
#include <string>
#include <stdio.h>
#include <inttypes.h>
#include <ctype.h>

static inline void to_hexstring(std::string& str, const NTL::GF2X& poly)
{
    long d = deg(poly);
    long m = d / 64;
    char buff[200];
    str = "";
    for (long c = m; c >= 0; c--) {
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

static inline void to_binstring(std::string& str, const NTL::GF2X& poly)
{
    long d = deg(poly);
    long m = d / 64;
    str = "";
    for (int c = 0; c <= m; c++) {
        for(int i = 0; i < 64; i++) {
            if (IsOne(coeff(poly, 64 * c + i))) {
                str += '1';
            } else {
                str += '0';
            }
        }
    }
    reverse(str.begin(), str.end());
}

static inline void hexto_poly(NTL::GF2X& poly, const std::string& str)
{
    NTL::GF2X x1;
    SetCoeff(x1, 1, 1);
    unsigned long len = str.length();
    poly = 0;
    for (unsigned int i = 0; i < len; i++) {
        char c = static_cast<char>(tolower(static_cast<char>(str[i])));
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

static inline void ntl_print(const char *s, const NTL::GF2X& x)
{
    std::string str;
    to_binstring(str, x);
    long d = deg(x);
    printf("%s %ld,%s\n", s, d, str.c_str());
    fflush(stdout);
}

#if defined(DEBUG)
#define PRT(s, x) ntl_print(s, x)
#define PUTS(s) puts(s)
#else
#define PRT(s, x)
#define PUTS(s)
#endif

#endif // TEST_NTL_HPP
