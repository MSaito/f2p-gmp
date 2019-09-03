/**
 * @file test_irreducible_ntl.hpp
 *
 * @test program for f2p_is_irreducible.
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

#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include "test_ntl.hpp"
#include "f2p_gmp.h"
#include <string>
#include <stdio.h>
#include <inttypes.h>

using namespace NTL;
using namespace std;

int test_irreducible(const string& hexstr)
{
    GF2X poly;
    hexto_poly(poly, hexstr);
    bool irre_ntl = IterIrredTest(poly);
    mpz_t gmppoly;
    mpz_init(gmppoly);
    f2p_set_hexstr(gmppoly, hexstr.c_str());
    bool irre_f2p = f2p_is_irreducible(gmppoly);
    int ok = 1;
    if (irre_ntl == irre_f2p) {
        ok = 1;
    } else {
        printf("irrentl = %d, irre_f2p = %d\n", irre_ntl, irre_f2p);
        ok = 0;
    }
    mpz_clear(gmppoly);
    return ok;
}

int test_irreducible(int n, const bool verbose)
{
    GF2X poly;
    GF2X ranpoly;
    int ok = 1;
    string str;
    char okstr[2] = {'x','o'};
    BuildIrred(poly, n);
    for (int i = 0; i < 100; i++) {
        BuildRandomIrred(ranpoly, poly);
        to_string(str, ranpoly);
        int r = test_irreducible(str);
        if (verbose) {
            printf("%c", okstr[r & 1]);
        }
        ok &= r;
    }
    for (int i = 0; i < 100; i++) {
        random(ranpoly, n);
        to_string(str, ranpoly);
        int r = test_irreducible(str);
        if (verbose) {
            printf("%c", okstr[r & 1]);
        }
        ok &= r;
    }
    if (verbose) {
        printf("\n");
    }
    if (ok) {
        printf("irreducible OK\n");
        return 0;
    } else {
        printf("irreducible NG\n");
        return 1;
    }
}

int main(int argc, char * argv[])
{
    bool verbose = false;
    int r = 0;
    int n = 200;
    if (argc > 1 && argv[1][0] == 'v') {
        verbose = true;
    }
    r += test_irreducible(n, verbose);
    if (r == 0) {
        return 0;
    } else {
        return -1;
    }
}
