#pragma once
#ifndef DEBUG_F2P_H
#define DEBUG_F2P_H

/**
 * @file debug_f2p.h
 *
 * @brief debug functions.
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
#include <stdio.h>
#include <inttypes.h>

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(DEBUG)
static inline void f2p_print(char *s, mpz_t x)
{
    static char buff[20000];
    f2p_get_binstr(buff, x);
    int d = f2p_degree(x);
    printf("%s %d,%s\n", s, d, buff);
    fflush(stdout);
}
#define PRT(s, x) f2p_print(s, x)
#define PUTS(s) puts(s)
#else
#define PRT(s, x)
#define PUTS(s)
#endif

#if defined(__cplusplus)
}
#endif

#endif // DEBUG_F2P_H
