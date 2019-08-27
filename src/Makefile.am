TESTS = test_f2p test_f2p_exeuclid test_minpoly

noinst_PROGRAMS = test_f2p test_f2p_exeuclid test_minpoly test_minpoly_ntl

test_f2p_SOURCES = test_f2p.c f2p-gmp.c
test_f2p_exeuclid_SOURCES = test_f2p_exeuclid.c f2p-gmp.c
test_minpoly_SOURCES = test_minpoly.c f2p-gmp.c tinymt32.c
test_minpoly_ntl_SOURCES = test_minpoly_ntl.cpp tinymt32.c

test_minpoly_ntl_CXXFLAGS = -Wall -O2 -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS #-lntl

AM_CFLAGS = -std=c99 -Wall -Wextra #-DDEBUG
