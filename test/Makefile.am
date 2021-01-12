AUTOMAKE_OPTIONS = subdir-objects

TESTS = test_f2p test_f2p_exeuclid test_minpoly test_minpoly_ntl \
 test_irreducible_ntl test_jump_ntl

check_PROGRAMS = test_f2p test_f2p_exeuclid test_minpoly test_minpoly_ntl \
 test_irreducible_ntl test_jump_ntl

EXTRA_DIST = tinymt32.c tinymt32.h test_ntl.hpp

#noinst_PROGRAMS = test_f2p test_f2p_exeuclid test_minpoly test_minpoly_ntl \
#test_irreducible_ntl test_jump_ntl

test_f2p_SOURCES = test_f2p.c ../src/f2p_gmp.c
test_f2p_exeuclid_SOURCES = test_f2p_exeuclid.c ../src/f2p_gmp.c
test_minpoly_SOURCES = test_minpoly.c ../src/f2p_gmp.c tinymt32.c
test_minpoly_ntl_SOURCES = test_minpoly_ntl.cpp tinymt32.c
test_irreducible_ntl_SOURCES = test_irreducible_ntl.cpp ../src/f2p_gmp.c
test_jump_ntl_SOURCES = test_jump_ntl.cpp ../src/f2p_gmp.c tinymt32.c

AM_CXXFLAGS = -Wall -O2 -Wextra -Wsign-compare -Wconversion \
	      -D__STDC_CONSTANT_MACROS \
	      -D__STDC_FORMAT_MACROS -I../src #-DDEBUG

AM_CFLAGS = -std=c99 -Wall -Wextra -Wsign-compare -Wconversion -I../src #-DDEBUG
