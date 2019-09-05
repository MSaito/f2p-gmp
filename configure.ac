#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([f2p-gmp], [0.1], [sai10@hiroshima-u.ac.jp])
#LT_INIT
AC_CONFIG_SRCDIR([src/Makefile.am])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE

AC_CANONICAL_HOST

# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
#AC_PROG_LIBTOOL

# keep this order
AC_CHECK_LIB(gmp, __gmpz_init)
AC_CHECK_LIB(gf2x, gf2x_mul)
AC_CHECK_LIB(ntl, main)

AC_LANG_PUSH(C++)
AX_CXX_COMPILE_STDCXX_0X

AX_CXX_COMPILE_STDCXX_14(noext, optional) # keep this order
AX_CXX_COMPILE_STDCXX_11(noext, optional)
AC_LANG_POP

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h stdint.h stdlib.h])
#AX_EXT

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
