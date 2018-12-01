#!/usr/bin/env python2

env = Environment (FORTRAN="gfortran",
                   F90FLAGS="-Wall -g -fcheck=all -fbacktrace -O0 -ffpe-trap=\"invalid,zero,overflow\"",
                   LINKFLAGS="",
                   # Import: We have to set the include path to the module files for the Scons scanner, without we will not trigger the scanner correctly.
                   F90PATH="#/build/mod",
                   FORTRANMODDIR="#/build/mod",
                   FORTRANMODDIRPREFIX="-J",
                   F90COMSTR="Compiling $TARGET",
                   ARCOMSTR="Archiving $TARGET",
                   LINKCOMSTR="Linking $TARGET",
                   LIBPATH=["#/build"]
                   )

SConscript(['src/SConscript'], exports='env', variant_dir='build')
SConscript(['test/SConscript'], exports='env', variant_dir='build/test')
