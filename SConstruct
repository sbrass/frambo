#!/usr/bin/env python2

AddOption('--prefix',
          default='/usr/local',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix')

env = Environment (PREFIX=GetOption('prefix'),
                   FORTRAN="gfortran",
                   F90FLAGS="-Wall -g -fcheck=all -fbacktrace -O0 -ffpe-trap=\"invalid,zero,overflow\"",
                   LINKFLAGS="",
                   CPPPATH="#/include",
                   # Import: We have to set the include path to the module files for the Scons scanner, without we will not trigger the scanner correctly.
                   F90PATH="#/build/mod",
                   FORTRANMODDIR="#/build/mod",
                   FORTRANMODDIRPREFIX="-J",
                   FORTRANCOMSTR="Compiling $TARGET",
                   F90COMSTR="Compiling $TARGET",
                   ARCOMSTR="Archiving $TARGET",
                   # LINKCOMSTR="Linking $TARGET",
                   LIBPATH=["#/build"]
                   )

env.Alias("install", ["$PREFIX/lib", "$PREFIX/include"])

SConscript(['include/SConscript'], exports='env')
SConscript(['src/SConscript'], exports='env', variant_dir='build')
SConscript(['test/SConscript'], exports='env', variant_dir='build/test')
