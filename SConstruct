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
                   # Gfortran correctly include all flags.
                   CC="gfortran",
                   FORTRAN="gfortran",
                   F90FLAGS="-Wall -O2",
                   LINKFLAGS="",
                   CPPPATH="#/include",
                   # Import: We have to set the include path to the module files for the Scons scanner, without we will not trigger the scanner correctly.
                   F90PATH="#/build/mod",
                   FORTRANMODDIR="#/build/mod",
                   FORTRANMODDIRPREFIX="-J",
                   FORTRANCOMSTR="Compiling $TARGET",
                   CCCOMSTR="Compiling $TARGET",
                   F90COMSTR="Compiling $TARGET",
                   ARCOMSTR="Archiving $TARGET",
                   LINKCOMSTR="Linking $TARGET",
                   LINK='gfortran',
                   LIBPATH=["#/build"]
                   )

env.Alias("install", ["$PREFIX/lib", "$PREFIX/include"])

SConscript(['include/SConscript'], exports='env')
SConscript(['src/SConscript'], exports='env', variant_dir='build')
SConscript(['test/SConscript'], exports='env', variant_dir='build/test')
