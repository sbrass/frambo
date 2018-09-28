#!/usr/bin/env python2

env = Environment (FORTRAN="gfortran",
                   F08="gfortran",
                   F08FLAGS="-Wall -g -fcheck=all -fbacktrace -O3 -ffast-math -ffpe-trap=\"invalid,zero,overflow\"",
                   LINKFLAGS="",
                   FORTRANMODDIR="mod",
                   FORTRANMODDIRPREFIX="-J",
                   F08PATH="mod",
                   F08COMSTR="Compiling $TARGET",
                   ARCOMSTR = "Archiving $TARGET",
                   LINKCOMSTR="Linking $TARGET",
                   LIBPATH=["."]
                   )

env.Library("rambo", Glob("*.f08"))

env.Program("test", "test.f08", LIBS=["rambo"])
