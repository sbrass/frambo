# Fortran RAMBO #

We have implemented the "RAMBO on diet" by S. Plätzer[^1] in modern Fortran 2003.

For details on the algorithm, see the manual or the original paper by S. Plätzer.

In addition to the Fortran API, we have implemented a C interface.

## Compilation ##

Prerequisite: scons

TODO Deliver frambo with scons.

Install SCons with `pip`:
```
pip install --user scons
```
Then, run `scons`  inside the main directory.

Recommended: Run `scons test` afterwards, all test shall pass, before using `frambo`.

## Installation ##

All files are installed by:
```
scons [--prefix=<user dir>] install
```
The default prefix is `/usr/local`.

## Deinstallation ##

All files including the installed files can be removed with:
```
scons -c install
```

[^1]: PLÄTZER, Simon. RAMBO on diet. arXiv preprint arXiv:1308.2922, 2013.
