#!/bin/bash

#foptflags='-Ofast'
#gfortran -fdefault-real-8 $foptflags test_inline.f90 -o test_inline_f

# -fast causes link time warnings, but the executable ran
ifort -fpp -O3 -xHost test_inline.f90 -o test_inline_f
