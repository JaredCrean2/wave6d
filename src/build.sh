#!/bin/bash

foptflags='-Ofast'
gfortran -fdefault-real-8 $foptflags test_inline.f90 -o test_inline_f
