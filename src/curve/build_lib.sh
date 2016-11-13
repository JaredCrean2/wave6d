#!/bin/bash

gcc -std=c99 --shared -fPIC test.c hilbert.c -I`pwd` -o hilbert.so

