#!/bin/bash

gcc -std=c99 --shared -fPIC -O3 test.c hilbert.c -I`pwd` -o hilbert.so

