#!/bin/bash

gcc -std=c99 -D BUILD_EXE test.c hilbert.c -I`pwd` -o test
