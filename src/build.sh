#!/bin/bash
gcc -O3 -march=native -mavx512f -mavx main.c  -lraylib -lm -o ../game
