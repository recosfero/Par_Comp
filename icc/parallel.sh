#!/bin/sh

#Variables
MKL_LIBS_PARALLEL="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"

#ICC compilation with openmp (parallel processing) and MKL
icc -c -std=c99 -openmp $MKL_INC -o mmult_icc_parallel.o mmult.c

#Link with MKL Libs and build executable binary
icc -L$MKL_LIBDIR $MKL_LIBS_PARALLEL -liomp5 -lpthread -o mmult_icc_parallel mmult_icc_parallel.o -lm

rm mmult_icc_parallel.o
