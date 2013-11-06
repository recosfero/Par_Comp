all : mmult

MKLSEQUENTIAL=-L$(MKL_LIBDIR) $(MKL_INC) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

mmult : mmult.c
	gcc -o mmult mmult.c $(MKLSEQUENTIAL)
