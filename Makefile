CC=g++
CFLAGS=-g3 -DSYSTIME
# Compiler Options
OPTFLAGS=-Ofast -march=native -ffast-math -O3 -fopenmp -D_OMP
ANALYSIS=-Rpass-analysis=loop-vectorize

CFLAGS+=$(OPTFLAGS) $(ANALYSIS)

all:
	$(CC) $(CFLAGS) adi_main.cpp -o DOUGLAS_ADI
