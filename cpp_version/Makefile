# Automatically switching between mac and Totient
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
		PLATFORM=mac
else
		PLATFORM=icc
endif

include Makefile.in.$(PLATFORM)

# ---
# Main function
all:
	$(CC) $(OMP_CFLAGS) adi_main.cpp -o douglas-adi

douglas-adi.dSYM: douglas-adi
	dsymutil douglas-adi -o douglas-adi.dSYM

# ---
# Rules for testing runs
.PHONY: run run-local
run:
	qsub -l nodes=1:ppn=24 qrun.pbs

run-local: 
	./douglas-adi -n 40 -p 2

# ---
# Rules for profiling using iprofiler
.PHONY: iprofile
iprofile: douglas-adi douglas-adi.dSYM
	iprofiler -timeprofiler -o douglas-adi_perf ./douglas-adi
	open douglas-adi_perf.dtps

# ---
# Rules for Cleaning up
.PHONY: clean realclean
clean:
	rm -f douglas-adi *.o
	rm -f *.optrpt
	rm -rf *.dSYM
	rm -rf douglas-adi_perf.dtps
realclean: clean
	rm Z_douglas-adi*
