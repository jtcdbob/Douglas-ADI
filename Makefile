# Automatically switching between mac and Totient
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
		PLATFORM=mac
else
		PLATFORM=icc
endif

include Makefile.in.$(PLATFORM)

# Main function
all:
	$(CC) $(CFLAGS) adi_main.cpp -o douglas-adi

douglas-adi.dSYM: douglas-adi
	dsymutil douglas-adi -o douglas-adi.dSYM

.PHONY: iprofile
iprofile: douglas-adi douglas-adi.dSYM
	iprofiler -timeprofiler -o douglas-adi_perf ./douglas-adi
	open douglas-adi_perf.dtps

# Clean up
.PHONY: clean
clean:
	rm -f douglas-adi *.o
	rm -f *.optrpt
	rm -rf *.dSYM
	rm -rf douglas-adi_perf.dtps
