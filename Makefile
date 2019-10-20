
# Example command to build: make all -j9 && make dup

# Multi threads and single thread debug
release : uvc debarcode
debug-mt : uvc.mt.out
debug-st : uvc.st.out
test-cppt : uvc.cppt.out
all : debug-mt debug-st release test-cppt

HDR=CLI11-1.7.1/CLI11.hpp logging.hpp consensus.hpp CmdLineArgs.hpp 
SRC=main.cpp  bcf_formats.step1.c conversion.hpp grouping.hpp grouping.cpp utils.hpp CmdLineArgs.cpp

HTSFLAGS=ext/lib/libhts.a -lm -lz -lcurl -lbz2 # -llzma -lcrypto # can be changed depending on the specific installed components of htslib (please refer to the INSTALL file in htslib)
CC=gcc  # can be changed to clang or other compilers as needed
CXX=g++ # can be changed to clang or other compilers as needed
CXXFLAGS=-std=c++14 -static-libstdc++
COMMIT_VERSION=$(shell git rev-parse HEAD)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" 

debarcode  : debarcode_main.c version.h
	$(CC) -O3 -o debarcode debarcode_main.c ext/lib/libhts.a -lz
	
# the main executable, uses OpenMP for multi-threading
uvc        : $(HDR) $(SRC) instcode.hpp Makefile
	$(CXX) -O3 -DNDEBUG -o uvc          $(CXXFLAGS) $(VERFLAGS) main.cpp grouping.cpp CmdLineArgs.cpp logging.cpp $(HTSFLAGS) -fopenmp # -l htslib

# the main executable, use C++ standard template library thread for multi-threading, useful if OpenMP runtime is not available
uvc.cppt.out : $(HDR) $(SRC) instcode.hpp Makefile
	$(CXX) -O3 -DNDEBUG -o uvc.cppt.out $(CXXFLAGS) $(VERFLAGS) main.cpp grouping.cpp CmdLineArgs.cpp logging.cpp $(HTSFLAGS) -DUSE_STDLIB_THREAD # -l htslib

# single-thread executable with runtime assertions and debug symbols, very useful for debugging
uvc.st.out : $(HDR) $(SRC) instcode.hpp Makefile
	$(CXX) -O2 -g -p    -o uvc.st.out   $(CXXFLAGS) $(VERFLAGS) main.cpp grouping.cpp CmdLineArgs.cpp logging.cpp $(HTSFLAGS) 

# multi-thread executable with runtime assertions and debug symbols, useful for debugging
uvc.mt.out : $(HDR) $(SRC) instcode.hpp Makefile
	$(CXX) -O2 -g -p    -o uvc.mt.out   $(CXXFLAGS) $(VERFLAGS) main.cpp grouping.cpp CmdLineArgs.cpp logging.cpp $(HTSFLAGS) -fopenmp

# generator for bcf templates
bcf_formats_generator1.out : bcf_formats_generator1.cpp version.h
	$(CXX) -o bcf_formats_generator1.out bcf_formats_generator1.cpp

# step1.c is the template code auto-generated by the generator for bcf templates
bcf_formats.step1.c : bcf_formats_generator1.out
	./bcf_formats_generator1.out > bcf_formats.step1.c

.PHONY: clean dup

clean:
	rm bcf_formats_generator1.out bcf_formats.step1.c *.o *.out *.gch uvc debarcode || true

# uvc1 is used by uvcTN.sh
dup:
	cp uvc uvc1

