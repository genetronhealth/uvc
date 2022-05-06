# Example command to build: make all -j9 && make deploy

release  = debarcode uvc-1-fopenmp-thread             uvc-1-cpp-std-thread 
release  : $(release)
debug    = debarcode uvc-2-fopenmp-thread-asan.debug  uvc-2-cpp-std-thread-asan.debug  uvc-3-asan.debug  uvc-4.debug 
debug    : $(debug)
debug-ub : debarcode uvc-2-fopenmp-thread-ubsan.debug uvc-2-cpp-std-thread-ubsan.debug uvc-3-ubsan.debug 
all      : release debug
ALL      : all     debug-ub

HDR=CLI11-1.7.1/CLI11.hpp Hash.hpp main_conversion.hpp main_consensus.hpp \
    CmdLineArgs.hpp common.hpp grouping.hpp iohts.hpp logging.hpp main.hpp MolecularID.hpp version.h
SRC=CmdLineArgs.cpp common.cpp grouping.cpp iohts.cpp logging.cpp main.cpp MolecularID.cpp version.cpp 
DEP=bcf_formats.step1.hpp instcode.hpp Makefile

HTSPATH=ext/htslib-1.11-lowdep/libhts.a
HTSFLAGS=$(HTSPATH) -I ext/htslib-1.11-lowdep/ -pthread -lm -lz -lbz2 -llzma # -lcurl -lcrypto # can be changed depending on the specific installed components of htslib (please refer to the INSTALL file in htslib)
CC=gcc  # can be changed to clang or other compilers as needed
CXX=g++ # can be changed to clang or other compilers as needed
CXXFLAGS=-std=c++14 -static-libstdc++ -Wall 
COMMIT_VERSION=$(shell git rev-parse HEAD | head -c 7)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
COMMIT_DIFF_FULL=$(shell echo "R\"ZXF_specQUOTE(\n $$(git diff HEAD | sed 's/ZXF_specQUOTE/ZXF_specquote/g') \n)ZXF_specQUOTE\"" > gitdiff.txt)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" -DCOMMIT_DIFF_FULL="\"$(COMMIT_DIFF_FULL)\""
# UVC_IN_DEBUG_MODE enables locus-specific diagnostic debugging info
DEBUG_OPTS=-DUVC_IN_DEBUG_MODE -static-libasan 
UBSAN=--param=max-vartrack-size=640000000 -fsanitize=undefine 

debarcode  : debarcode_main.c version.h Makefile
	$(CC) -O3 -o debarcode $(VERFLAGS) debarcode_main.c ${HTSFLAGS}
	
# the main executable, uses OpenMP for multi-threading
uvc-1-fopenmp-thread            : $(HDR) $(SRC) $(DEP)
	$(CXX) -O3 -DNDEBUG -o uvc-1-fopenmp-thread             $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) -fopenmp # -l htslib

# the main executable, uses C++ standard template library thread for multi-threading, useful if OpenMP runtime is not available
uvc-1-cpp-std-thread            : $(HDR) $(SRC) $(DEP)
	$(CXX) -O3 -DNDEBUG -o uvc-1-cpp-std-thread             $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) -DUSE_STDLIB_THREAD # -l htslib

uvc-2-fopenmp-thread-asan.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O2 -g -p    -o uvc-2-fopenmp-thread-asan.debug  $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) -fopenmp -fsanitize=address
uvc-2-fopenmp-thread-ubsan.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O2 -g -p    -o uvc-2-fopenmp-thread-ubsan.debug $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) -fopenmp $(UBSAN)
uvc-2-cpp-std-thread-asan.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O2 -g -p    -o uvc-2-cpp-std-thread-asan.debug  $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) -DUSE_STDLIB_THREAD -fsanitize=address
uvc-2-cpp-std-thread-ubsan.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O2 -g -p    -o uvc-2-cpp-std-thread-ubsan.debug $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) -DUSE_STDLIB_THREAD $(UBSAN)
uvc-3-asan.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O2 -g -p    -o uvc-3-asan.debug                 $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) -fsanitize=address
uvc-3-ubsan.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O2 -g -p    -o uvc-3-ubsan.debug                $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) $(UBSAN)
uvc-4.debug : $(HDR) $(SRC) $(DEP)
	$(CXX) -O0 -g -p    -o uvc-4.debug                      $(CXXFLAGS) $(VERFLAGS) $(SRC) $(HTSFLAGS) $(DEBUG_OPTS) -Wextra -DENABLE_ASSERT_IN_UVC

bcf_formats_generator1.out : bcf_formats_generator1.cpp version.h 
	$(CXX) -o bcf_formats_generator1.out $(CXXFLAGS) bcf_formats_generator1.cpp

bcf_formats.step1.hpp : bcf_formats_generator1.out
	./bcf_formats_generator1.out > bcf_formats.step1.hpp # auto-generate the C++ code from the BCF-template generator

.PHONY: release all debug debug-ub ALL clean deploy

clean:
	rm bcf_formats_generator1.out bcf_formats.step1.hpp *.o *.debug uvc-1-fopenmp-thread uvc-1-cpp-std-thread *.gch debarcode || true
	
deploy:
	cp uvc-1-fopenmp-thread bin/uvc1 # The default binary executable uses OpenML thread, and uvc1 is used by uvcTN.sh
	cp $(release) bin/ || true
	cp $(debug)   bin/ || true

