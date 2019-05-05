
# multi threads and single thread debug
debug-mt : uvc.mt.out 
debug-st : uvc.st.out
release : uvc

all : debug-mt debug-st release

HDR=CLI11-1.7.1/CLI11.hpp logging.hpp consensus.hpp CmdLineArgs.hpp 
SRC=main.cpp  bcf_formats.step1.c conversion.hpp grouping.hpp utils.hpp CmdLineArgs.cpp 

CXXFLAGS=-std=c++14
COMMIT_VERSION=$(shell git rev-parse HEAD)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" 

uvc        : $(HDR) $(SRC) instcode.hpp Makefile
	g++ -Ofast -DNDEBUG -o uvc        $(CXXFLAGS) $(SRC) $(VERFLAGS) -L/biocluster/data/bioexec/software/htslib-1.6/ -lhts -fopenmp # -l htslib

uvc.st.out : $(HDR) $(SRC) instcode.hpp Makefile
	g++ -O2 -pg         -o uvc.st.out $(CXXFLAGS) $(SRC) $(VERFLAGS) -L/biocluster/data/bioexec/software/htslib-1.6/ -lhts	

uvc.mt.out : $(HDR) $(SRC) instcode.hpp Makefile
	g++ -O2 -pg         -o uvc.mt.out $(CXXFLAGS) $(SRC) $(VERFLAGS) -L/biocluster/data/bioexec/software/htslib-1.6/ -lhts -fopenmp

bcf_formats_generator1.out : bcf_formats_generator1.cpp version.h
	g++ -o bcf_formats_generator1.out bcf_formats_generator1.cpp

bcf_formats.step1.c : bcf_formats_generator1.out
	./bcf_formats_generator1.out > bcf_formats.step1.c

#bcf_formats.o      : bcf_formats.step1.c
#	g++ -c -o bcf_formats.o bcf_formats.step1.c

.PHONY: clean dup

clean:
	rm bcf_formats_generator1.out bcf_formats.step1.c *.o *.out *.gch uvc

dup:
	cp uvc uvc1

