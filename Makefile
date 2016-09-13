IS_ICPC:= $(shell command -v icpc 2> /dev/null)

# STATIC_FLAG:= -static -static-intel


ifdef IS_ICPC
	CXX:=icpc
	CXXFLAGS:= -O3 -std=c++11 -xAVX -qopenmp -funroll-loops -ggdb -pg
else
	CXX:=g++
	CXXFLAGS:= -O3 -std=c++11 -fopenmp -mtune=native -ggdb
endif

LD_INC_FLAGS:= -Imkmh -Imkmh/murmur3 -I. #-Ikseq
LD_LIB_FLAGS:= -Lmkmh/murmur3 -Lmkmh -lmkmh -lz -lmurmur3

rkmh: rkmh.o equiv.hpp mkmh/libmkmh.a HASHTCounter.o
	$(CXX) $(CXXFLAGS) -o $@ $< HASHTCounter.o $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

rkmh.o: rkmh.cpp equiv.hpp mkmh/libmkmh.a HASHTCounter.o
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

static: rkmh.cpp equiv.hpp mkmh/libmkmh.a
	$(CXX) $(CXXFLAGS) $(STATIC_FLAG) -o rkmh $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

mkmh/libmkmh.a: mkmh/mkmh.cpp mkmh/mkmh.hpp
	cd mkmh && $(MAKE) clean && $(MAKE)

HASHTCounter.o: HASHTCounter.cpp HASHTCounter.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)


.PHONY: clean clobber lib static

clean:
	$(RM) *.o
	$(RM) rkmh
