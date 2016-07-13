IS_ICPC:= $(shell command -v icpc 2> /dev/null)

ifdef IS_ICPC
	CXX:=icpc
	CXXFLAGS:= -O3 -std=c++11 -xAVX -qopenmp -funroll-loops -ggdb
else
	CXX:=g++
	CXXFLAGS:= -O3 -std=c++11 -fopenmp -mtune=native -ggdb
endif

LD_INC_FLAGS:= -Imkmh -I. -Imkmh/murmur3 -Ikseq
LD_LIB_FLAGS:= -Lmkmh/murmur3 -Lmkmh -lmkmh -lz -lmurmur3


rkmh: rkmh.cpp equiv.hpp mkmh/libmkmh.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

mkmh/libmkmh.a: mkmh/mkmh.cpp mkmh/mkmh.hpp
	cd mkmh && $(MAKE) clean && $(MAKE)


.PHONY: clean clobber lib

clean:
	$(RM) *.o
	$(RM) rkmh
