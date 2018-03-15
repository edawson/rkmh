#IS_ICPC:= $(shell command -v icpc 2> /dev/null)

# STATIC_FLAG:= -static -static-intel


ifdef IS_ICPC
	CXX:=icpc
	#CXXFLAGS:= -O0 -std=c++11 -funroll-loops -ggdb -pg -qopenmp
	CXXFLAGS:= -O3 -xAVX -std=c++11 -qopenmp -funroll-loops
else
	CXX:=g++
	CXXFLAGS:= -O3 -std=c++11 -fopenmp -mtune=native -ggdb -g
endif

SRC_DIR:=src

LD_INC_FLAGS:= -Imkmh -Imkmh/murmur3 -I. -Ikseq_reader
LD_LIB_FLAGS:= -Lmkmh/murmur3 -Lmkmh -L. -Lkseq_reader -lmkmh -lz -lmurmur3 -lksr

rkmh: $(SRC_DIR)/rkmh.o $(SRC_DIR)/equiv.hpp mkmh/libmkmh.a kseq_reader/libksr.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

$(SRC_DIR)/rkmh.o: $(SRC_DIR)/rkmh.cpp $(SRC_DIR)/equiv.hpp mkmh/libmkmh.a
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

kseq_reader/libksr.a: kseq_reader/kseq_reader.cpp kseq_reader/kseq_reader.hpp
	cd kseq_reader && $(MAKE)

mkmh/libmkmh.a: mkmh/mkmh.cpp mkmh/mkmh.hpp
	cd mkmh && $(MAKE)

.PHONY: clean clobber lib static

clean:
	$(RM) $(SRC_DIR)/*.o
	cd mkmh && $(MAKE) clean
	cd kseq_reader && $(MAKE) clean
	$(RM) rkmh
