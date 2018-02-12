#IS_ICPC:= $(shell command -v icpc 2> /dev/null)

# STATIC_FLAG:= -static -static-intel


ifdef IS_ICPC
	CXX:=icpc
	#CXXFLAGS:= -O0 -std=c++11 -funroll-loops -ggdb -pg -qopenmp
	CXXFLAGS:= -O3 -xAVX -std=c++11 -qopenmp -funroll-loops
else
	CXX:=g++
	CXXFLAGS:= -O3 -std=c++11 -fopenmp -mtune=native -ggdb
endif

SRC_DIR:=src

LD_INC_FLAGS:= -Imkmh -Imkmh/murmur3 -I. -Ikseq_reader
LD_LIB_FLAGS:= -Lmkmh/murmur3 -Lmkmh -L. -lmkmh -lz -lmurmur3 -lksr

rkmh: $(SRC_DIR)/rkmh.o $(SRC_DIR)/equiv.hpp mkmh/libmkmh.a $(SRC_DIR)/HASHTCounter.o libksr.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(SRC_DIR)/HASHTCounter.o $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

$(SRC_DIR)/rkmh.o: $(SRC_DIR)/rkmh.cpp $(SRC_DIR)/equiv.hpp mkmh/libmkmh.a $(SRC_DIR)/HASHTCounter.o
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

libksr.a:
	cd kseq_reader && $(MAKE) && cp libksr.a ../ && cp *.hpp ../

mkmh/libmkmh.a: mkmh/mkmh.cpp mkmh/mkmh.hpp
	cd mkmh && $(MAKE) clean && $(MAKE)

$(SRC_DIR)/HASHTCounter.o: $(SRC_DIR)/HASHTCounter.cpp $(SRC_DIR)/HASHTCounter.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)


.PHONY: clean clobber lib static

clean:
	$(RM) $(SRC_DIR)/*.o
	cd mkmh && $(MAKE) clean
	$(RM) rkmh
