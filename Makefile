CXX:=g++
CXXFLAGS:= -std=c++11 -O3 -mtune=native -fopenmp
LD_INC_FLAGS:= -Imkmh -I. -Imkmh/murmur3 -Ikseq
LD_LIB_FLAGS:= -Lmkmh/murmur3 -Lmkmh -lmkmh -lz -lmurmur3


ekmmh: ekmmh.cpp equiv.hpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

#equiv.hpp:

.PHONY: clean clobber

clean:
	$(RM) *.o
	$(RM) ekmmh
