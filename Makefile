CXX:=g++
CXXFLAGS:= -std=c++11 -O3 -mtune=native
LD_INC_FLAGS:= -Imkmh -I. -Imkmh/murmur3 -Ikseq
LD_LIB_FLAGS:= -Lmkmh -lmkmh -lz


ekmmh: ekmmh.cpp equiv.hpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

#equiv.hpp:

.PHONY: clean clobber

clean:
	$(RM) *.o
	$(RM) ekmmmh
