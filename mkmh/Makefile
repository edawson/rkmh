CXX:=g++
CXXFLAGS:= -O3 -mtune=native -std=c++11
LD_LIB_FLAGS:= -Lmurmur3 -L.
LD_INC_FLAGS:= -I. -Imurmur3

libmkmh.a: mkmh.o murmur3/libmurmur3.a Makefile
	ar -rs $@ $< murmur3/libmurmur3.a

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lmkmh -lmurmur3

test: test.cpp libmkmh.a murmur3/libmurmur3.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lmkmh -lmurmur3
	./test

mkmh.o: mkmh.cpp mkmh.hpp murmur3/libmurmur3.a murmur3/murmur3.hpp
	$(CXX) $(CXXFLAGS) -c $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lmurmur3

murmur3/libmurmur3.a:
	+cd murmur3 && $(MAKE) lib

.PHONY: clean clobber

clean:
	$(RM) *.o

clobber: clean
	$(RM) libmkmh
	$(RM) *.a
	cd murmur3 && $(MAKE) clean
	$(RM) test
	$(RM) example
