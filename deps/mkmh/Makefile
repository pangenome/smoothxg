#IS_ICPC:= $(shell command -v icpc 2> /dev/null)

ifdef IS_ICPC
	CXX:=icpc
	CXXFLAGS:= -O3 -std=c++11 -xAVX -qopenmp -funroll-loops -ggdb -pg
else
	CXX:=g++
	CXXFLAGS:= -O3 -std=c++11 -fopenmp -mtune=native -ggdb
endif

LD_LIB_FLAGS:= -Lmurmur3 -L. -Lxxhash
LD_INC_FLAGS:= -I. -Imurmur3 -IxxHash


example: example.cpp mkmh.hpp murmur3/libmurmur3.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lmurmur3

test: test-exe
	./test-exe

test-exe: mkmh_test.cpp mkmh.hpp murmur3/libmurmur3.a
	$(CXX) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lmurmur3

fast_test: test.cpp mkmh.hpp murmur3/libmurmur3.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lmurmur3

libmkmh.a: HASHTCounter.o mkmh.hpp murmur3/libmurmur3.a
	ar -rs $@ $^ 

HASHTCounter.o: HASHTCounter.cpp HASHTCounter.hpp
	$(CXX) $(CXXFLAGS) -c $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

murmur3/libmurmur3.a:
	+cd murmur3 && $(MAKE) lib

xxHash/libxxhash.a:
	+cd xxHash && $(MAKE)

.PHONY: clean clobber test

clean:
	$(RM) *.o
	$(RM) libmkmh.a
	cd murmur3 && $(MAKE) clean
	cd xxHash && $(MAKE) clean

clobber: clean
	$(RM) *.a
	cd murmur3 && $(MAKE) clean
	cd xxHash && $(MAKE) clean
	$(RM) test
	$(RM) example
