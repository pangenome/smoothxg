CC=g++
CFLAGS = -O3 -Wall -mtune=native

.PHONY: all clean tests lib

all: example lib

example: example.o murmur3.o
tests: test.o murmur3.o
	$(CC) $^ -o $@
	./tests

lib: libmurmur3.a


libmurmur3.a: murmur3.o
	ar -rs libmurmur3.a murmur3.o

murmur3.o: murmur3.cpp murmur3.hpp
	$(CC) $(CFLAGS) -c $<

shared: murmur3.cpp murmur3.hpp
	$(CC) -fPIC -O3 -c murmur3.cpp
	$(CC) -shared -Wl,--export-dynamic murmur3.o -o libmurmur3.so

clean:
	$(RM) *.o
	$(RM) *.so
	$(RM) *.a
