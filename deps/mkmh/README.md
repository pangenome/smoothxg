# mkmh
Make kmers, minimizers, hashes, and MinHash sketches (with multiple k), and compare them. 

![C/C++ CI for mkmh](https://github.com/edawson/mkmh/workflows/C/C++%20CI/badge.svg)

## Usage
To use mkmh functions in your code:  
1. Include the header file in your code  
    ```#include "mkmh.hpp"```      
2. Compile the library:  
    `` cd mkmh && make lib``  
3. Make sure the lib and header are on the LD include/lib paths (e.g. in your makefile):  
    `` gcc -o my_code my_code.cpp -L/path/to/mkmh -I/path/to/mkmh -lmkmh  
4. That's it!

## Available functionality
Convenience functions:  
    - Reverse complement a string  
    - Reverse a string  
    - Capitalize the characters of a string
    - Check if a string contains only canonical DNA letters ("A", "a", "C", "c", "T", "t", "G", "g")


Substrings and transforms:  
    - Get the forward shingles of a string  
    - Get the kmers size *k* of a string  
    - For multiple *k*, Get the kmers of a string for all *k*  
    - Get the (*w*, *k*) minimizers of a string  
    - Calculate the 64-bit hashes of the kmers of a string (with either single or multiple *k* values)  
    - Get the MinHash sketch of a string (from either single or multiple *k* values), using either the top *s* hashes or the bottom *s* hashes.  


Compare sets of shingles / kmers / minimizers / hashes:  
    - Take the union of two sets of kmers or hashes.  
    - Take the intersection of two sets of kmers or hashes.


Fun extras:  
    - Given a string and a set of query strings, sort the queries in order
    of percent similarity.

## Getting help
Please reach out through [github](https://github.com/edawson/mkmh) by posting an issue (even if it's just feedback). Email is acceptable as a secondary medium.
