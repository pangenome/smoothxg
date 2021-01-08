#include "mkmh.hpp"
#include <string>
#include <iostream>

using namespace std;
using namespace mkmh;

int main(int argc, char **argv) {
    int hashSize = 1000;
    int kSize = 10;

    if (argc < 2) {
        cerr << "Usage:" << endl
             << "example <fileToHash> [<filetoCompare>]" << endl;
        exit(1);
    } else if (argc == 2) {
        cerr << "Hashing " << argv[1] << "..." << endl;
        string input = string(argv[1]);

    } else if (argc == 3) {
        cerr << "Hashing " << argv[1] << " and " << argv[2] <<
             " and calculating their similarity..." << endl;
    } else {
        cerr << "Usage:" << endl
             << "example <fileToHash> [<filetoCompare>]" << endl;
        exit(1);

    }


    return 0;
}
