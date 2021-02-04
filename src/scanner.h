#ifndef SIMPLE_EXAMPLE_SCANNER_H
#define SIMPLE_EXAMPLE_SCANNER_H

#include <string>

using namespace std;

class scanner {

public:
    scanner();
    bool scan_directory(string path);
    bool scan_file(string path);
};


#endif //SIMPLE_EXAMPLE_SCANNER_H
