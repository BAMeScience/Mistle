#ifndef SIMPLE_EXAMPLE_MGF_READER_H
#define SIMPLE_EXAMPLE_MGF_READER_H
#include "spectrum.h"

using namespace std;

class mgf_reader {
public:
    mgf_reader();

    static vector<spectrum*> read_file(string path);
};


#endif //SIMPLE_EXAMPLE_MGF_READER_H
