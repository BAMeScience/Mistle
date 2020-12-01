#ifndef SIMPLE_EXAMPLE_MSP_READER_H
#define SIMPLE_EXAMPLE_MSP_READER_H

#include "spectrum.h"

using namespace std;

enum msp_read_mode {
    DETAILED = 0,
    PEAKS,
    BINNED_PEAKS,
    SPARSE
};

class msp_reader {

public:
    msp_reader();

    static vector<spectrum*> read_file(string path, msp_read_mode read_mode=DETAILED);


};


#endif //SIMPLE_EXAMPLE_MSP_READER_H
