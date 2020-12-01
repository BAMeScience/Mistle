#include "msp_reader.h"
#include <iostream>

msp_reader::msp_reader() {

}

vector<spectrum *> msp_reader::read_file(string path, msp_read_mode read_mode) {
    path = "test";

    if (read_mode == msp_read_mode::DETAILED) {
        cout << "Reading in detailed mode" << endl;
    }
    return vector<spectrum *>();
}
