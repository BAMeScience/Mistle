#include "library.h"

#include <utility>
#include <iostream>
#include "msp_reader.h"
#include "mgf_reader.h"

library::library() {

}

library::library(string &path) {
    load_library_from_file(path);
}

bool library::load_library_from_file(string &path) {

    //todo multiple file extensions, if endswith .msp do ...
    string extension = path.substr(path.rfind('.') + 1, string::npos);
    cout << "Loading library from file" << endl;
    if (extension == "msp") {
        spectrum_list = msp_reader::read_file(path);
    }
    else if (extension == "mgf") {
        spectrum_list = mgf_reader::read_file(path);
    }
    else {
        cerr << "Library creation not possible, unknown file extension" << endl;
    }

    return true;
}

library::~library() {
    /*for (int i = 0; i < spectrum_list.size(); ++i) {
        delete spectrum_list[i];
    }*/
    spectrum_list.clear();
}
