#include "library.h"

#include <utility>
#include <iostream>
#include <filesystem>
#include "msp_reader.h"
#include "mgf_reader.h"

library::library() {

}

library::~library() {
    /*for (int i = 0; i < spectrum_list.size(); ++i) {
        delete spectrum_list[i];
    }*/
    spectrum_list.clear();
}


library::library(string &path) {
    load_spectra_from_file(path);
}

bool library::load_spectra_from_file(string &path) {

    string extension = path.substr(path.rfind('.') + 1, string::npos);
    cout << "Loading library from file" << endl;
    if (extension == "msp") {
        if (!msp_reader::read_file(path, spectrum_list)) {
            cerr << "Error reading file: " << path << endl;
        }
    }
    else if (extension == "mgf") {
        if (!mgf_reader::read_file(path, spectrum_list)) {
            cerr << "Error reading file: " << path << endl;
        }
    }
    else {
        cerr << "unknown file extension" << endl;
    }

    return true;
}


bool library::load_library_from_directory(string &path) {

    cout << "Loading library from directory:" << endl;
    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        cout << entry << endl;
        load_spectra_from_file((string &) entry);
    }
    return true;
}

bool library::build_library_index() {
    cout << "Building precursor index" << endl;
    precursor_index = new class precursor_index(spectrum_list);
    cout << "Building fragment ion index" << endl;
    fragment_ion_index = new class fragment_ion_index(precursor_index);
    is_indexed = true;
    return false;
}

