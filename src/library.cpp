#include "library.h"

#include <utility>
#include <iostream>
#include <filesystem>
#include "msp_reader.h"
#include "mgf_reader.h"

using namespace std;

library::library() {

}

library::~library() {
    /*for (int i = 0; i < spectrum_list.size(); ++i) {
        delete spectrum_list[i];
    }*/
    spectrum_list.clear();
}


library::library(string &path) {
    if (path[path.length() - 1] == '/' || path[path.length() - 1] == '\\') {
        cout << "\nLoading library from directory:" << endl;
        load_library_from_directory(path);
    }
    cout << "\nLoading library from single file:" << endl;
    load_spectra_from_file(path);
}

bool library::load_library_from_directory(string &path) {

    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        load_spectra_from_file(entry.path().string());
        //cout << spectrum_list.size() << endl;
    }
    return true;
}

bool library::load_spectra_from_file(string path) {

    string extension = path.substr(path.rfind('.') + 1, string::npos);
    cout << "Loading spectra from file: " << path << endl;
    if (extension == "msp") {
        if (!msp_reader::read_file(path, spectrum_list)) {
            cout << "Error reading file: " << path << endl;
            return false;
        }
    }
    else if (extension == "mgf") {
        if (!mgf_reader::read_file(path, spectrum_list)) {
            cout << "Error reading file: " << path << endl;
            return false;
        }
    }
    else {
        cout << "Unknown file extension" << endl;
        return false;
    }

    return true;
}




bool library::build_library_index() {
    cout << "Building precursor index" << endl;
    precursor_idx = new class precursor_index(spectrum_list);
    cout << "Building fragment ion index" << endl;
    fragment_ion_idx = new class fragment_ion_index(precursor_idx);
    is_indexed = true;
    return false;
}

library::library(vector<spectrum *> &spectra) {
    spectrum_list = spectra;
}

