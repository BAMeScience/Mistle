#include "library.h"

#include <utility>
#include <iostream>
#include <filesystem>
#include "msp_reader.h"
#include "mgf_reader.h"
#include "settings.h"

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
        cout << "Loading library from directory:" << endl;
        load_library_from_directory(path);
    }
    cout << "Loading library from single file:" << endl;
    load_spectra_from_file(path);
}

bool library::construct(string &path) {
    if (path[path.length() - 1] == '/' || path[path.length() - 1] == '\\') {
        cout << "Loading library from directory:" << endl;
        load_library_from_directory(path);
    }
    cout << "Loading library from single file:" << endl;
    load_spectra_from_file(path);

    return true;
}

bool library::load_library_from_directory(string &path) {

    if (settings::load_batches) {
        cerr << "Warning: Batch search for multiple files not fully implemented" << endl;
        exit(1);
    }
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
        if (settings::load_batches) {
            file_stream.open(path, ios::in);
            last_batch = mgf_reader::read_file_batch(file_stream, spectrum_list, settings::batch_size);
        }
        else if (!mgf_reader::read_file(path, spectrum_list)) {
            cout << "Error reading file: " << path << endl;
            return false;
        }
    }
    else {
        cout << "Unknown file extension" << endl;
        return false;
    }
    std::cout << "\t" << spectrum_list.size() << " scans loaded" << std::endl;

    return true;
}




bool library::build_library_index() {
    cout << "Building precursor index" << endl;
    precursor_idx = new precursor_index();
    //TODO refactoring get up to date
    cout << "Building fragment ion index" << endl;
    fragment_ion_idx = new class fragment_ion_index(precursor_idx);
    is_indexed = true;
    return false;
}

library::library(vector<std::shared_ptr<spectrum>> &spectra) {
    spectrum_list = spectra;
}

bool library::load_next_batch() {
    last_batch = mgf_reader::read_file_batch(file_stream, spectrum_list, settings::batch_size);
    return last_batch;
}

