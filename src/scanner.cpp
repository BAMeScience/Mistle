#include <filesystem>
#include <iostream>
#include "scanner.h"
#include "msp_reader.h"
#include "mgf_reader.h"

scanner::scanner() {

}

bool scanner::scan_directory(string path) {
    cout << "Scanning directory: " << path << endl;
    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        scan_file(entry.path().string());
    }
    return true;
}

bool scanner::scan_file(string path) {
    cout << "Scanning: " << path << endl;

    string extension = path.substr(path.rfind('.') + 1, string::npos);
    vector<spectrum *> &soondeleted();
    if (extension == "msp") {
        /*if (!msp_reader::read_file(path, soondeleted())) {
            cout << "Error reading file: " << path << endl;
            return false;
        }*/
    }
    else if (extension == "mgf") {
        /*if (!mgf_reader::read_file(path, soondeleted())) {
            cout << "Error reading file: " << path << endl;
            return false;
        }*/
    }
    else {
        cout << "Unknown file extension" << endl;
        return false;
    }

    return true;
    
    
    return false;
}
