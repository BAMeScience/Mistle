#include <filesystem>
#include <iostream>
#include <fstream>
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

    if (extension == "msp") {
        if (!msp_reader::read_file_precursors(path, parents)) {
            cout << "Error reading file: " << path << endl;
            return false;
        }

        ++no_files_read;
        int size = std::filesystem::file_size(path) / 1024;
        cout << "File size " << size  << " KB" << endl;
        kb_lib_size += size;
    }
    else if (extension == "mgf") {
        /*if (!mgf_reader::read_file(path, soondeleted())) {
            cout << "Error reading file: " << path << endl;
            return false;
        }*/
        cout << ".mgf quick scan not implemented yet" << endl;
        return false;
    }
    else {
        cout << "Unknown file extension" << endl;
        return false;
    }

    return true;
    
    
}

bool scanner::analyze() {

    sort(parents.begin(), parents.end(), [](const precursor *a, const precursor *b) {
        return *a < *b;
    });

    return true;
}

bool scanner::save_precursor_distribution_to_file(string path, string delimiter) {

    fstream outfile;
    outfile.open(path, ios::out);

    if (!outfile.good())
        return false;

    // Add header
    outfile << "mass"+delimiter+"charge" << endl;

    // Go through matches and parse relevant information for each
    for (int i = 0; i < parents.size(); ++i) {
        precursor *p = parents[i];
        outfile << p->mass << delimiter << p->charge << endl;
    }

    outfile.close();
    return true;

    return false;
}

bool scanner::print_scan_results() {

    cout << "Readable files detected: " << no_files_read << endl;
    cout << "Total size: " << kb_lib_size / 1024 << " MB (" << float(kb_lib_size) / float(1024*1024) << " GB)" << endl;

    return false;
}
