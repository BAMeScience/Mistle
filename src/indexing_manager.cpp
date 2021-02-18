#include "indexing_manager.h"
#include <iostream>
#include <utility>
#include <filesystem>
#include <fstream>
#include "msp_reader.h"

using namespace std;

indexing_manager::indexing_manager() {

}

indexing_manager::indexing_manager(string path) : path(path) {
    cout << "Calling Manager" << endl;

    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".msp") {
            lib_files.push_back(entry);
        }
    }


}

bool indexing_manager::build_indices() {
    cout << "Manager too busy to answer phone. Calling back later" << endl;

    /*
     * Prepare in/output
     */
    set_up_output_streams();


    for (int i = 0; i < lib_files.size(); ++i) {
        cout << "Parsing library file no. " << i << " (" << lib_files[i].path().filename() << ")" << endl;

        parse_file(i);

    }


    return true;
}

unsigned int indexing_manager::assign_to_index(float mz) {
    for (int i = 0; i < (num_indices - 1); ++i) {
        if (mz < idx_limits[i]) {
            return i;
        }
    }
    return num_indices;
}


bool indexing_manager::set_up_output_streams() {

    for (int i = 0; i < num_indices; ++i) {
        string file_name = idx_path + "frag_idx_" + to_string(i) + ".csv";
        cout << file_name << endl;
        output_streams.emplace_back(ofstream(file_name, std::ofstream::out));
    }


    return true;
}

bool indexing_manager::parse_file(unsigned int file_num) {
    string file_path = lib_files[file_num].path().string();

    ifstream f(file_path);
    string buffer;
    spectrum * tmp_spectrum;

    while (true) {

        msp_reader::read_next_entry_into_buffer(f, buffer);
        tmp_spectrum = msp_reader::read_spectrum_from_buffer(buffer);
        //TODO PROCESS spectrum
        //Stream peaks into corresponding bin and sub-index

        //Save fingerprint into precursor index
        delete tmp_spectrum;
        break;
    }

    return false;
}
