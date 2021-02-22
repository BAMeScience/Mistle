#include "indexing_manager.h"
#include <iostream>
#include <utility>
#include <filesystem>
#include <fstream>
#include "msp_reader.h"
#include "index_file_writer.h"
#include "DefineConstants.h"
#include "fragment_ion_index.h"

using namespace std;

indexing_manager::indexing_manager() {
    cout << "Empty Constructor not in use" << endl;
}

indexing_manager::indexing_manager(string path) : path(path) {
    cout << "Calling Manager" << endl;

    /*
     * Init
     */

    precursorIndex = make_unique<precursor_index>();

    sub_idx_range = (STANDARD_PARENT_UPPER_MZ - STANDARD_PARENT_LOWER_MZ) / num_indices;
    for (int i = 1; i < num_indices; ++i) { //Starting from 1
        cout << "LIMIT: " << STANDARD_PARENT_LOWER_MZ + sub_idx_range * i << endl;
        sub_idx_limits.push_back(STANDARD_PARENT_LOWER_MZ + sub_idx_range * i);
    }


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


    /*
     * Parsing files and creating preliminary indices
     */

    for (int i = 0; i < lib_files.size(); ++i) {
        cout << "Parsing library file no. " << i << " (" << lib_files[i].path().filename() << ")" << endl;

        parse_file(i);

    }

    /*
     * Storing and rebuilding indices
     */

    //Precursor index
    precursorIndex->sort_index();
    //TODO SAVE

    //Closing output streams and reopening them as input streams
    for (int i = 0; i < output_streams.size(); ++i) {
        output_streams[i].close();
    }

    for (int i = 0; i < sub_idx_file_names.size(); ++i) {
        string file_name = idx_path + sub_idx_file_names[i];

        cout << "Loading unsorted index " << sub_idx_file_names[i] << endl;
        fragment_ion_index frag_index(file_name);
        cout << "Sorting ..." << endl;
        frag_index.sort_index(precursorIndex);
        cout << "Saving..." << endl << endl;
        frag_index.save_index_to_file(file_name);
    }


    return true;
}

unsigned int indexing_manager::assign_to_index(float mz) {
    for (int i = 0; i < (num_indices - 1); ++i) {
        if (mz < sub_idx_limits[i]) {
            return i;
        }
    }
    return num_indices - 1;
}


bool indexing_manager::set_up_output_streams() {

    for (int i = 0; i < num_indices; ++i) {
        string file_name = idx_path + "frag_idx_" + to_string(i) + ".csv";
        cout << file_name << endl;
        sub_idx_file_names.push_back("frag_idx_" + to_string(i) + ".csv");
        output_streams.emplace_back(fstream(file_name, std::ofstream::out));
    }


    return true;
}

bool indexing_manager::parse_file(unsigned int file_num) {
    string file_path = lib_files[file_num].path().string();

    ifstream f(file_path, ios::in);
    string buffer;


    /*
     * Main loop reading library and creating preliminary indices on the fly
     */

    while (!f.eof()) {

        //Read spectrum from file and pre-processing
        msp_reader::read_next_entry_into_buffer(f, buffer);
        shared_ptr<spectrum> tmp_spectrum = msp_reader::read_spectrum_from_buffer(buffer);

        //Save bookmark in precursor index
        precursor &bookmark = precursorIndex->record_new_precursor(tmp_spectrum);

        //Stream (binned) peaks into corresponding sub-index file
        unsigned int idx_num = assign_to_index(bookmark.mass);
        index_file_writer::stream_peaks_to_file(output_streams[idx_num], bookmark.id, tmp_spectrum);


    }
    return true;
}
