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
    exit(1);
}

indexing_manager::indexing_manager(string path) : path(path) {
    cout << "! Using IndexingManager without config parameter is deprecated" << endl;

    /*
     * Init
     * TODO delete (eventually)
     */

    precursorIndex = make_unique<precursor_index>();

    config->sub_idx_range = (STANDARD_PARENT_UPPER_MZ - STANDARD_PARENT_LOWER_MZ) / config->num_indices;
    for (int i = 1; i < config->num_indices; ++i) { //Starting from 1
        cout << "LIMIT: " << STANDARD_PARENT_LOWER_MZ + config->sub_idx_range * i << endl;
        config->sub_idx_limits.push_back(STANDARD_PARENT_LOWER_MZ + config->sub_idx_range * i);
    }


    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".msp") {
            lib_files.push_back(entry);
        }
    }
}


indexing_manager::indexing_manager(std::string path, std::shared_ptr<configuration> config) : path(path), config(config) {
    cout << "Calling Manager" << endl;

    /*
     * Init
     */

    precursorIndex = make_unique<precursor_index>();

    config->sub_idx_range = (STANDARD_PARENT_UPPER_MZ - STANDARD_PARENT_LOWER_MZ) / config->num_indices;
    for (int i = 1; i < config->num_indices; ++i) { //Starting from 1
        cout << "LIMIT: " << STANDARD_PARENT_LOWER_MZ + config->sub_idx_range * i << endl;
        config->sub_idx_limits.push_back(STANDARD_PARENT_LOWER_MZ + config->sub_idx_range * i);
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

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < lib_files.size(); ++i) {
        cout << "Parsing library file no. " << i << " (" << lib_files[i].path().filename() << ")" << endl;

        parse_file_buffered(i); //TODO choose best

    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;

    /*
     * Storing and rebuilding indices
     */

    //Precursor index
    cout << "Sorting precursors index" << endl;
    precursorIndex->sort_index();
    cout << "Saving ..." << endl;
    precursorIndex->save_index_to_file(config->idx_path + "precursor_idx.csv");
    config->save_configuration_to_file(config->idx_path + "config.txt");

    //Closing output streams and reopening them as input streams
    for (int i = 0; i < output_streams.size(); ++i) {
        output_streams[i].close();
    }

    cout << "Sorting fragment ion indices" << endl;
    for (int i = 0; i < config->sub_idx_file_names.size(); ++i) {
        string file_name = config->idx_path + config->sub_idx_file_names[i];

        //cout << "Loading fragment index " << config->sub_idx_file_names[i] << endl;
        fragment_ion_index frag_index(file_name);
        //cout << "Sorting ..." << endl;
        frag_index.sort_index(precursorIndex);
        //cout << "Saving..." << endl << endl;
        frag_index.save_index_to_binary_file(file_name);
    }
    cout << "Done" << endl;

    return true;
}

bool indexing_manager::set_up_output_streams() {

    for (int i = 0; i < config->num_indices; ++i) {
        string file_name = config->idx_path + "frag_idx_" + to_string(i) + ".bin";
        cout << file_name << endl;
        config->sub_idx_file_names.push_back("frag_idx_" + to_string(i) + ".bin");
        output_streams.emplace_back(fstream(file_name, std::ios::binary | std::ofstream::out));
    }


    return true;
}

bool indexing_manager::parse_file(unsigned int file_num) {
    string file_path = lib_files[file_num].path().string();

    ifstream f(file_path, ios::in);
    f.precision(FLOAT_OUTPUT_PRECISION);

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
        unsigned int idx_num = config->assign_to_index(bookmark.mz);
        index_file_writer::stream_peaks_to_file(output_streams[idx_num], bookmark.id, tmp_spectrum);


    }
    return true;
}

bool indexing_manager::parse_file_buffered(unsigned int file_num) {
    string file_path = lib_files[file_num].path().string();

    ifstream f(file_path, ios::in);
    f.precision(FLOAT_OUTPUT_PRECISION);

    unsigned int buffer_size = 40960;//1048576; //Byte //TODO fix if too small
    unsigned int carryover_pos = 0;
    string buffer;
    buffer.resize(buffer_size);

    /*
     * Main loop reading library and creating preliminary indices on the fly
     */
    while (!f.eof()) {

        //Read large char buffer
        f.read(&buffer[carryover_pos], buffer_size - carryover_pos);
        unsigned int last_pos = buffer.rfind("Name:");
        unsigned int current_pos = buffer.find("Name:");

        // Parse spectra within the buffer
        while (current_pos != last_pos) {
            unsigned int next_pos = buffer.find("Name:", current_pos + 1);
            shared_ptr<spectrum> tmp_spectrum = msp_reader::read_spectrum_from_buffer(buffer.substr(current_pos, next_pos - current_pos));

            //Save bookmark in precursor index
            precursor &bookmark = precursorIndex->record_new_precursor(tmp_spectrum);

            //Stream (binned) peaks into corresponding sub-index file
            unsigned int idx_num = config->assign_to_index(bookmark.mz);
            index_file_writer::stream_peaks_to_binary_file(output_streams[idx_num], bookmark.id, tmp_spectrum); //TODO made binary
            current_pos = next_pos;
        }

        carryover_pos = buffer_size - current_pos;
        buffer.replace(0, carryover_pos, buffer.substr(current_pos, std::string::npos)); //keep ms2 carryover in the buffer
        current_pos = 0; //TODO why is this 0 in the first place
        //buffer = buffer.substr(current_pos, std::string::npos));
        //buffer.resize(buffer_size);
    }
    return true;
}
