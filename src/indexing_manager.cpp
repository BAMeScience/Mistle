#include "indexing_manager.h"
#include <iostream>
#include <utility>
#include <filesystem>
#include <fstream>
#include "msp_reader.h"
#include "mgf_reader.h"
#include "index_file_writer.h"
#include "DefineConstants.h"
#include "fragment_ion_index.h"

using namespace std;

indexing_manager::indexing_manager() {
    cout << "Empty Constructor not in use" << endl;
    exit(1);
}

indexing_manager::indexing_manager(string path) {
    cout << "! Using IndexingManager without config parameter is deprecated" << endl;
    exit(1);
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


indexing_manager::indexing_manager(std::vector<std::string> &input_paths, std::shared_ptr<configuration> config) : input_paths(input_paths), config(config) {

    /*
     * Init
     */

    precursorIndex = make_unique<precursor_index>();
    pool = std::make_shared<thread_pool>(config->num_build_threads - 1);

    config->sub_idx_range = (STANDARD_PARENT_UPPER_MZ - STANDARD_PARENT_LOWER_MZ) / config->num_indices;
    for (int i = 1; i < config->num_indices; ++i) { //Starting from 1
        //cout << "LIMIT: " << STANDARD_PARENT_LOWER_MZ + config->sub_idx_range * i << endl;
        config->sub_idx_limits.push_back(STANDARD_PARENT_LOWER_MZ + config->sub_idx_range * i);
    }


    for (std::string &path : input_paths) {

        if (!std::filesystem::exists(path)) {
            std::cerr << "Bad file" << std::endl;
            std::cerr << path << " is broken or does not exist!" << std::endl;
            exit(1);
        }

        if (std::filesystem::is_directory(path)) {
            for (const auto & entry : std::filesystem::directory_iterator(path)) {
                if (entry.path().extension() == ".msp") {
                    lib_files.push_back(entry.path());
                }
            }
        } else {
            std::filesystem::path file_path = path;
            if (file_path.extension() == ".msp") {
                file_format = MSP;
                lib_files.push_back(file_path);
            } else if (file_path.extension() == ".mgf") {
                file_format = MGF;
                lib_files.push_back(file_path);
            } else {
                std::cerr << "Unsupported file extension" << std::endl;
                std::cerr << file_path << " file is not supported" << std::endl;
                exit(1);
            }
        }

    }
}


bool indexing_manager::build_indices() {

    /*
     * Prepare in/output
     */
    set_up_output_streams();


    /*
     * Parsing files and creating preliminary indices
     */

    cout << "Parsing "<< lib_files.size() <<" library files ..." << endl;
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < lib_files.size(); ++i) {
        //cout << "Parsing library file no. " << i << " (" << lib_files[i].path().filename() << ")" << endl;
        parse_file(i); //TODO has multi-threading (experimental)
    }
    if (pool->get_size() > 0) {
        std::cout << "Waiting for threads to finish processing" << std::endl;
        pool->add_thread(); //Have "main" thread help out with the computation
        pool->wait_for_all_threads();
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
    //precursorIndex->save_index_to_file(config->idx_path + "precursor_idx.csv");
    precursorIndex->save_index_to_binary_file(config->idx_path + "precursor_idx.bin");

    config->save_configuration_to_file(config->idx_path + "config.txt");

    //Closing output streams and reopening them as input streams
    for (int i = 0; i < output_streams.size(); ++i) {
        output_streams[i].close();
    }

    cout << "Sorting fragment ion indices" << endl;
    for (int i = 0; i < config->sub_idx_file_names.size(); ++i) {
        string file_name = config->idx_path + config->sub_idx_file_names[i];

        fragment_ion_index frag_index;
        frag_index.load_preliminary_index_from_binary_file(file_name);
        frag_index.sort_index(precursorIndex);
        frag_index.save_index_to_binary_file(file_name);
    }
    cout << "Done" << endl;

    return true;
}

bool indexing_manager::set_up_output_streams() {

    for (int i = 0; i < config->num_indices; ++i) {
        string file_name = config->idx_path + "frag_idx_" + to_string(i) + ".bin";
        //cout << file_name << endl;
        config->sub_idx_file_names.push_back("frag_idx_" + to_string(i) + ".bin");
        output_streams.emplace_back(fstream(file_name, std::ios::binary | std::ofstream::out));
    }


    return true;
}

bool indexing_manager::parse_file(unsigned int file_num) {
    string file_path = lib_files[file_num].string();

    ifstream f(file_path, ios::in);
    //f.precision(FLOAT_OUTPUT_PRECISION);

    string buffer;


    /*
     * Main loop reading library and creating preliminary indices on the fly
     */

    while (!f.eof()) {

        //Read spectrum from file and pre-processing
        bool read_successfully = false;
        if (file_format == MSP) {
            read_successfully = msp_reader::read_next_entry_into_buffer(f, buffer);
        } else if (file_format == MGF) {
            read_successfully = mgf_reader::read_next_entry_into_buffer(f, buffer);
        }
        if (!read_successfully)
            continue;

        auto read_and_stream = [this, buffer]() {

            shared_ptr<spectrum> tmp_spectrum;
            if (file_format == MSP) {
                tmp_spectrum = msp_reader::read_spectrum_from_buffer(buffer);
            } else if (file_format == MGF) {
                tmp_spectrum = mgf_reader::read_spectrum_from_buffer(buffer);
            }

            if (tmp_spectrum->peptide.length() < config->minimum_peptide_length) {
                return;
            }


            //Lock for recording and streaming
            std::lock_guard<std::mutex> guard(pool->mtx);

            //Save bookmark in precursor index
            precursor &bookmark = precursorIndex->record_new_precursor(tmp_spectrum);

            //Stream (binned) peaks into corresponding sub-index file

            unsigned int idx_num = config->assign_to_index(bookmark.mz);
            index_file_writer::stream_peaks_to_binary_file(output_streams[idx_num], bookmark.id, tmp_spectrum);
            //std::cout << bookmark.mz << std::endl;

        };


        if (pool->get_size() > 0) {
            pool->enqueue(read_and_stream);
        } else {
            read_and_stream();
        }
    }

    return true;
}

bool indexing_manager::parse_file_buffered(unsigned int file_num) {
    string file_path = lib_files[file_num].string();

    ifstream f(file_path, ios::in);
    //f.precision(FLOAT_OUTPUT_PRECISION); //TODO

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
            index_file_writer::stream_peaks_to_binary_file(output_streams[idx_num], bookmark.id, tmp_spectrum);
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
