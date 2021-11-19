#ifndef SIMPLE_EXAMPLE_CONFIGURATION_H
#define SIMPLE_EXAMPLE_CONFIGURATION_H


#include <string>
#include <vector>

/*
 * Config class
 *
 * Handles all (meta-)information about index set-up and configuration
 * Can be set by arguments or loaded from file.
 */

class configuration {
public:

    std::string idx_path = "./test/";
    std::string precursor_index_path;
    unsigned int num_indices = 24;

    unsigned int sub_idx_range;
    std::vector<unsigned int> sub_idx_limits;
    std::vector<std::string> sub_idx_file_names;
    unsigned int minimum_peptide_length;
    std::string build_command;

    //TODO parse more info and move to file_writer/reader
    bool save_configuration_to_file(const std::string& config_file_path);
    bool load_configuration_from_file(const std::string& config_file_path);


    unsigned int assign_to_index(float mz);


    /*
     * Build only
     */

    int num_build_threads = 1;

};


#endif //SIMPLE_EXAMPLE_CONFIGURATION_H
