#ifndef SIMPLE_EXAMPLE_CONFIGURATION_H
#define SIMPLE_EXAMPLE_CONFIGURATION_H


#include <string>
#include <vector>

class configuration {
public:

    std::string idx_path = "./test/";
    std::string precursor_index_path;
    unsigned int num_indices = 24;

    unsigned int sub_idx_range;
    std::vector<unsigned int> sub_idx_limits;
    std::vector<std::string> sub_idx_file_names;

    //TODO parse more info and move to file_writer/reader
    bool save_configuration_to_file(const std::string& config_file_path);
    bool load_configuration_from_file(const std::string& config_file_path);


    unsigned int assign_to_index(float mz);


};


#endif //SIMPLE_EXAMPLE_CONFIGURATION_H
