#ifndef SIMPLE_EXAMPLE_CONFIGURATION_H
#define SIMPLE_EXAMPLE_CONFIGURATION_H


#include <string>
#include <vector>

class configuration {
public:

    std::string idx_path = "./test/";
    unsigned int num_indices = 8;

    unsigned int sub_idx_range;
    std::vector<unsigned int> sub_idx_limits;
    std::vector<std::string> sub_idx_file_names;

    bool save_configuration(std::string config_file_path);
    bool load_configuration(std::string config_file_path);
    

    unsigned int assign_to_index(float mz);


};


#endif //SIMPLE_EXAMPLE_CONFIGURATION_H
