#ifndef SIMPLE_EXAMPLE_INDEXING_MANAGER_H
#define SIMPLE_EXAMPLE_INDEXING_MANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>

class indexing_manager {

    std::string path;
    std::vector<std::filesystem::directory_entry> lib_files;


    /*
     * (Sub-) Indices
     */
    std::string idx_path = "./test/";
    unsigned int num_indices = 1;
    std::vector<unsigned int> idx_limits;
    std::vector<std::ofstream> output_streams;

public:
    indexing_manager();
    explicit indexing_manager(std::string path);


    bool build_indices();
    bool set_up_output_streams();
    bool parse_file(unsigned int file_num);


    unsigned int assign_to_index(float mz);


};


#endif //SIMPLE_EXAMPLE_INDEXING_MANAGER_H
