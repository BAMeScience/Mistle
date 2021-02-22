#ifndef SIMPLE_EXAMPLE_INDEXING_MANAGER_H
#define SIMPLE_EXAMPLE_INDEXING_MANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include "precursor_index.h"

class indexing_manager {

    std::string path;
    std::vector<std::filesystem::directory_entry> lib_files;


    //Precursor Index
    std::unique_ptr<precursor_index> precursorIndex;

    /*
     * (Sub-) Indices
     */
    std::string idx_path = "./test/";
    unsigned int num_indices = 8;

    unsigned int sub_idx_range;
    std::vector<unsigned int> sub_idx_limits;
    std::vector<std::fstream> output_streams;
public:
    indexing_manager();
    explicit indexing_manager(std::string path);


    bool build_indices();
    bool set_up_output_streams();
    bool parse_file(unsigned int file_num);


    unsigned int assign_to_index(float mz);


};


#endif //SIMPLE_EXAMPLE_INDEXING_MANAGER_H
