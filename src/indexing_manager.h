#ifndef SIMPLE_EXAMPLE_INDEXING_MANAGER_H
#define SIMPLE_EXAMPLE_INDEXING_MANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include "configuration.h"
#include "precursor_index.h"


class indexing_manager {

    std::string path;
    std::vector<std::filesystem::directory_entry> lib_files;


    //Precursor Index
    std::unique_ptr<precursor_index> precursorIndex;

    /*
     * (Sub-) Indices
     */

    std::shared_ptr<configuration> config = std::make_shared<configuration>();
    std::vector<std::fstream> output_streams;


public:
    indexing_manager();
    explicit indexing_manager(std::string path);
    indexing_manager(std::string path, std::shared_ptr<configuration> config);


    bool build_indices();
    bool set_up_output_streams();
    bool parse_file(unsigned int file_num);
    bool parse_file_buffered(unsigned int file_num);




};


#endif //SIMPLE_EXAMPLE_INDEXING_MANAGER_H
