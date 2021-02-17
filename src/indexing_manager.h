#ifndef SIMPLE_EXAMPLE_INDEXING_MANAGER_H
#define SIMPLE_EXAMPLE_INDEXING_MANAGER_H

#include <string>
#include <vector>
#include <filesystem>

class indexing_manager {

    std::string path;
    std::vector<std::filesystem::directory_entry> lib_files;


    /*
     * Indices
     */
    unsigned int num_indices = 1;
    std::vector<unsigned int> idx_limits;


public:
    indexing_manager();
    explicit indexing_manager(std::string path);


    bool build_indices();
    bool parse_file(std::string file_path);


    unsigned int assign_to_index(float mz);


};


#endif //SIMPLE_EXAMPLE_INDEXING_MANAGER_H
