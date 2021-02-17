#ifndef SIMPLE_EXAMPLE_INDEXING_MANAGER_H
#define SIMPLE_EXAMPLE_INDEXING_MANAGER_H

#include <string>
#include <vector>
#include <filesystem>

class indexing_manager {

    std::string path;
    std::vector<std::filesystem::directory_entry> lib_files;


public:
    indexing_manager();
    explicit indexing_manager(std::string path);


    bool build_indices();


};


#endif //SIMPLE_EXAMPLE_INDEXING_MANAGER_H
