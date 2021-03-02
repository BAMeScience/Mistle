#ifndef SIMPLE_EXAMPLE_SEARCH_MANAGER_H
#define SIMPLE_EXAMPLE_SEARCH_MANAGER_H

#include <string>
#include "library.h"

class search_manager {

    std::string search_file_path;
    std::string index_directory_path;

    library search_library;

    std::vector<std::vector<unsigned int>> mapped_search_ids; //TODO name right (bucket = subindex)

public:

    search_manager(std::string search_file_path, std::string index_directory_path);

    bool configure();
    bool prepare_search_library();


};


#endif //SIMPLE_EXAMPLE_SEARCH_MANAGER_H
