#ifndef SIMPLE_EXAMPLE_SEARCH_MANAGER_H
#define SIMPLE_EXAMPLE_SEARCH_MANAGER_H

#include <string>
#include "library.h"
#include "configuration.h"

class search_manager {

    std::string search_file_path;
    std::string index_directory_path;
    std::shared_ptr<configuration> config;

    library search_library;

    std::vector<std::vector<unsigned int>> mapped_search_ids; //TODO name right (bucket = subindex)
    float mz_tolerance = 3.0;

public:

    search_manager(std::string search_file_path, std::string index_directory_path);

    bool prepare_search_library();
    bool schedule_searches();
    bool merge_matches(); //todo probably going over ids back to front and popping matches in the back


};


#endif //SIMPLE_EXAMPLE_SEARCH_MANAGER_H
