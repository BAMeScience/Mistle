#include <iostream>
#include "search_manager.h"

search_manager::search_manager(std::string search_file_path, std::string index_directory_path) : search_file_path(search_file_path), index_directory_path(index_directory_path) {
    std::cout << "Calling search manager" << std::endl;

}

bool search_manager::configure() {
    return true;
}

bool search_manager::prepare_search_library() {

    // Load search library
    search_library = library(search_file_path);

    // Divide ms2 into search container (corresponding to sub-index)
    // TODO

    return true;
}
