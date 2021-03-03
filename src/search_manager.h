#ifndef SIMPLE_EXAMPLE_SEARCH_MANAGER_H
#define SIMPLE_EXAMPLE_SEARCH_MANAGER_H

#include <string>
#include "library.h"
#include "configuration.h"


class search_manager {

    std::string search_file_path;
    std::string index_directory_path;
    std::shared_ptr<configuration> config;

    library search_library; //TODO probably replace by simple list of spectra

    //Mapped ms2 ids to sub-index where they might occur
    std::vector<std::vector<unsigned int>> mapped_search_ids; //TODO name right (bucket = subindex)
    float mz_tolerance = 3.0;

    /*
     * Indices
     */
    std::shared_ptr<precursor_index> precursor_idx;
    std::shared_ptr<fragment_ion_index> frag_idx;

public:

    search_manager(std::string search_file_path, std::string index_directory_path);

    bool prepare_search_library();
    bool prepare_precursor_index();
    bool perform_searches();
    bool merge_matches(); //todo probably going over ids back to front and popping matches in the back

private:
    bool search_spectrum(std::shared_ptr<spectrum> &spec);
};


#endif //SIMPLE_EXAMPLE_SEARCH_MANAGER_H
