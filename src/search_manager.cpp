#include <iostream>
#include "search_manager.h"

search_manager::search_manager(std::string search_file_path, std::string index_directory_path) : search_file_path(search_file_path), index_directory_path(index_directory_path) {
    std::cout << "Calling MS2 recruiter" << std::endl;

    std::cout << "Configuring ... " << std::endl;
    config = std::make_shared<configuration>();
    config->load_configuration(index_directory_path + "config.txt");


}


bool search_manager::prepare_search_library() {

    // Load search library
    search_library = library(search_file_path);

    // Divide ms2 into search container (corresponding to sub-index)

    for (int i = 0; i < search_library.spectrum_list.size(); ++i) {

        float mz = search_library.spectrum_list[i]->precursor_mass;
        float min_mz = mz - mz_tolerance;
        float max_mz = mz + mz_tolerance;

        // Map id (i) to every sub-index where min/max bounds fall into
        for (int idx_num = 0; idx_num < (config->num_indices - 1); ++idx_num) {
            if (min_mz < config->sub_idx_limits[idx_num]) {
                if (idx_num == 0 || max_mz > config->sub_idx_limits[idx_num - 1]) {
                    mapped_search_ids[idx_num].push_back(i);
                } else {
                    break; // because max_mz is below lower border of that and every subsequent sub-indices
                }
            }
        }
        if (max_mz > config->sub_idx_limits.back()) { // check if it falls into last sub index (without upper border)
            mapped_search_ids.back().push_back(i);
        }



    }

    return true;
}
