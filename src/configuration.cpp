#include "configuration.h"
#include <iostream>
#include <fstream>

unsigned int configuration::assign_to_index(float mz) {
    for (int i = 0; i < (num_indices - 1); ++i) {
        if (mz < sub_idx_limits[i]) {
            return i;
        }
    }
    return num_indices - 1;
}

bool configuration::save_configuration(const std::string& config_file_path) {

    std::ofstream f(config_file_path, std::ios::out);
    std::string delimiter = ";";

    f << "Num indices: " << num_indices << "\n";
    f << "Index limits: ";
    for (unsigned int lim : sub_idx_limits) {
        f << lim << delimiter;
    }
    f << "\n";

    f.close();
    return true;
}

bool configuration::load_configuration(const std::string& config_file_path) {

    return true;
}
