#include "configuration.h"

unsigned int configuration::assign_to_index(float mz) {
    for (int i = 0; i < (num_indices - 1); ++i) {
        if (mz < sub_idx_limits[i]) {
            return i;
        }
    }
    return num_indices - 1;
}

bool configuration::save_configuration(std::string config_file_path) {
    return false;
}

bool configuration::load_configuration(std::string config_file_path) {
    return false;
}
