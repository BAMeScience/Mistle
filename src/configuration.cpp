#include "configuration.h"
#include <vector>
#include <string>
#include <sstream>
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

bool configuration::save_configuration_to_file(const std::string& config_file_path) {

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

bool configuration::load_configuration_from_file(const std::string& config_file_path) {

    std::ifstream f(config_file_path, std::ios::in);
    std::string delimiter = ";";

    //First line
    std::string line;
    getline(f, line);

    if(line.rfind("Num indices: ", 0) == 0) {
        num_indices = std::stoi(line.substr(13, std::string::npos));
    } else {
        std::cerr << "Wrong config format" << std::endl;
        return false;
    }

    //Second line
    getline(f, line);
    if(line.rfind("Index limits: ", 0) == 0) {
        std::stringstream ss(line.substr(14, std::string::npos));
        std::string str;
        while(getline(ss, str, ';')) {
            sub_idx_limits.push_back(std::stoi(str));
        }
        if (sub_idx_limits.size() != num_indices - 1) {
            std::cerr << "Num sub idx not matching" << std::endl;
            return false;
        }
    } else {
        std::cerr << "Wrong config format" << std::endl;
        return false;
    }

    f.close();

    idx_path = config_file_path.substr(0,config_file_path.rfind('/') + 1);
    precursor_index_path = idx_path + "precursor_idx.csv";
    for (int i = 0; i < num_indices; ++i) {
        sub_idx_file_names.push_back(idx_path + "frag_idx_" + std::to_string(i) + ".bin");
    }



    return true;
}
