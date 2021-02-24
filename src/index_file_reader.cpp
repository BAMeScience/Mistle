#include "index_file_reader.h"
#include <fstream>
#include <iostream>


bool index_file_reader::read_file_into_precursor_index(const std::string &file_path,
                                                       const std::shared_ptr<precursor_index>& precursor_idx) {


    std::ifstream f(file_path, std::ios::in);
    std::string line;


    if (!getline(f, line)) {
        return false;
    }

    if (line.rfind("Num: ", 0) == 0) {
        std::cerr << "Incorrect file format" << std::endl;
        return false;
    }

    //Read header
    unsigned int size = std::stoi(line.substr(5, std::string::npos)); //TODO check 4 or 5
    precursor_idx->set_size(size);

    //
    //precursor_idx->add_precursor_record(p);


    return false;
}
