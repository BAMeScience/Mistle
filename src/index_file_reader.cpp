#include "index_file_reader.h"
#include <fstream>
#include <iostream>


bool index_file_reader::read_file_into_precursor_index(const std::string &file_path,
                                                       const std::shared_ptr<precursor_index>& precursor_idx) {


    std::ifstream f(file_path, std::ios::in);
    std::string delimiter = ";";
    std::string line;

    if (!getline(f, line)) {
        return false;
    }

    if (line.rfind("Num: ", 0) != 0) {
        std::cerr << "Incorrect file format" << std::endl;
        return false;
    }

    //Read header
    unsigned int size = std::stoi(line.substr(5, std::string::npos)); //TODO check 4 or 5
    precursor_idx->set_size(size);

    // Parse precursors line by line
    //precursor_idx->add_precursor_record(p);
    while (getline(f, line)) {

        size_t delim_pos = line.find(delimiter);
        unsigned int id = std::stoi(line.substr(0, delim_pos));
        
        size_t length = line.find(delimiter, delim_pos + 1) - delim_pos;
        unsigned int rank = std::stoi(line.substr(delim_pos + 1, length - 1));

        delim_pos = delim_pos + length;
        length = line.find(delimiter, delim_pos + 1) - delim_pos;
        float mz = std::stof(line.substr(delim_pos + 1, length - 1));

        delim_pos = delim_pos + length;
        length = line.find(delimiter, delim_pos + 1) - delim_pos;
        int charge = std::stoi(line.substr(delim_pos + 1, length - 1));
        std::string peptide = line.substr(delim_pos + length + 1, std::string::npos);

        precursor_idx->add_precursor_record(precursor(id, rank, mz, charge, peptide));
    }
    if (precursor_idx->get_size() != size) {
        std::cerr << "Wrong number of precursors" << std::endl;
    }

    return true;
}
