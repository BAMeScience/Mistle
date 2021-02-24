#ifndef SIMPLE_EXAMPLE_INDEX_FILE_READER_H
#define SIMPLE_EXAMPLE_INDEX_FILE_READER_H


#include <string>
#include "precursor_index.h"

class index_file_reader {
public:
    static bool read_file_into_precursor_index(const std::string &file_path, const std::shared_ptr<precursor_index>& precursor_idx);
};


#endif //SIMPLE_EXAMPLE_INDEX_FILE_READER_H
