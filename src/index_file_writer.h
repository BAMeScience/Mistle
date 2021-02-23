#ifndef SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H
#define SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H


#include <memory>
#include "spectrum.h"
#include "precursor_index.h"

class index_file_writer {


public:
    static bool stream_peaks_to_file(std::fstream &f, unsigned int parent_id, const std::shared_ptr<spectrum>& spec);
    static bool save_precursor_index(const std::string& file_path, std::vector<precursor> &precursors);
};


#endif //SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H
