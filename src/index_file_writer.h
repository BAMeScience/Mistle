#ifndef SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H
#define SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H


#include <memory>
#include "spectrum.h"

class index_file_writer {


public:
    static bool stream_peaks_to_file(std::fstream &f, unsigned int parent_id, const std::shared_ptr<spectrum>& spec);
};


#endif //SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H
