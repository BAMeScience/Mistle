#ifndef SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H
#define SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H


#include <memory>
#include "spectrum.h"
#include "precursor_index.h"
#include "match.h"

class index_file_writer {
public:

    //TODO static std::string delimiter;

    static bool stream_peaks_to_file(std::fstream &f, unsigned int parent_id, const std::shared_ptr<spectrum>& spec);
    static bool stream_peaks_to_binary_file(std::fstream &f, unsigned int parent_id, const std::shared_ptr<spectrum>& spec);
    static bool save_precursor_index(const std::string& file_path, std::vector<precursor> &precursors);
    static bool save_precursor_index_to_binary_file(const std::string& file_path, std::vector<precursor> &precursors);
    static bool save_matches_to_file(const std::string& file_path, std::vector<match> &matches);
};


#endif //SIMPLE_EXAMPLE_INDEX_FILE_WRITER_H
