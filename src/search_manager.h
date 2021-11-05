#ifndef SIMPLE_EXAMPLE_SEARCH_MANAGER_H
#define SIMPLE_EXAMPLE_SEARCH_MANAGER_H

#include <string>
#include <chrono>
#include "library.h"
#include "configuration.h"
#include "match.h"
#include "thread_pool.h"


class search_manager {

    std::string search_file_path;
    std::string index_directory_path;
    std::shared_ptr<configuration> config;

    bool last_batch = true;
    library search_library; //TODO probably replace by simple list of spectra

    //Mapped ms2 ids to sub-index where they might occur
    std::vector<std::vector<unsigned int>> mapped_search_ids; //TODO name right (bucket = subindex)

    /*
     * Indices
     */
    std::shared_ptr<precursor_index> precursor_idx;
    std::shared_ptr<fragment_ion_index> frag_idx;

    /*
     * Threading
     */

    std::shared_ptr<thread_pool> pool;

    /*
     * Results
     */
    long starting_time;
    long total_time_elapsed;
    std::vector<match> matches;

    /*
     * Timer
     */

    std::chrono::duration<double> inner_search_duration;

public:

    search_manager(std::string search_file_path, std::string index_directory_path);

    bool prepare_search_library();
    bool prepare_precursor_index();
    bool perform_searches();
    bool perform_searches_parallel();
    bool merge_matches(); //todo probably going over ids back to front and popping matches in the back
    bool save_search_results_to_file(const std::string &file_path);

    long get_time_spent_in_inner_search();

private:
    /*
     * Search: 3-fold implementation, depending on cpu instruction level (STANDARD, AVX2, AVX512)
     */
    bool search_spectrum(unsigned int search_id);
    float rescore_spectrum(unsigned int search_id, unsigned int target_id);
    bool rescore_match_old(match &psm);
    bool rescore_match(match &psm);
    bool prepare_next_batch();


    float sigma;
    float max_normal;
    static float normal_pdf(float x, float mean, float standard_deviation);
    static int factorial(int n);
};


#endif //SIMPLE_EXAMPLE_SEARCH_MANAGER_H
