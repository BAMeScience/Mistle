#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>
//#include <avx512vlintrin.h>
#include <x86intrin.h>

//#include <x86intrin.h>
#include <immintrin.h>
#include <complex>
#include <cmath>

#include "search_manager.h"
#include "settings.h"
#include <set>
search_manager::search_manager(std::string search_file_path, std::string index_directory_path) : search_file_path(search_file_path), index_directory_path(index_directory_path) {
    std::cout << "Configuring ... " << std::endl;
    config = std::make_shared<configuration>();
    config->load_configuration_from_file(index_directory_path + "config.txt");
    pool = std::make_shared<thread_pool>(settings::num_threads);

    //Setting up scoring parameters
    sigma = settings::bin_size; // Alternative: / 2.f;
    max_normal = normal_pdf(0,0, sigma);

}


bool search_manager::prepare_search_library() {

    // Load search library
    search_library.construct(search_file_path);
    last_batch = search_library.last_batch;

    // Map ms2 ids to search container (corresponding to sub-index)
    mapped_search_ids = std::vector<std::vector<unsigned int>>(config->num_indices);

    if (config->num_indices == 1) {
        mapped_search_ids[0].resize(search_library.spectrum_list.size());
        std::iota(mapped_search_ids[0].begin(), mapped_search_ids[0].end(), 0);
        return true;
    }

    for (unsigned int i = 0; i < search_library.spectrum_list.size(); ++i) {

        float mz = search_library.spectrum_list[i]->precursor_mass;
        float min_mz, max_mz;
        if (settings::use_ppm_tolerance) {
            float mz_deviation = (mz * settings::ppm_factor);
            min_mz = mz - mz_deviation;
            max_mz = mz + mz_deviation;
        } else {
            min_mz = mz - settings::mz_tolerance;
            max_mz = mz + settings::mz_tolerance;
        }

        // Map id (i) to every sub-index where min/max bounds fall into
        for (int idx_num = 0; idx_num < (config->num_indices - 1); ++idx_num) {
            if (min_mz <= config->sub_idx_limits[idx_num]) {
                if (idx_num == 0 || max_mz >= config->sub_idx_limits[idx_num - 1]) {
                    mapped_search_ids[idx_num].push_back(i);
                    ++search_library.spectrum_list[i]->search_counter;
                } else {
                    break; // because max_mz is below lower border of that and every subsequent sub-indices todo >(equality)< not needed for both I think, just there to be sure
                }
            }
        }
        if (max_mz >= config->sub_idx_limits.back()) { // check if it falls into last sub index (without upper border)
            mapped_search_ids.back().push_back(i);
            ++search_library.spectrum_list[i]->search_counter;
        }

    }

    return true;
}


bool search_manager::prepare_next_batch() {

    std::cout << "Loading next batch" << std::endl;
    //Delete Spectra not being searched anymore
    /*for (int i = 0; i < search_library.spectrum_list.size(); ++i) {
        if (search_library.spectrum_list[i]->search_counter == 0) {
            search_library.spectrum_list[i] = nullptr;
        }
    }*/


    // Load
    int k = search_library.spectrum_list.size();
    last_batch = search_library.load_next_batch();

    if (config->num_indices == 1) { //TODO TEST THIS
        mapped_search_ids[0].resize(search_library.spectrum_list.size() - k);
        std::iota(mapped_search_ids[0].begin(), mapped_search_ids[0].end(), k);
        return true;
    }

    for (unsigned int i = k; i < search_library.spectrum_list.size(); ++i) {

        float mz = search_library.spectrum_list[i]->precursor_mass;
        float min_mz, max_mz;
        if (settings::use_ppm_tolerance) {
            float mz_deviation = (mz * settings::ppm_factor);
            min_mz = mz - mz_deviation;
            max_mz = mz + mz_deviation;
        } else {
            min_mz = mz - settings::mz_tolerance;
            max_mz = mz + settings::mz_tolerance;
        }

        // Map id (i) to every sub-index where min/max bounds fall into
        for (int idx_num = 0; idx_num < (config->num_indices - 1); ++idx_num) {
            if (min_mz <= config->sub_idx_limits[idx_num]) {
                if (idx_num == 0 || max_mz >= config->sub_idx_limits[idx_num - 1]) {
                    mapped_search_ids[idx_num].push_back(i);
                    ++search_library.spectrum_list[i]->search_counter;
                } else {
                    break; // because max_mz is below lower border of that and every subsequent sub-indices todo >(equality)< not needed for both I think, just there to be sure
                }
            }
        }
        if (max_mz >= config->sub_idx_limits.back()) { // check if it falls into last sub index (without upper border)
            mapped_search_ids.back().push_back(i);
            ++search_library.spectrum_list[i]->search_counter;
        }

    }

    return true;
}

bool search_manager::prepare_precursor_index() {
    precursor_idx = std::make_shared<precursor_index>();
    precursor_idx->load_index_from_binary_file(config->precursor_index_path);
    return true;
}

bool search_manager::perform_searches() {

    if (settings::parallel) {
        return perform_searches_parallel();
    }
    frag_idx = std::make_shared<fragment_ion_index>();
    frag_idx->precursor_idx = precursor_idx;
    for (int i = 0; i < config->num_indices; ++i) {
        //std::cout << "Loading index number " << i << std::endl;
        if (!last_batch && mapped_search_ids[i].size() < settings::batch_size / config->num_indices)
            continue;
        if (mapped_search_ids[i].empty()) {
            continue;
        }
        frag_idx->load_index_from_binary_file(config->sub_idx_file_names[i]);
        frag_idx->prepare_axv_access();
        //TODO set precursor index limits by subindex borders... has to be properly implemented
        //std::cout << "Searching ... " << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        std::vector<unsigned int> &search_ids = mapped_search_ids[i];
        for (unsigned int &s_id : search_ids) {
            search_spectrum(s_id);
        }
        mapped_search_ids[i].clear();
        auto t2 = std::chrono::high_resolution_clock::now();
        inner_search_duration += (t2 - t1);
    }

    /*
     * Load next batch
     */
    if (!last_batch) {
        prepare_next_batch();
        perform_searches();
    }

    return true;
}


bool search_manager::perform_searches_parallel() {
    frag_idx = std::make_shared<fragment_ion_index>();
    frag_idx->precursor_idx = precursor_idx;

    for (int i = 0; i < config->num_indices; ++i) {
        //std::cout << "Loading index number " << i << std::endl;
        if (!last_batch && mapped_search_ids[i].size() < settings::batch_size / config->num_indices)
            continue;
        if (mapped_search_ids[i].empty()) {
            continue;
        }
        frag_idx->load_index_from_binary_file(config->sub_idx_file_names[i]);
        frag_idx->prepare_axv_access();
        //TODO set precursor index limits by subindex borders... has to be properly implemented
        //std::cout << "Searching ... " << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        std::vector<unsigned int> &search_ids = mapped_search_ids[i];
        for (unsigned int &s_id : search_ids) {
            std::shared_ptr<spectrum> spec = search_library.spectrum_list[s_id];

            //auto f = [&](int s) {test(s);};
            //pool->enqueue([f, s_id] { return f(s_id); });

            pool->enqueue([this, s_id] () {search_spectrum(s_id);});
        }
        pool->wait_for_all_threads();

        mapped_search_ids[i].clear();
        auto t2 = std::chrono::high_resolution_clock::now();
        inner_search_duration += (t2 - t1);
        /*
         * Load next batch
         */
        if (!last_batch) {
            prepare_next_batch();
            perform_searches_parallel(); //TODO correct?
        }
    }
    return true;
}

#if USE_AVX_512
bool search_manager::search_spectrum(unsigned int search_id) {

    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];

    // Determine range of candidate spectra
    float mz_tolerance = settings::mz_tolerance;
    if (settings::use_ppm_tolerance) {
        mz_tolerance = settings::ppm_factor * spec->precursor_mass;
    }
    int lower_rank = precursor_idx->get_lower_bound(spec->charge,spec->precursor_mass - mz_tolerance);
    int upper_rank = precursor_idx->get_upper_bound(spec->charge,spec->precursor_mass + mz_tolerance);

    if (lower_rank < 0 || upper_rank < 0 || lower_rank > upper_rank) {
        return false;
    }


    // Init candidate scores
    std::vector<float> dot_scores(upper_rank - lower_rank + 1, 0.f);

    // Update scores by matching all peaks using the fragment ion index
    for (int j = 0; j < spec->binned_peaks.size(); ++j) {

        // Open ion mz bin for corresponding peak
        fragment_bin &ion_bin = frag_idx->fragment_bins[spec->binned_peaks[j]];
        fragment_binn &bin = frag_idx->frag_bins[spec->binned_peaks[j]];


        // Determine starting point of lowest (candidate) parent index inside bin
        int starting_point_inside_bin = std::lower_bound(ion_bin.begin(), ion_bin.end(), lower_rank, [&](fragment &f, int rank) {
            return precursor_idx->get_rank(f.parent_id) < rank;
        }) - ion_bin.begin();
        int end_point = std::upper_bound(ion_bin.begin() + starting_point_inside_bin, ion_bin.end(), upper_rank, [&](int rank, fragment &f) {
            return precursor_idx->get_rank(f.parent_id) > rank;
        }) - ion_bin.begin();

        //Obtain dot-scores iff in valid range (0 otherwise
        auto valid_score = [&](int ii) {
            if (ii < end_point)
                return dot_scores[precursor_idx->get_rank(ion_bin[ii].parent_id) - lower_rank];
            return 0.f;
        };
        auto score_pos_secure = [&](int ii) {
            if (ii < end_point)
                return precursor_idx->get_rank(ion_bin[ii].parent_id) - lower_rank;
            return (unsigned) 0;
        };
        auto score_pos = [&](int ii) {
            return precursor_idx->get_rank(ion_bin[ii].parent_id);
        };



        /*
        __m256 _scalar = _mm256_set1_ps(spec->binned_intensities[j]);
        __m256i _lower_rank = _mm256_set1_epi32(lower_rank);
         */
        //Update scores for all parents with fragments in the range
        /*for (int k = starting_point_inside_bin; k < end_point; k+=8) {

            //Fill vectors with 8 float values
            __m256 _mini_vector = _mm256_loadu_ps(&bin.intensities[k]);
            __m256i _score_pos = _mm256_setr_epi32(score_pos_secure(k), score_pos_secure(k+1), score_pos_secure(k+2), score_pos_secure(k+3),
                                                  score_pos_secure(k+4), score_pos_secure(k+5), score_pos_secure(k+6), score_pos_secure(k+7));
            __m256 _scores = _mm256_i32gather_ps(&dot_scores[0], _score_pos, 4); //TODO why 4????

            //__m256 _scoresBefore = {valid_score(k), valid_score(k+1), valid_score(k+2), valid_score(k+3),
            //                     valid_score(k+4), valid_score(k+5), valid_score(k+6), valid_score(k+7)};


            __m256 _result = _mm256_fmadd_ps(_scalar, _mini_vector, _scores);

            //_mm256_i32scatter_ps(&dot_scores[0], _score_pos, _result, 4);
            for (int l = 0; l < 8 && k + l < end_point; ++l) {
                dot_scores[precursor_idx->get_rank(ion_bin[k + l].parent_id) - lower_rank] = _result[l];
            }
        }*/

        __m512 _scalar = _mm512_set1_ps(spec->binned_intensities[j]);
        __m512i _lower_rank = _mm512_set1_epi32(lower_rank);
        __m512 _mini_vector, _scores, _result;
        __m512i _score_pos0, _score_pos1, _score_pos2;

        int k = starting_point_inside_bin;
        //First reach k mod 16 = 0
        for (; k % 16 != 0 && k < end_point; ++k) {
            fragment &f = ion_bin[k];
            dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
        }
        for (; k+16 < end_point; k+=16) {
            //Fill vectors with 8 float values
            _mini_vector = bin._intensities[k/16];
            //__m256 _mini_vector = _mm256_load_ps(&bin.intensities[k]);
            _score_pos0 = bin._parent_ids[k/16];
            _score_pos1 = _mm512_i32gather_epi32(_score_pos0, &precursor_idx->to_rank[0], 4);
                    //_mm512_i32gather_epi32(&precursor_idx->to_rank[0], _score_pos0, 4); //TODO why 4????
            _score_pos2 = _mm512_sub_epi32(_score_pos1, _lower_rank);
            _scores = _mm512_i32gather_ps(_score_pos2, &dot_scores[0], 4); //TODO why 4????

            //__m256i _score_pos1 = _mm256_setr_epi32(score_pos(k), score_pos(k + 1), score_pos(k + 2), score_pos(k + 3),
            //score_pos(k + 4), score_pos(k + 5), score_pos(k + 6), score_pos(k + 7));

            //TODO working: _scores = _mm512_setr_ps(valid_score(k), valid_score(k+1), valid_score(k+2), valid_score(k+3),
            //                  valid_score(k+4), valid_score(k+5), valid_score(k+6), valid_score(k+7),
            //                  valid_score(k+8), valid_score(k+9), valid_score(k+10), valid_score(k+11),
            //                  valid_score(k+12), valid_score(k+13), valid_score(k+14), valid_score(k+15));


            _result = _mm512_fmadd_ps(_scalar, _mini_vector, _scores);

            _mm512_i32scatter_ps(&dot_scores[0], _score_pos2, _result, 4);
            /*for (int l = 0; l < 16; ++l) {
                dot_scores[precursor_idx->get_rank(ion_bin[k + l].parent_id) - lower_rank] = _result[l];
            }*/

        }

        //Computing remainder
        for (; k < end_point; ++k) {
            fragment &f = ion_bin[k];
            dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
        }
    }

    /*
     * Prepare and rescore top matches
     */

    std::vector<int> order = order_of_scores(dot_scores);

    int num_required_rescores = std::max(2, settings::num_hit_ranks); //min 2 PSM to rescore (for delta dot)
    std::vector<match> top_matches;
    int duplicates = 0;
    for (int i = 0; i < (num_required_rescores + duplicates) && i < order.size(); ++i) {
        int elem_idx = order[i];
        int target_rank = elem_idx + lower_rank;
        unsigned int target_id = precursor_idx->get_precursor_by_rank(target_rank).id;

        //Determine if matched peptide is already in the list
        bool is_duplicate = false;
        for (match &m : top_matches) {
            std::string peptide = precursor_idx->get_precursor(target_id).peptide;
            if (is_peptide_isomer(peptide, precursor_idx->get_precursor(m.target_id).peptide)) {

                is_duplicate = true;
                if (peptide!=precursor_idx->get_precursor(m.target_id).peptide) { // Not identical (actually an isomer)
                    bool is_new_isomer = true;
                    for (std::string &isomer : m.isomers) {
                        if (isomer == peptide) { // Check if peptide is already in isomer list
                            is_new_isomer = false;
                            break;
                        }
                    }
                    if (is_new_isomer)
                        m.isomers.push_back(peptide); //keeping track of matched isomers
                }
                break;

            }
        }
        if (is_duplicate) {
            ++duplicates;
            continue;
        }


        // Create and rescore the current match
        match current_match = match(search_id, target_id);
        current_match.dot_product =  dot_scores[elem_idx];
        current_match.mass_difference = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
        current_match.abs_mass_difference = abs(current_match.mass_difference);
        current_match.ppm_difference = current_match.abs_mass_difference / (spec->precursor_mass / 1000000.f);

        rescore_match(current_match);
        top_matches.push_back(current_match);
    }


    // Record matches
    std::lock_guard<std::mutex> guard(pool->mtx);
    for (auto &psm : top_matches) {
        matches.push_back(psm);
    }
    --spec->search_counter;
    return true;
}
#elif USE_AVX_2
bool search_manager::search_spectrum(unsigned int search_id) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];

    // Determine range of candidate spectra
    float mz_tolerance = settings::mz_tolerance;
    if (settings::use_ppm_tolerance) {
        mz_tolerance = settings::ppm_factor * spec->precursor_mass;
    }
    int lower_rank = precursor_idx->get_lower_bound(spec->charge,spec->precursor_mass - mz_tolerance);
    int upper_rank = precursor_idx->get_upper_bound(spec->charge,spec->precursor_mass + mz_tolerance);

    if (lower_rank < 0 || upper_rank < 0 || lower_rank > upper_rank) {
        return false;
    }


    // Init candidate scores
    std::vector<float> dot_scores(upper_rank - lower_rank + 1, 0.f);

    // Update scores by matching all peaks using the fragment ion index
    for (int j = 0; j < spec->binned_peaks.size(); ++j) {

        // Open ion mz bin for corresponding peak
        fragment_bin &ion_bin = frag_idx->fragment_bins[spec->binned_peaks[j]];
        fragment_binn &bin = frag_idx->frag_bins[spec->binned_peaks[j]];


        // Determine starting point of lowest (candidate) parent index inside bin
        int starting_point_inside_bin = std::lower_bound(ion_bin.begin(), ion_bin.end(), lower_rank, [&](fragment &f, int rank) {
            return precursor_idx->get_rank(f.parent_id) < rank;
        }) - ion_bin.begin();
        int end_point = std::upper_bound(ion_bin.begin() + starting_point_inside_bin, ion_bin.end(), upper_rank, [&](int rank, fragment &f) {
            return precursor_idx->get_rank(f.parent_id) > rank;
        }) - ion_bin.begin();

        //Obtain dot-scores iff in valid range (0 otherwise
        auto valid_score = [&](int ii) {
            if (ii < end_point)
                return dot_scores[precursor_idx->get_rank(ion_bin[ii].parent_id) - lower_rank];
            return 0.f;
        };
        auto score_pos_secure = [&](int ii) {
            if (ii < end_point)
                return precursor_idx->get_rank(ion_bin[ii].parent_id) - lower_rank;
            return (unsigned) 0;
        };
        auto score_pos = [&](int ii) {
            return precursor_idx->get_rank(ion_bin[ii].parent_id);
        };



        /*
        __m256 _scalar = _mm256_set1_ps(spec->binned_intensities[j]);
        __m256i _lower_rank = _mm256_set1_epi32(lower_rank);
         */
        //Update scores for all parents with fragments in the range
        /*for (int k = starting_point_inside_bin; k < end_point; k+=8) {

            //Fill vectors with 8 float values
            __m256 _mini_vector = _mm256_loadu_ps(&bin.intensities[k]);
            __m256i _score_pos = _mm256_setr_epi32(score_pos_secure(k), score_pos_secure(k+1), score_pos_secure(k+2), score_pos_secure(k+3),
                                                  score_pos_secure(k+4), score_pos_secure(k+5), score_pos_secure(k+6), score_pos_secure(k+7));
            __m256 _scores = _mm256_i32gather_ps(&dot_scores[0], _score_pos, 4); //TODO why 4????

            //__m256 _scoresBefore = {valid_score(k), valid_score(k+1), valid_score(k+2), valid_score(k+3),
            //                     valid_score(k+4), valid_score(k+5), valid_score(k+6), valid_score(k+7)};


            __m256 _result = _mm256_fmadd_ps(_scalar, _mini_vector, _scores);

            //_mm256_i32scatter_ps(&dot_scores[0], _score_pos, _result, 4);
            for (int l = 0; l < 8 && k + l < end_point; ++l) {
                dot_scores[precursor_idx->get_rank(ion_bin[k + l].parent_id) - lower_rank] = _result[l];
            }
        }*/

        __m256 _scalar = _mm256_set1_ps(spec->binned_intensities[j]);
        __m256i _lower_rank = _mm256_set1_epi32(lower_rank);
        __m256 _mini_vector, _scores, _result;
        int k = starting_point_inside_bin;
        //First reach k mod 8 = 0
        for (; k % 8 != 0 && k < end_point; ++k) {
            fragment &f = ion_bin[k];
            dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
        }
        for (; k+8 < end_point; k+=8) {
            //Fill vectors with 8 float values
            _mini_vector = bin._intensities[k/8];
            //__m256 _mini_vector = _mm256_load_ps(&bin.intensities[k]);
            /*__m256i _score_pos0 = bin._parent_ids[k/8];
            __m256i _score_pos1 = _mm256_i32gather_epi32(&precursor_idx->to_rank[0], _score_pos0, 4); //TODO why 4????
            __m256i _score_pos2 = _mm256_sub_epi32(_score_pos1, _lower_rank);*/
            //__m256 _scores = _mm256_i32gather_ps(&dot_scores[0], _mm256_sub_epi32(bin._parent_ranks[k/8], _lower_rank), 4); //TODO why 4????

            //__m256i _score_pos1 = _mm256_setr_epi32(score_pos(k), score_pos(k + 1), score_pos(k + 2), score_pos(k + 3),
            //score_pos(k + 4), score_pos(k + 5), score_pos(k + 6), score_pos(k + 7));

            _scores = _mm256_setr_ps(valid_score(k), valid_score(k+1), valid_score(k+2), valid_score(k+3),
                              valid_score(k+4), valid_score(k+5), valid_score(k+6), valid_score(k+7));


            _result = _mm256_fmadd_ps(_scalar, _mini_vector, _scores);

            //_mm256_i32scatter_ps(&dot_scores[0], _score_pos, _result, 4);
            for (int l = 0; l < 8; ++l) {
                dot_scores[precursor_idx->get_rank(ion_bin[k + l].parent_id) - lower_rank] = _result[l];
            }

        }

        //Computing remainder
        for (; k < end_point; ++k) {
            fragment &f = ion_bin[k];
            dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
        }
    }

    /*
     * Prepare and rescore top matches
     */

    std::vector<int> order = order_of_scores(dot_scores);

    int num_required_rescores = std::max(2, settings::num_hit_ranks); //min 2 PSM to rescore (for delta dot)
    std::vector<match> top_matches;
    int duplicates = 0;
    for (int i = 0; i < (num_required_rescores + duplicates) && i < order.size(); ++i) {
        int elem_idx = order[i];
        int target_rank = elem_idx + lower_rank;
        unsigned int target_id = precursor_idx->get_precursor_by_rank(target_rank).id;

        //Determine if matched peptide is already in the list
        bool is_duplicate = false;
        for (match &m : top_matches) {
            std::string peptide = precursor_idx->get_precursor(target_id).peptide;
            if (is_peptide_isomer(peptide, precursor_idx->get_precursor(m.target_id).peptide)) {

                is_duplicate = true;
                if (peptide!=precursor_idx->get_precursor(m.target_id).peptide) { // Not identical (actually an isomer)
                    bool is_new_isomer = true;
                    for (std::string &isomer : m.isomers) {
                        if (isomer == peptide) { // Check if peptide is already in isomer list
                            is_new_isomer = false;
                            break;
                        }
                    }
                    if (is_new_isomer)
                        m.isomers.push_back(peptide); //keeping track of matched isomers
                }
                break;

            }
        }
        if (is_duplicate) {
            ++duplicates;
            continue;
        }


        // Create and rescore the current match
        match current_match = match(search_id, target_id);
        current_match.dot_product =  dot_scores[elem_idx];
        current_match.mass_difference = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
        current_match.abs_mass_difference = abs(current_match.mass_difference);
        current_match.ppm_difference = current_match.abs_mass_difference / (spec->precursor_mass / 1000000.f);

        rescore_match(current_match);
        top_matches.push_back(current_match);
    }


    // Record matches
    std::lock_guard<std::mutex> guard(pool->mtx);
    for (auto &psm : top_matches) {
        matches.push_back(psm);
    }
    --spec->search_counter;
    return true;
}

#else
bool search_manager::search_spectrum(unsigned int search_id) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];

    // Determine range of candidate spectra
    float mz_tolerance = settings::mz_tolerance;
    if (settings::use_ppm_tolerance) {
        mz_tolerance = settings::ppm_factor * spec->precursor_mass;
    }
    int lower_rank = precursor_idx->get_lower_bound(spec->charge,spec->precursor_mass - mz_tolerance);
    int upper_rank = precursor_idx->get_upper_bound(spec->charge,spec->precursor_mass + mz_tolerance);

    if (lower_rank < 0 || upper_rank < 0 || lower_rank > upper_rank) { // todo necessary? No matching precursor masses
        return false;
    }


    // Init candidate scores
    std::vector<float> dot_scores(upper_rank - lower_rank + 1, 0.f);

    // Update scores by matching all peaks using the fragment ion index
    for (int j = 0; j < spec->binned_peaks.size(); ++j) {

        // Open ion mz bin for corresponding peak
        fragment_bin &ion_bin = frag_idx->fragment_bins[spec->binned_peaks[j]];

        // Determine starting point of lowest (candidate) parent index inside bin
        int starting_point_inside_bin = std::lower_bound(ion_bin.begin(), ion_bin.end(), lower_rank, [&](fragment &f, int rank) {
            return precursor_idx->get_rank(f.parent_id) < rank;
        }) - ion_bin.begin();
        int end_point = std::upper_bound(ion_bin.begin() + starting_point_inside_bin, ion_bin.end(), upper_rank, [&](int rank, fragment &f) {
            return precursor_idx->get_rank(f.parent_id) > rank;
        }) - ion_bin.begin();

        //Update scores for all parents with fragments in the range
        for (int k = starting_point_inside_bin; k < end_point; ++k) {
            fragment &f = ion_bin[k];
            dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
        }

    }

    /*
     * Prepare and rescore top matches
     */

    std::vector<int> order = order_of_scores(dot_scores);

    int num_required_rescores = std::max(2, settings::num_hit_ranks); //min 2 PSM to rescore (for delta dot)
    std::vector<match> top_matches;
    int duplicates = 0;
    for (int i = 0; i < (num_required_rescores + duplicates) && i < order.size(); ++i) {
        int elem_idx = order[i];
        int target_rank = elem_idx + lower_rank;
        unsigned int target_id = precursor_idx->get_precursor_by_rank(target_rank).id;

        //Determine if matched peptide is already in the list
        bool is_duplicate = false;
        for (match &m : top_matches) {
            std::string peptide = precursor_idx->get_precursor(target_id).peptide;
            if (is_peptide_isomer(peptide, precursor_idx->get_precursor(m.target_id).peptide)) {

                is_duplicate = true;
                if (peptide!=precursor_idx->get_precursor(m.target_id).peptide) { // Not identical (actually an isomer)
                    bool is_new_isomer = true;
                    for (std::string &isomer : m.isomers) {
                        if (isomer == peptide) { // Check if peptide is already in isomer list
                            is_new_isomer = false;
                            break;
                        }
                    }
                    if (is_new_isomer)
                        m.isomers.push_back(peptide); //keeping track of matched isomers
                }
                break;

            }
        }
        if (is_duplicate) {
            ++duplicates;
            continue;
        }


        // Create and rescore the current match
        match current_match = match(search_id, target_id);
        current_match.dot_product =  dot_scores[elem_idx];
        current_match.mass_difference = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
        current_match.abs_mass_difference = abs(current_match.mass_difference);
        current_match.ppm_difference = current_match.abs_mass_difference / (spec->precursor_mass / 1000000.f);

        rescore_match(current_match);
        top_matches.push_back(current_match);
    }


    // Record matches
    std::lock_guard<std::mutex> guard(pool->mtx);
    for (auto &psm : top_matches) {
        matches.push_back(psm);
    }
    --spec->search_counter;
    return true;
}
#endif

bool search_manager::merge_matches() {

    /*
     * Sort by id ascending and then score descending
     */

    //std::cout << matches.size() << std::endl;
    std::sort(matches.begin(), matches.end(), [](match &a, match &b) {
        return (a.query_id == b.query_id && a.dot_product > b.dot_product) || a.query_id < b.query_id;
    });

    /*
     * Sort by id and determine delta scores, then sort by discriminant
     */
    float spectraST_minimum = 0.00001f; //positive and negative (for invalid scores)
    auto start = matches.begin();
    for (auto iter = matches.begin(); iter != matches.end() + 1; ++iter) {
        if (iter == matches.end() || iter->query_id != start->query_id) {
            auto end = iter - 1;

            //Assess dot difference
            float reference_dot = 0.f;
            if (start != end) {
                reference_dot = (start + 1)->dot_product;
            }
            for (auto iiter=start; iiter != (end+1); ++iiter) {
                iiter->delta_dot = iiter->dot_product - reference_dot;
                if (iiter->dot_product < spectraST_minimum) {
                    iiter->spectraST_score_dot = -spectraST_minimum;
                    continue;
                }
                float delta_norm = std::min(iiter->delta_dot / iiter->dot_product, 2.f * iiter->dot_product);
                iiter->spectraST_score_dot += 0.4f * delta_norm;
                iiter->spectraST_score_dot = std::max(iiter->spectraST_score_dot, -spectraST_minimum);
            }


            //Assess sim difference
            std::sort(start, end+1, [](match &a, match &b) {
                return a.similarity > b.similarity;
            });
            float reference_sim = 0.f;
            if (start != end) {
                reference_sim = (start + 1)->similarity; //2nd of the same query id
            }
            for (auto iiter=start; iiter != end+1; ++iiter) {
                iiter->delta_similarity = iiter->similarity - reference_sim;
                if (iiter->similarity < spectraST_minimum) {
                    iiter->spectraST_score = -spectraST_minimum;
                    continue;
                }
                float delta_norm = std::min(iiter->delta_similarity / iiter->similarity, 2.f * iiter->similarity);
                iiter->spectraST_score += 0.4f * delta_norm;
                iiter->spectraST_score = std::max(iiter->spectraST_score, -spectraST_minimum);
            }

            //Sort by sim2 and calculate sim2 dif
            std::sort(start, end+1, [](match &a, match &b) {
                return a.sim2 > b.sim2;
            });

            float reference_sim2 = 0.f;
            if (start != end) {
                reference_sim2 = (start + 1)->sim2; //2nd of the same query id
            }
            for (auto iiter=start; iiter != end+1; ++iiter) {
                iiter->delta_sim2 = iiter->sim2 - reference_sim2;
            }

            //Sort annotation sim
            std::sort(start, end+1, [](match &a, match &b) {
                return a.annotation_similarity > b.annotation_similarity;
            });

            float reference_annotation_sim = 0.f;
            if (start != end) {
                reference_annotation_sim = (start + 1)->annotation_similarity; //2nd of the same query id
            }
            for (auto iiter=start; iiter != end+1; ++iiter) {
                iiter->delta_annotation_sim = iiter->annotation_similarity - reference_annotation_sim;
            }

            // Finally
            // Sort by discriminant scoring metric (here: avg bias adjusted similarity, if equal then standard similarity))
            std::sort(start, end+1, [](match &a, match &b) {
                return a.avg_bias_adj_similarity > b.avg_bias_adj_similarity || (a.avg_bias_adj_similarity == b.avg_bias_adj_similarity && a.similarity > b.similarity);
            });

            float reference_avg = 0.f;
            if (start != end) {
                reference_avg = (start + 1)->avg_bias_adj_similarity; //2nd of the same query id
            }

            int rank = 1;
            for (auto iiter=start; iiter != end+1; ++iiter) {
                iiter->delta_avg = iiter->avg_bias_adj_similarity - reference_avg;
                iiter->hit_rank = rank;
                ++rank;
            }
            start = iter;
        }
    }

    //std::cout << matches.size() << std::endl;
    return true;
}

bool search_manager::save_search_results_to_file(const std::string &file_path) {
    std::ofstream outfile;
    std::string delim = "\t";

    outfile.open(file_path, std::ios::out);
    if (!outfile.good())
        return false;

    // Add header
    if (settings::save_search_command) {
        outfile << "#" << settings::search_command << "\n";
    }
    outfile << "id" + delim + "spectrum" + delim + "charge" + delim + "hit_rank" + delim + "match" + delim + "peptide" + delim + "isomers" + delim + "similarity" + delim + "bias" + delim + "annotation_similarity"+ delim + "annotation_bias" + delim + "avg_bias_adjusted_similarity" + delim + "dot_product" + delim + "delta_dot" + delim + "delta_similarity" + delim + "delta_annotaion_similarity" + delim + "delta_sim2" + delim + "delta_avg" + delim + "dot_contrast_angle" + delim + "similarity_contrast_angle" + delim + "annotation_contrast_angle" + delim + "mass_difference" + delim + "abs_mass_difference" + delim + "ppm_difference" + delim + "peptide_length" + delim + "precursor_mz" + delim + "peak_count_query" + delim + "peak_count_ref" + delim + "sim2" + delim + "annotation_sim2" + delim + "x_score" + delim + "x_score_dot" + delim + "x_lgamma" + delim + "x_lgamma_dot" + delim + "st_score" + delim + "st_score_dot\n";

    // Go through matches and parse relevant information for each
    for (int i = 0; i < matches.size(); ++i) {
        match &psm = matches[i];
        if (psm.hit_rank <= settings::num_hit_ranks) {
            precursor &target = precursor_idx->get_precursor(psm.target_id);
            std::string name = search_library.spectrum_list[psm.query_id]->name;
            std::string id = name + "/" + std::to_string(psm.hit_rank);
            std::string iso;
            for (std::string &s : psm.isomers) {
                iso += s + ";";
            }
            if (!iso.empty())
                iso.pop_back();
            outfile << id << delim << name << delim << psm.charge << delim << psm.hit_rank << delim << target.id << delim << target.peptide << delim << iso << delim << psm.similarity << delim << psm.bias << delim << psm.annotation_similarity << delim << psm.annotation_bias << delim << psm.avg_bias_adj_similarity << delim << psm.dot_product << delim << psm.delta_dot << delim << psm.delta_similarity << delim << psm.delta_annotation_sim << delim << psm.delta_sim2 << delim << psm.delta_avg << delim << psm.dot_contrast_angle << delim << psm.similarity_contrast_angle << delim << psm.annotation_contrast_angle << delim << psm.mass_difference << delim << psm.abs_mass_difference << delim << psm.ppm_difference << delim << target.peptide.length() << delim << search_library.spectrum_list[psm.query_id]->precursor_mass << delim << psm.peak_count_query << delim << psm.peak_count_target << delim << psm.sim2 << delim << psm.annotation_sim2 << delim << psm.x_hunter_score << delim << psm.x_hunter_score_dot << delim << psm.x_lgamma << delim << psm.x_lgamma_dot << delim << psm.spectraST_score << delim << psm.spectraST_score_dot << "\n";
        }

    }

    outfile.close();
    return true;
}


bool search_manager::save_search_results_in_pin_format(const std::string &file_path) {
    std::ofstream outfile;
    std::string delim = "\t";

    outfile.open(file_path, std::ios::out);
    if (!outfile.good())
        return false;

    // Add header
    //Percolator format SpecId	Label	ScanNr	lnrSp	deltLCn	deltCn	Xcorr	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	enzN	enzC	enzInt	lnNumSP	dM	absdM	Peptide	Proteins
    outfile << "PSMId" + delim + "Label" + delim + "ScanNr" + delim + "charge" + delim + "similarity" + delim + "bias" + delim + "annotation_similarity"+ delim + "annotation_bias" + delim + "avg_bias_adjusted_similarity" + delim + "dot_product" + delim + "delta_dot" + delim + "delta_similarity" + delim + "delta_annotation_similarity" + delim + "delta_sim2" + delim + "delta_avg" + delim + "dot_contrast_angle" + delim + "similarity_contrast_angle" + delim + "annotation_contrast_angle" + delim + "mass_difference" + delim + "abs_mass_difference" + delim + "ppm_difference" + delim + "peptide_length" + delim + "precursor_mz" + delim + "peak_count_query" + delim + "peak_count_ref" + delim + "sim2" + delim + "annotation_sim2" + delim + "x_score" + delim + "x_score_dot" + delim + "x_lgamma" + delim + "x_lgamma_dot" + delim + "st_score" + delim + "st_score_dot" + delim + "Peptide" + delim + "Proteins" + "\n";

    // Go through matches and parse relevant information for each
    for (int i = 0; i < matches.size(); ++i) {
        match &psm = matches[i];
        if (psm.hit_rank <= 1) { // Only output rank-1 psms for percolator
            precursor &target = precursor_idx->get_precursor(psm.target_id);
            std::string name = search_library.spectrum_list[psm.query_id]->name;
            std::string id = name + "/" + std::to_string(psm.hit_rank);
            std::string name_cropped = name.substr(0, name.rfind('.'));
            std::string scannum = name_cropped.substr(name_cropped.rfind('.') + 1, std::string::npos);
            std::string peptide = "X." + target.peptide + ".X"; //Placeholders for flanking amino acids (needed in pin-tab format)
            std::string iso;
            for (std::string &s : psm.isomers) {
                iso += s + ";";
            }
            if (!iso.empty())
                iso.pop_back();
            outfile << i << delim << config->label << delim << scannum << delim << psm.charge << delim << psm.similarity << delim << psm.bias << delim << psm.annotation_similarity << delim << psm.annotation_bias << delim << psm.avg_bias_adj_similarity << delim << psm.dot_product << delim << psm.delta_dot << delim << psm.delta_similarity << delim << psm.delta_annotation_sim << delim << psm.delta_sim2 << delim << psm.delta_avg << delim << psm.dot_contrast_angle << delim << psm.similarity_contrast_angle << delim << psm.annotation_contrast_angle << delim << psm.mass_difference << delim << psm.abs_mass_difference << delim << psm.ppm_difference << delim << target.peptide.length() << delim << search_library.spectrum_list[psm.query_id]->precursor_mass << delim << psm.peak_count_query << delim << psm.peak_count_target << delim << psm.sim2 << delim << psm.annotation_sim2 << delim << psm.x_hunter_score << delim << psm.x_hunter_score_dot << delim << psm.x_lgamma << delim << psm.x_lgamma_dot << delim << psm.spectraST_score << delim << psm.spectraST_score_dot << delim << peptide << delim << "Unknown" << "\n";
        }

    }

    outfile.close();
    return true;
}

//__m256 _mini_vector = {ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity}; //_mm256_set_ps(ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity);//_mm256_load_ps(&vec[i]);


float search_manager::rescore_spectrum(unsigned int search_id, unsigned int target_id) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];

    float score = 0.f;

    for (int i = 0; i < spec->peak_positions.size(); ++i) {
        float mz = spec->peak_positions[i];
        float intensity = spec->intensities[i];

        /*
         * Extract all matching peaks within a range of +-5 sigma
         */

        std::vector<std::pair<float, float>> peaks;
        int lower_bin = spectrum::get_mz_bin(mz - 5 * sigma);
        int upper_bin = spectrum::get_mz_bin(mz + 5 * sigma);

        for (int bin = lower_bin; bin <= upper_bin; ++bin) {
            if (bin < 0 || bin >= spec->num_bins) {
                continue;
            }
            fragment_bin &ion_bin = frag_idx->fragment_bins[bin];
            if (ion_bin.empty())
                continue;

            fragment &f = *std::lower_bound(ion_bin.begin(), ion_bin.end(), precursor_idx->get_rank(target_id), [&](fragment &f, int rank) {
                return precursor_idx->get_rank(f.parent_id) < rank;
            });

            if (f.parent_id != target_id) {
                continue;
            }


            if (!f.peak_composition.empty()) {
                for (std::pair<float, float> p : f.peak_composition) {
                    peaks.push_back(p);
                }
            } else {
                peaks.emplace_back(f.mz, f.intensity);
            }

        }

        /*
         * Score spectrum peak to all peaks of target using normal distribution for intensity fall-off
         */

        for (auto &peak : peaks) {
            float distance = mz - peak.first;
            float normal_factor = normal_pdf(distance, 0, sigma) / max_normal;

            score += intensity * peak.second * normal_factor;
        }
    }

    return score;






}


bool search_manager::rescore_match_old(match &psm) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[psm.query_id];

    float score = 0.f;
    int peak_count_query = 0;
    std::set<float> matched_peaks; //on target side

    for (int i = 0; i < spec->peak_positions.size(); ++i) {
        float mz = spec->peak_positions[i];
        float intensity = spec->intensities[i];

        /*
         * Extract all matching peaks within a range of +-5 sigma
         */

        std::vector<std::pair<float, float>> peaks;
        int lower_bin = spectrum::get_mz_bin(mz - 5 * sigma);
        int upper_bin = spectrum::get_mz_bin(mz + 5 * sigma);

        for (int bin = lower_bin; bin <= upper_bin; ++bin) {
            if (bin < 0 || bin >= spec->num_bins) {
                continue;
            }
            fragment_bin &ion_bin = frag_idx->fragment_bins[bin];
            if (ion_bin.empty())
                continue;

            fragment &f = *std::lower_bound(ion_bin.begin(), ion_bin.end(), precursor_idx->get_rank(psm.target_id), [&](fragment &f, int rank) {
                return precursor_idx->get_rank(f.parent_id) < rank;
            });

            if (f.parent_id != psm.target_id) {
                continue;
            }


            if (!f.peak_composition.empty()) {
                for (std::pair<float, float> p : f.peak_composition) {
                    peaks.push_back(p);
                }
            } else {
                peaks.emplace_back(f.mz, f.intensity);
            }

        }

        /*
         * Score spectrum peak to all peaks of target using normal distribution for intensity fall-off
         */

        bool peak_matched = false;
        for (auto &peak : peaks) {
            float distance = mz - peak.first;
            float normal_factor = normal_pdf(distance, 0, sigma) / max_normal;

            score += intensity * peak.second * normal_factor;

            if (distance <= sigma) {
                matched_peaks.insert(peak.first);
                peak_matched = true;
            }
        }
        if (peak_matched) {
            ++peak_count_query;
        }
    }


    /*
     * Update Match
     */

    psm.peak_count_query = peak_count_query;
    psm.peak_count_target = matched_peaks.size();
    psm.similarity = score;

    return true;
}

bool search_manager::rescore_match(match &psm) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[psm.query_id];

    float score = 0.f;
    float bias = 0.f;
    float min_counter_score = 0.f;
    std::set<std::pair<float, float>> target_peaks;

    /*
     * Extract all matching peaks within a range of +-5 sigma
     */

    for (int i = 0; i < spec->peak_positions.size(); ++i) {
        float mz = spec->peak_positions[i];


        std::vector<std::pair<float, float>> peaks;
        int lower_bin = spectrum::get_mz_bin(mz - 5 * sigma);
        int upper_bin = spectrum::get_mz_bin(mz + 5 * sigma);

        for (int bin = lower_bin; bin <= upper_bin; ++bin) {
            if (bin < 0 || bin >= spec->num_bins) {
                continue;
            }
            fragment_bin &ion_bin = frag_idx->fragment_bins[bin];


            if (ion_bin.empty()) //TODO maybe not necessary due to next if (!= end())
                continue;

            auto lower = std::lower_bound(ion_bin.begin(), ion_bin.end(), precursor_idx->get_rank(psm.target_id), [&](fragment &f, int rank) {
                return precursor_idx->get_rank(f.parent_id) < rank;
            });
            if (lower == ion_bin.end())
                continue;

            fragment &f = *lower;
            if (f.parent_id != psm.target_id)
                continue;


            if (!f.peak_composition.empty()) {
                for (std::pair<float, float> p : f.peak_composition) {
                    target_peaks.insert(p);
                }
            } else {
                target_peaks.insert(std::make_pair(f.mz, f.intensity));
            }

        }


    }

    /*
     * Score peaks (max score per peak match)
     *
     * From Target (reference) side
     */

    std::vector<std::pair<float, float>> paired_query_peaks(target_peaks.size());
    std::vector<bool> is_paired(target_peaks.size(), false);
    int i = 0;
    int peak_count_ref = 0;
    for (auto &peak : target_peaks) {
        float mz = peak.first;
        float intensity = peak.second;
        float peak_score = 0.f;
        bool peak_matched = false;

        for (int j = 0; j < spec->peak_positions.size(); ++j) {
            float s_mz = spec->peak_positions[j];

            float distance = mz - s_mz;
            if (abs(distance) < 5 * sigma) {

                float normal_factor = normal_pdf_scaled(distance, 0, sigma);
                float new_score = intensity * spec->intensities[j] * normal_factor;
                if (new_score > peak_score) {
                    paired_query_peaks[i] = std::make_pair(spec->peak_positions[j], spec->intensities[j]);
                    is_paired[i] = true;
                    peak_score = new_score;
                }
                if (abs(distance) < sigma && new_score >= min_counter_score) {
                    peak_matched = true;
                }

            }
        }

        score += peak_score;
        bias += (peak_score * peak_score); //TODO dot bias
        if (peak_matched) {
            ++peak_count_ref;
        }
        ++i;
    }


    /*
     * Prosit b/y-ion annotation score (disregarding all other peaks)
     */

    //Re-Normalize peak intensities (this is all that changes)
    spectrum::normalize_intensity_vector(paired_query_peaks);

    i = 0;
    float annotation_score = 0.f;
    float annotation_bias = 0.f;

    for (auto &peak : target_peaks) {
        if (is_paired[i]) {
            float mz = peak.first;
            float intensity = peak.second;
            float distance = mz - paired_query_peaks[i].first;

            float normal_factor = normal_pdf_scaled(distance, 0, sigma);
            float peak_score = intensity * paired_query_peaks[i].second * normal_factor;

            annotation_score += peak_score;
            annotation_bias += (peak_score * peak_score); //TODO dot bias
        }
        ++i;
    }


    /*
     * Update Match
     */

    psm.charge = spec->charge;
    psm.peak_count_query = 1000;
    psm.peak_count_target = peak_count_ref;
    psm.similarity = score;
    psm.annotation_similarity = annotation_score;
    if (score > 0) {
        psm.bias = std::sqrt(bias) / score;
        psm.annotation_bias = std::sqrt(annotation_bias) / annotation_score;
    }
    else {
        psm.bias = 0.f;
        psm.annotation_bias = 0.f;
        psm.annotation_sim2 = 0.f;
    }

    // Advanced scores
    psm.sim2 = psm.similarity * (1.f - psm.bias);
    psm.annotation_sim2 = psm.annotation_similarity * (1.f - psm.annotation_bias);
    psm.avg_bias_adj_similarity = (psm.sim2 + psm.annotation_sim2) / 2.f;

    auto scp_factorial = (double) factorial(psm.peak_count_target);

    psm.x_hunter_score = psm.similarity * scp_factorial;
    psm.x_hunter_score_dot = psm.dot_product * scp_factorial;

    // Spectral contrast angles
    psm.similarity_contrast_angle = contrast_angle(psm.similarity);
    psm.dot_contrast_angle = contrast_angle(psm.dot_product);
    psm.annotation_contrast_angle = contrast_angle(psm.annotation_similarity);


    /*
     * lgamma of the x!score
     */
    if (psm.similarity < 0.00001f || psm.dot_product < 0.00001f || psm.peak_count_target == 0) { //catch log(0)
        psm.x_lgamma = -100.f;
        psm.x_lgamma_dot = -100.f;
    } else {
        psm.x_lgamma = log((double) psm.similarity) + lgamma(psm.peak_count_target);
        psm.x_lgamma_dot = log((double) psm.dot_product) + lgamma(psm.peak_count_target);
    }

    auto bias_penalty = [](float bias) {
        if (bias < 0.1 || bias > 0.35 && bias <= 0.4)
            return 0.12f;
        else if(bias > 0.4 && bias <= 0.45)
            return 0.18f;
        else if (bias > 0.45) {
            return 0.24f;
        }
        return 0.f;
    };
    psm.spectraST_score_dot = 0.6f * psm.dot_product - bias_penalty(psm.bias); //TODO this should be dot bias
    psm.spectraST_score = 0.6f * psm.similarity - bias_penalty(psm.bias); //delta sim is added when merging
    psm.hit_rank = 0;

    return true;
}

long search_manager::get_time_spent_in_inner_search() {
    return std::chrono::duration_cast<std::chrono::seconds>(inner_search_duration).count();
}

float search_manager::normal_pdf(float x, float mean, float standard_deviation) {
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float x_deviation = (x - mean) / standard_deviation;

    return inv_sqrt_2pi / standard_deviation * std::exp(-0.5f * x_deviation * x_deviation);
}

float search_manager::normal_pdf_scaled(float x, float mean, float standard_deviation) {
    float x_deviation = (x - mean) / standard_deviation;
    return std::exp(-0.5f * x_deviation * x_deviation);
}

float search_manager::contrast_angle(float dot) {
    return 1 - (2 * std::acos(dot)) / M_PI;
}

long long unsigned int search_manager::factorial(int n) {
    long long unsigned int fact = 1;
    for (int i = 1; i <= n; ++i) {
        fact *= i;
    }
    return fact;
}

std::vector<int> search_manager::order_of_scores(std::vector<float> &scores) {
    std::vector<int> indices(scores.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(), [&scores](int a, int b) {
        return scores[a] > scores[b];
    });

    return indices;
}

bool search_manager::is_peptide_isomer(std::string &peptide, std::string &other) {
    if (peptide.length() != other.length())
        return false;
    for (int i = 0; i < peptide.length(); ++i) {
        if (peptide[i] != other[i]) {
            if (peptide[i] == 'L' && other[i] == 'I' || peptide[i] == 'I' && other[i] == 'L') {
                continue;
            }
            return false;
        }
    }

    return true;
}

