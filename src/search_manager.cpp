#include <iostream>
#include <fstream>
#include <numeric>
#include <x86intrin.h>
//#include <x86intrin.h>
//#include <immintrin.h>

#include "search_manager.h"

search_manager::search_manager(std::string search_file_path, std::string index_directory_path) : search_file_path(search_file_path), index_directory_path(index_directory_path) {
    std::cout << "Calling MS2 recruiter" << std::endl;

    std::cout << "Configuring ... " << std::endl;
    config = std::make_shared<configuration>();
    config->load_configuration_from_file(index_directory_path + "config.txt");
    pool = std::make_shared<thread_pool>(4);
}


bool search_manager::prepare_search_library() {

    // Load search library
    search_library = library(search_file_path);

    // Map ms2 ids to search container (corresponding to sub-index)
    mapped_search_ids = std::vector<std::vector<unsigned int>>(config->num_indices);

    if (config->num_indices == 1) {
        mapped_search_ids[0].resize(search_library.spectrum_list.size());
        std::iota(mapped_search_ids[0].begin(), mapped_search_ids[0].end(), 0);
        return true;
    }

    for (unsigned int i = 0; i < search_library.spectrum_list.size(); ++i) {

        float mz = search_library.spectrum_list[i]->precursor_mass;
        float min_mz = mz - mz_tolerance;
        float max_mz = mz + mz_tolerance;

        // Map id (i) to every sub-index where min/max bounds fall into
        for (int idx_num = 0; idx_num < (config->num_indices - 1); ++idx_num) {
            if (min_mz <= config->sub_idx_limits[idx_num]) {
                if (idx_num == 0 || max_mz >= config->sub_idx_limits[idx_num - 1]) {
                    mapped_search_ids[idx_num].push_back(i);
                } else {
                    break; // because max_mz is below lower border of that and every subsequent sub-indices todo >(equality)< not needed for both I think, just there to be sure
                }
            }
        }
        if (max_mz >= config->sub_idx_limits.back()) { // check if it falls into last sub index (without upper border)
            mapped_search_ids.back().push_back(i);
        }

    }

    return true;
}

bool search_manager::prepare_precursor_index() {
    precursor_idx = std::make_shared<precursor_index>();
    precursor_idx->load_index_from_file(config->precursor_index_path);
    return true;
}

bool search_manager::perform_searches() {

    frag_idx = std::make_shared<fragment_ion_index>();

    for (int i = 0; i < config->num_indices; ++i) {
        std::cout << "Loading index number " << i << std::endl;
        frag_idx->load_index_from_binary_file(config->sub_idx_file_names[i]);
        frag_idx->update_intensities();
        //TODO set precursor index limits by subindex borders... has to be properly implemented
        std::cout << "Searching ... " << std::endl;
        std::vector<unsigned int> &search_ids = mapped_search_ids[i];

        for (unsigned int &s_id : search_ids) {
            //search_spectrum(s_id);
            //search_spectrum_avx(s_id);
            search_spectrum_avx2(s_id);
        }
    }
    return true;
}


bool search_manager::perform_searches_parallel() {
    frag_idx = std::make_shared<fragment_ion_index>();

    for (int i = 0; i < config->num_indices; ++i) {
        std::cout << "Loading index number " << i << std::endl;
        frag_idx->load_index_from_binary_file(config->sub_idx_file_names[i]);
        //TODO set precursor index limits by subindex borders... has to be properly implemented
        std::cout << "Searching ... " << std::endl;
        std::vector<unsigned int> &search_ids = mapped_search_ids[i];
        for (unsigned int &s_id : search_ids) {
            std::shared_ptr<spectrum> spec = search_library.spectrum_list[s_id];

            //auto f = [&](int s) {test(s);};
            //pool->enqueue([f, s_id] { return f(s_id); });

            pool->enqueue([this, s_id] () {task_search_spectrum(s_id);}); //TODO THIS WORKS
        }
        pool->wait_for_all_threads();
    }
    return true;
}

bool search_manager::search_spectrum(unsigned int search_id) {

    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];


    // Determine range of candidate spectra
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
            //std::cout << k << std::endl;
        }
        //std::cout << end_point << std::endl;
        //exit(12);
    }

    // Prepare best-scoring PSM
    int max_elem = max_element(dot_scores.begin(), dot_scores.end()) - dot_scores.begin();
    int target_rank = max_elem + lower_rank;
    float dot = dot_scores[max_elem];

    float mass_diff = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
    // Record match
    matches.emplace_back(match(search_id, precursor_idx->get_precursor_by_rank(target_rank).id, dot, mass_diff, 1));

    return true;
}

bool search_manager::task_search_spectrum(unsigned int search_id) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];
    // Determine range of candidate spectra
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

    // Prepare best-scoring PSM
    int max_elem = max_element(dot_scores.begin(), dot_scores.end()) - dot_scores.begin();
    int target_rank = max_elem + lower_rank;
    float dot = dot_scores[max_elem];

    float mass_diff = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
    // Record match

    std::lock_guard<std::mutex> guard(pool->mtx);
    matches.emplace_back(search_id, precursor_idx->get_precursor_by_rank(target_rank).id, dot, mass_diff, 1);
    return true;

}

bool search_manager::merge_matches() {

    /*
     * Sort by id ascending and then score descending
     */
    std::cout << matches.size() << std::endl;
    std::sort(matches.begin(), matches.end(), [](match a, match b) {
        return (a.query_id == b.query_id && a.dot_product > b.dot_product) || a.query_id < b.query_id;
    });

    /*
     * Remove subsequent entries with same query id (allowed by ordering based on id and score)
     */
    auto it = matches.begin() + 1;
    while(it != matches.end()) {
        if(it->query_id == (it - 1)->query_id) {
            it = matches.erase(it);
        }
        else ++it;
    }
    std::cout << matches.size() << std::endl;
    return true;
}

bool search_manager::save_search_results_to_file(const std::string &file_path) {
    std::ofstream outfile;
    std::string delimiter = "\t";


    outfile.open(file_path, std::ios::out);
    if (!outfile.good())
        return false;

    // Add header
    outfile << "spectrum"+delimiter+"match"+delimiter+"peptide"+delimiter+"dot-product"+delimiter+"mass-difference\n";

    // Go through matches and parse relevant information for each
    for (int i = 0; i < matches.size(); ++i) {
        match psm = matches[i];
        precursor &target = precursor_idx->get_precursor(psm.target_id);
        outfile << search_library.spectrum_list[psm.query_id]->name << delimiter << target.id << delimiter << target.peptide << delimiter << psm.dot_product << delimiter << psm.mass_difference << "\n";
    }

    outfile.close();
    return true;
}


bool search_manager::search_spectrum_avx(unsigned int search_id) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];


    // Determine range of candidate spectra
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
        int starting_point_inside_bin = std::lower_bound(ion_bin.begin(), ion_bin.end(), lower_rank, [&](fragment f, int rank) {
            return precursor_idx->get_rank(f.parent_id) < rank;
        }) - ion_bin.begin();



        __m256 _scalar = _mm256_set1_ps(spec->binned_intensities[j]);
        //float res[8];

        //Update scores for all parents with fragments in the range
        for (int k = starting_point_inside_bin; k < ion_bin.size() && precursor_idx->get_rank(ion_bin[k].parent_id) <= upper_rank; k+=8) {

            //Fill vector with 8 float values
            //__m256 _mini_vector = _mm256_setr_ps(ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity);
            __m256 _mini_vector = {ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity}; //_mm256_set_ps(ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity);//_mm256_load_ps(&vec[i]);
            __m256 _result = _mm256_mul_ps(_scalar, _mini_vector);
            //_mm256_store_ps(res, _result);
            for (int l = 0; (l < 8) && (k + l < ion_bin.size()) && (precursor_idx->get_rank(ion_bin[k + l].parent_id) <= upper_rank); ++l) {
                dot_scores[precursor_idx->get_rank(ion_bin[k + l].parent_id) - lower_rank] += _result[l];
            }
            //std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " << res[4] << " " << res[5] << " " << res[6] << " " << res[7] << std::endl;
            //break;
            //dot_scores[precursor_idx->get_rank(ion_bin[k].parent_id) - lower_rank] += res[0];
            /*dot_scores[precursor_idx->get_rank(ion_bin[k+1].parent_id) - lower_rank] += res[1];
            dot_scores[precursor_idx->get_rank(ion_bin[k+2].parent_id) - lower_rank] += res[2];
            dot_scores[precursor_idx->get_rank(ion_bin[k+3].parent_id) - lower_rank] += res[3];
            dot_scores[precursor_idx->get_rank(ion_bin[k+4].parent_id) - lower_rank] += res[4];
            dot_scores[precursor_idx->get_rank(ion_bin[k+5].parent_id) - lower_rank] += res[5];
            dot_scores[precursor_idx->get_rank(ion_bin[k+6].parent_id) - lower_rank] += res[6];
            dot_scores[precursor_idx->get_rank(ion_bin[k+7].parent_id) - lower_rank] += res[7];
            //fragment &f = ion_bin[k];
            //dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];*/
        }
        //Update scores for all parents with fragments in the range
        /*for (int k = starting_point_inside_bin; k < ion_bin.size() && precursor_idx->get_rank(ion_bin[k].parent_id) <= upper_rank; ++k) {
            fragment &f = ion_bin[k]; //todo SIMD
            dot_scores_old[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
            std::cout << f.intensity * spec->binned_intensities[j] << " " << std::flush;
            if (k == starting_point_inside_bin + 8)
                exit(12);
        }

        //TODO compare dot with old

        exit(12);
        */
    }

    // Prepare best-scoring PSM
    int max_elem = max_element(dot_scores.begin(), dot_scores.end()) - dot_scores.begin();
    int target_rank = max_elem + lower_rank;
    float dot = dot_scores[max_elem];

    float mass_diff = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
    // Record match
    matches.emplace_back(match(search_id, precursor_idx->get_precursor_by_rank(target_rank).id, dot, mass_diff, 1));

    return true;

}

bool search_manager::search_spectrum_avx2(unsigned int search_id) {
    std::shared_ptr<spectrum> spec = search_library.spectrum_list[search_id];


    // Determine range of candidate spectra
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

        //Update scores for all parents with fragments in the range
        __m256 _scalar = _mm256_set1_ps(spec->binned_intensities[j]);
        for (int k = starting_point_inside_bin; k < end_point; k+=8) {

            //Fill vectors with 8 float values
            __m256 _mini_vector = _mm256_loadu_ps(&bin.intensities[k]);
            __m256 _scores = {valid_score(k), valid_score(k+1), valid_score(k+2), valid_score(k+3),
                              valid_score(k+4), valid_score(k+5), valid_score(k+6), valid_score(k+7)};
            __m256 _result = _mm256_fmadd_ps(_scalar, _mini_vector, _scores);

            for (int l = 0; l < 8 && k + l < end_point; ++l) {
                dot_scores[precursor_idx->get_rank(ion_bin[k + l].parent_id) - lower_rank] = _result[l];
            }
        }
    }

    // Prepare best-scoring PSM
    int max_elem = max_element(dot_scores.begin(), dot_scores.end()) - dot_scores.begin();
    int target_rank = max_elem + lower_rank;
    float dot = dot_scores[max_elem];

    float mass_diff = precursor_idx->get_precursor_by_rank(target_rank).mz - spec->precursor_mass;
    // Record match
    matches.emplace_back(match(search_id, precursor_idx->get_precursor_by_rank(target_rank).id, dot, mass_diff, 1));

    return true;
}

//__m256 _mini_vector = {ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity}; //_mm256_set_ps(ion_bin[k].intensity, ion_bin[k+1].intensity, ion_bin[k+2].intensity, ion_bin[k+3].intensity, ion_bin[k+4].intensity, ion_bin[k+5].intensity, ion_bin[k+6].intensity, ion_bin[k+7].intensity);//_mm256_load_ps(&vec[i]);
