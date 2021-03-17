#include <iostream>
#include <fstream>
#include <numeric>
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
        //TODO set precursor index limits by subindex borders... has to be properly implemented
        std::cout << "Searching ... " << std::endl;
        std::vector<unsigned int> &search_ids = mapped_search_ids[i];

        for (unsigned int &s_id : search_ids) {
            std::shared_ptr<spectrum> spec = search_library.spectrum_list[s_id];
            search_spectrum(s_id, spec);
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

bool search_manager::search_spectrum(unsigned int search_id, std::shared_ptr<spectrum> &spec) {

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
        fragment_bin ion_bin = frag_idx->fragment_bins[spec->binned_peaks[j]];

        // Determine starting point of lowest (candidate) parent index inside bin
        int starting_point_inside_bin = std::lower_bound(ion_bin.begin(), ion_bin.end(), lower_rank, [&](fragment f, int rank) {
            return precursor_idx->get_rank(f.parent_id) < rank;
        }) - ion_bin.begin();

        //Update scores for all parents with fragments in the range
        for (int k = starting_point_inside_bin; k < ion_bin.size() && precursor_idx->get_rank(ion_bin[k].parent_id) <= upper_rank; ++k) {
            fragment &f = ion_bin[k]; //todo SIMD
            dot_scores[precursor_idx->get_rank(f.parent_id) - lower_rank] += f.intensity * spec->binned_intensities[j];
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
        int starting_point_inside_bin = std::lower_bound(ion_bin.begin(), ion_bin.end(), lower_rank, [&](fragment f, int rank) {
            return precursor_idx->get_rank(f.parent_id) < rank;
        }) - ion_bin.begin();

        //Update scores for all parents with fragments in the range
        for (int k = starting_point_inside_bin; k < ion_bin.size() && precursor_idx->get_rank(ion_bin[k].parent_id) <= upper_rank; ++k) {
            fragment &f = ion_bin[k]; //todo SIMD
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

void search_manager::task_search_spectrum() {


}

bool search_manager::test(unsigned int search_id) {

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

        //Update scores for all parents with fragments in the range
        for (int k = starting_point_inside_bin; k < ion_bin.size() && precursor_idx->get_rank(ion_bin[k].parent_id) <= upper_rank; ++k) {
            fragment &f = ion_bin[k]; //todo SIMD
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

