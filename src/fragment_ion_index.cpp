#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cassert>
#include "fragment_ion_index.h"
#include "DefineConstants.h"

using namespace std;


fragment_ion_index::fragment_ion_index() {

}


fragment_ion_index::fragment_ion_index(precursor_index *parent_index) {

    /*
    * todo

    fragment_bins = vector<fragment_bin>(BIN_MAX_MZ + 1); //TODO remove/determine actual max #bins
    for (int i = 0; i < parent_index->get_size(); ++i) {
        spectrum *c_spectrum = parent_index->get_spectrum(i);


        //Iterate all peaks and save them as fragments in the corresponding ion mz bin
        for (int j = 0; j < c_spectrum->binned_peaks.size(); ++j) {
            int bin = c_spectrum->binned_peaks[j];
            fragment frag(c_spectrum->id, c_spectrum->binned_intensities[j]);
            fragment_bins[bin].push_back(frag);
        }
    }
    */
}





fragment_ion_index::fragment_ion_index(string path) : file_path(path) {

    //load_index_from_file(file_path);
    load_index_from_binary_file(file_path);
    prepare_axv2_access();

}

bool fragment_ion_index::sort_index(std::unique_ptr<precursor_index>& parent_index) {

    /*
     * Sort all bins according to parent rankings
     */

    for (fragment_bin &bin : fragment_bins) {
        sort(bin.begin(), bin.end(), [&](fragment a,  fragment b){
            return parent_index->get_rank(a.parent_id) < parent_index->get_rank(b.parent_id);
        });
    }

    return true;
}

bool fragment_ion_index::load_index_from_file(const std::string& path) {

    /*
     * Read index from file
     */

    ifstream f(path, ios::in);
    string delimiter = ";";

    fragment_bins.clear();
    fragment_bins.resize(BIN_MAX_MZ + 1);

    string line;
    while (getline(f,line)) {
        size_t delim_pos = line.find(delimiter);
        size_t delim_right_pos = line.rfind(delimiter);

        unsigned int id = stoi(line.substr(0, delim_pos));
        int mz_bin = stoi(line.substr(delim_pos + 1, delim_right_pos - delim_pos - 1));
        float intensity = stof(line.substr(delim_right_pos + 1, string::npos));

        fragment_bins[mz_bin].emplace_back(fragment(id, intensity));

    }


    f.close();
    return true;
}

bool fragment_ion_index::save_index_to_file(const string &path) {

    /*
     * Save Index to file TODO use file_writer::
     */

    ofstream f(path, ios::out);
    f.precision(FLOAT_OUTPUT_PRECISION);

    string delimiter = ";";

    for (int i = 0; i < fragment_bins.size(); ++i) {

        fragment_bin  bin = fragment_bins[i];
        for (auto & j : bin) {
            f << j.parent_id << delimiter << i << delimiter << j.intensity << "\n";
        }
    }

    f.close();
    return true;
}

bool fragment_ion_index::load_index_from_binary_file(const string &path) {

    /*
     * Read index from binary file
     */

    ifstream f(path, ios::binary | ios::in);

    fragment_bins.clear();
    fragment_bins.resize(BIN_MAX_MZ + 1);

    string line;
    while (!f.eof()) {
        unsigned int id;
        int mz_bin;
        float intensity;

        f.read((char *) &id, sizeof(unsigned int));
        f.read((char *) &mz_bin, sizeof(int));
        f.read((char *) &intensity, sizeof(float));

        fragment_bins[mz_bin].emplace_back(fragment(id, intensity));

    }


    f.close();
    return true;
}

bool fragment_ion_index::save_index_to_binary_file(const string &path) {

    ofstream f(path, ios::binary | ios::out);


    for (int i = 0; i < fragment_bins.size(); ++i) {

        fragment_bin  bin = fragment_bins[i];
        for (auto & j : bin) {
            f.write((char *) &j.parent_id, sizeof(unsigned int)); //TODO
            f.write((char *) &i, sizeof(int));
            f.write((char *) &j.intensity, sizeof(float));
        }
    }

    f.close();
    return true;
}

bool fragment_ion_index::prepare_axv2_access() {
    frag_bins.clear();
    frag_bins.resize(fragment_bins.size());
    for (int i = 0; i < fragment_bins.size(); ++i) {
        __m256 intensity_x8;
        __m256i identity_x8;
        int ranks[8];
        for (int j = 0; j < fragment_bins[i].size(); ++j) {
            frag_bins[i].intensities.push_back(fragment_bins[i][j].intensity);
            frag_bins[i].parent_ids.push_back(fragment_bins[i][j].parent_id);
            if(j % 8 == 0 && j > 0) {
                intensity_x8 = _mm256_loadu_ps(&frag_bins[i].intensities[j-8]);
                //identity_x8 = _mm256_loadu_si256((__m256i*)&frag_bins[i].parent_ids[j-8]);
                frag_bins[i]._intensities.push_back(intensity_x8);
                //frag_bins[i]._parent_ids.push_back(identity_x8);
                //frag_bins[i]._parent_ranks.push_back(_mm256_load_si256((__m256i*)& ranks));
            }
            //ranks[j % 8] = (int) precursor_idx->get_rank(fragment_bins[i][j].parent_id);
        }
        assert(reinterpret_cast<uintptr_t>(frag_bins[i]._intensities.data()) % alignof(__m256) == 0);
        //assert(reinterpret_cast<uintptr_t>(frag_bins[i]._parent_ids.data()) % alignof(__m256i) == 0);
        //assert(reinterpret_cast<uintptr_t>(frag_bins[i]._parent_ranks.data()) % alignof(__m256i) == 0);
    }
    return true;
}
