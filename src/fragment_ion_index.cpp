#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "fragment_ion_index.h"
#include "DefineConstants.h"

using namespace std;

fragment_ion_index::fragment_ion_index(precursor_index *parent_index) {

    /*
    * todo

    fragment_bins = vector<fragment_bin>(BIN_MAX_MZ + 1); //TODO remove/determine actual max #bins
    for (int i = 0; i < parent_index->get_size(); ++i) {
        spectrum *c_spectrum = parent_index->get_spectrum(i);


        //Iterate all peaks and save them as fragments in the corresponding ion mass bin
        for (int j = 0; j < c_spectrum->binned_peaks.size(); ++j) {
            int bin = c_spectrum->binned_peaks[j];
            fragment frag(c_spectrum->id, c_spectrum->binned_intensities[j]);
            fragment_bins[bin].push_back(frag);
        }
    }
    */
}


fragment_ion_index::fragment_ion_index(string path) : file_path(path) {

    load_index_from_file(file_path);

}





bool fragment_ion_index::sort_index(std::unique_ptr<precursor_index>& parent_index) {

    /*
     * Sort all bins according to parent rankings
     */

    for (fragment_bin bin : fragment_bins) {
        sort(bin.begin(), bin.end(), [&](fragment a,  fragment b){
            return parent_index->get_rank(a.parent_id) < parent_index->get_rank(b.parent_id);
        });

    }

    return false;
}

bool fragment_ion_index::load_index_from_file(const std::string& path) {

    /*
     * Read index from file
     */

    ifstream f(path, ios::in);
    string delimiter = ";";

    fragment_bins.resize(BIN_MAX_MZ + 1);

    string line;
    int c = 0;
    while (getline(f,line)) {
        ++c;
        size_t delim_pos = line.find(delimiter);
        size_t delim_right_pos = line.rfind(delimiter);

        unsigned int id = stoi(line.substr(0, delim_pos));
        int mz_bin = stoi(line.substr(delim_pos + 1, delim_right_pos - delim_pos - 1));
        float intensity = stof(line.substr(delim_right_pos + 1, string::npos));

        fragment_bins[mz_bin].emplace_back(fragment(id, intensity));

    }


    f.close();

    return false;
}

bool fragment_ion_index::save_index_to_file(const string &path) {

    /*
     * Save Index to file TODO use file_writer::
     */

    ofstream f(path, ios::out);
    string delimiter = ";";

    for (int i = 0; i < fragment_bins.size(); ++i) {

        fragment_bin  bin = fragment_bins[i];
        for (fragment frag : bin) {
            f << frag.parent_id << delimiter << i << delimiter << frag.intensity << "\n";
        }
    }

    f.close();
    return false;
}
