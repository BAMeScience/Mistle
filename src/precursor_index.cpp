#include <iostream>
#include <utility>
#include <memory>
#include "precursor_index.h"
#include "index_file_writer.h"
#include "index_file_reader.h"


using namespace std;

int precursor_index::get_size() {
    return precursors.size();
}

/*spectrum *precursor_index::get_spectrum(int i) {
    return spectra[i];
}*/

float precursor_index::get_max_precursor_mass() {
    return precursors[ranking.back()].mass;
}
/*
//TODO add self chosen bounds
int precursor_index::get_lower_bound(int charge, float min_mass) {

    int lb = std::lower_bound(spectra.begin(), spectra.end(), make_pair(charge, min_mass), [](spectrum *s, pair<int, float> charge_mass_tuple) {
        return *s < charge_mass_tuple;
    }) - spectra.begin();

    return lb;
}


int precursor_index::get_upper_bound(int charge, float max_mass) {

    int ub = std::upper_bound(spectra.begin(), spectra.end(), make_pair(charge, max_mass), [](pair<int, float> charge_mass_tuple, spectrum *s) {
        return  !(*s <= charge_mass_tuple);
    }) - spectra.begin();

    return ub - 1;
}
*/

precursor_index::precursor_index() {

}

bool precursor_index::sort_index() {

    if (!(id_counter == precursors.size() && id_counter == ranking.size())) {
        cerr << "Number of recorded precursors does not match precursor ids :: Required to warrant correct mapping of id to precursors" << endl;
        return false;
    }


    sort(ranking.begin(), ranking.end(), [&](unsigned int a,  unsigned int b) {
        return precursors[a] < precursors[b];
    });

    for (int i = 0; i < ranking.size(); ++i) {
        unsigned int id = ranking[i];
        precursors[id].rank = i;
    }

    return true;
}


precursor &precursor_index::get_precursor(int i) {
    return precursors[i];
}

precursor &precursor_index::record_new_precursor(const shared_ptr<spectrum>& spec) {
    precursors.emplace_back(precursor(id_counter, spec->precursor_mass, spec->charge, spec->peptide));
    ranking.push_back(id_counter);
    ++id_counter;
    return precursors.back();
}

precursor &precursor_index::record_new_precursor(float mz, int charge, std::string peptide) {
    precursors.emplace_back(precursor(id_counter, mz, charge, peptide));
    ranking.push_back(id_counter);
    ++id_counter;
    return precursors.back();
}

unsigned int precursor_index::get_rank(unsigned int id) {
    return precursors[id].rank;
}

bool precursor_index::save_index_to_file(const string &file_path) {

    //Saving spectrum bookmarks (precursor info)
    index_file_writer::save_precursor_index(file_path, precursors);


    return true;
}

bool precursor_index::load_index_from_file(const string &file_path) {
    index_file_reader::read_file_into_precursor_index(file_path, make_shared<precursor_index>(*this));

    return true;
}

bool precursor_index::add_precursor_record(precursor &p) {
    precursors.emplace_back(p); //TODO test this
    ranking[p.rank] = p.id;
    return true;
}

bool precursor_index::set_size(unsigned int size) {
    precursors.reserve(size);
    ranking.reserve(size);

    return true;
}


