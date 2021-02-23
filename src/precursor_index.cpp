#include <iostream>
#include <utility>
#include <memory>
#include "precursor_index.h"


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

    id_to_rank = vector<unsigned int>(ranking.size(), 0);
    for (int i = 0; i < ranking.size(); ++i) {
        unsigned int id = ranking[i];
        id_to_rank[id] = i;
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
    return id_to_rank[id];
}
