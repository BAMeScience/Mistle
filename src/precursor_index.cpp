#include <iostream>
#include <utility>
#include <memory>
#include <fstream>
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
    return precursors[ranking.back()].mz;
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

int precursor_index::get_lower_bound(int charge, float min_mass) {
    int lb = std::lower_bound(ranking.begin(), ranking.end(), make_pair(charge, min_mass), [&](unsigned int rank, pair<int, float> charge_mass_tuple) {
        return precursors[rank] < charge_mass_tuple;
    }) - ranking.begin();

    return lb;
}


int precursor_index::get_upper_bound(int charge, float max_mass) {
    int ub = std::upper_bound(ranking.begin(), ranking.end(), make_pair(charge, max_mass), [&](pair<int, float> charge_mass_tuple, unsigned int rank) {
        return  !(precursors[rank] <= charge_mass_tuple);
    }) - ranking.begin();

    return ub - 1;
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

precursor &precursor_index::get_precursor(unsigned int id) {
    return precursors[id];
}

precursor &precursor_index::get_precursor_by_rank(unsigned int id) {
    return precursors[ranking[id]];
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
    //index_file_reader::read_file_into_precursor_index(file_path, make_shared<precursor_index>(*this));

    std::ifstream f(file_path, std::ios::in);
    std::string delimiter = ";";
    std::string line;

    if (!getline(f, line)) {
        return false;
    }

    if (line.rfind("Num: ", 0) != 0) {
        std::cerr << "Incorrect file format" << std::endl;
        return false;
    }

    //Read header
    unsigned int size = std::stoi(line.substr(5, std::string::npos)); //TODO check 4 or 5
    set_size(size);

    // Parse precursors line by line
    //precursor_idx->add_precursor_record(p);
    while (getline(f, line)) {

        size_t delim_pos = line.find(delimiter);
        unsigned int id = std::stoi(line.substr(0, delim_pos));

        size_t length = line.find(delimiter, delim_pos + 1) - delim_pos;
        unsigned int rank = std::stoi(line.substr(delim_pos + 1, length - 1));

        delim_pos = delim_pos + length;
        length = line.find(delimiter, delim_pos + 1) - delim_pos;
        float mz = std::stof(line.substr(delim_pos + 1, length - 1));

        delim_pos = delim_pos + length;
        length = line.find(delimiter, delim_pos + 1) - delim_pos;
        int charge = std::stoi(line.substr(delim_pos + 1, length - 1));
        std::string peptide = line.substr(delim_pos + length + 1, std::string::npos);

        add_precursor_record(precursor(id, rank, mz, charge, peptide));
    }

    f.close();

    if (get_size() != size) {
        std::cerr << "Wrong number of precursors" << std::endl;
    }

    //TODO delete this if not used anymore
    for (auto & precursor : precursors) {
        to_rank.push_back(precursor.rank);
    }

    return true;
}

bool precursor_index::add_precursor_record(const precursor& p) {
    precursors.emplace_back(p); //TODO test this
    ranking[p.rank] = p.id;
    return true;
}

bool precursor_index::set_size(unsigned int size) {
    precursors.reserve(size);
    ranking.resize(size);

    return true;
}

bool precursor_index::save_index_to_binary_file(const string &file_path) {

    //Saving spectrum bookmarks (precursor info)
    index_file_writer::save_precursor_index_to_binary_file(file_path, precursors);

    return true;

}

bool precursor_index::load_index_from_binary_file(const string &file_path) {
    ifstream f(file_path, ios::binary | ios::in);

    //Read header
    unsigned int size;
    f.read((char *) &size, sizeof(unsigned int));
    set_size(size);

    unsigned int id, rank;
    float mz;
    int charge;
    size_t pep_size;
    std::string peptide;
    while (f.read((char *) &id, sizeof(unsigned int))) {
        f.read((char *) &rank, sizeof(unsigned int));
        f.read((char *) &mz, sizeof(float));
        f.read((char *) &charge, sizeof(int));
        f.read((char *) &pep_size, sizeof(size_t));
        peptide.resize(pep_size);
        f.read((char *) &peptide[0], pep_size);

        add_precursor_record(precursor(id, rank, mz, charge, peptide));
    }

    f.close();

    if (get_size() != size) {
        std::cerr << "Wrong number of precursors" << std::endl;
        std::cout << get_size() << " " << size << endl;
    }

    //TODO delete this if not used anymore
    for (auto & precursor : precursors) {
        to_rank.push_back(precursor.rank);
    }

    return true;

}


