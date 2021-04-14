#ifndef SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#define SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#include <vector>
#include <memory>
#include "spectrum.h"



struct precursor {
    /*
     * Key values
     */
    unsigned int id;
    unsigned int rank; //unsure if needed here
    float mz;
    int charge;
    std::string peptide;

    /*
     * special cases
     */

    unsigned long offset_begin;
    unsigned long offset_end;
    std::string name;

    precursor() {};
    precursor(unsigned int id, float mass, int charge, std::string peptide="") : id(id), mz(mass), charge(charge), peptide(peptide) {};
    precursor(unsigned int id, unsigned int rank, float mass, int charge, std::string peptide="") : id(id), rank(rank), mz(mass), charge(charge), peptide(peptide) {};

    bool operator<(const precursor &other) const {
        return charge < other.charge || (charge == other.charge && mz < other.mz);
    };

    bool operator<(std::pair<int, float> charge_mass_tuple) const {
        return charge < charge_mass_tuple.first || (charge == charge_mass_tuple.first && mz < charge_mass_tuple.second);
    }

    bool operator<=(std::pair<int, float> charge_mass_tuple) const {
        return charge < charge_mass_tuple.first || (charge == charge_mass_tuple.first && mz <= charge_mass_tuple.second);
    }


};


class precursor_index {

    // Contains all spectrum bookmarks (precursors), sorted first by charge, then by precursor mz
    std::vector<precursor> precursors;
    std::vector<unsigned int> ranking;
    unsigned int id_counter = 0;


public:
    std::vector<int> to_rank;
    precursor_index(); //Init empty index
    precursor& record_new_precursor(const std::shared_ptr<spectrum>& spec);
    precursor& record_new_precursor(float mz, int charge, std::string peptide);
    bool add_precursor_record(const precursor& p);

    bool sort_index();
    bool save_index_to_file(const std::string &file_path);
    bool save_index_to_binary_file(const std::string &file_path);
    bool load_index_from_file(const std::string &file_path);
    bool load_index_from_binary_file(const std::string &file_path);


    int get_size();
    bool set_size(unsigned int size);
    int get_lower_bound(int charge, float min_mass);
    int get_upper_bound(int charge, float max_mass);
    float get_max_precursor_mass();

    precursor& get_precursor(unsigned int id);
    precursor& get_precursor_by_rank(unsigned int id);
    unsigned int get_rank(unsigned int id);

};


#endif //SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
