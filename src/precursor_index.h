#ifndef SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#define SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#include <vector>
#include "spectrum.h"



struct precursor {
    /*
     * Key values
     */
    unsigned int id;
    int rank; //unsure if needed here
    float mass;
    int charge;
    std::string peptide;

    /*
     * special cases
     */

    unsigned long offset_begin;
    unsigned long offset_end;
    std::string name;

    precursor() {};
    precursor(unsigned int id, float mass, int charge, std::string peptide="") : id(id), mass(mass), charge(charge), peptide(peptide) {};

    bool operator<(const precursor &other) const {
        return charge < other.charge || (charge == other.charge && mass < other.mass);
    };

};


class precursor_index {

    // Contains all spectrum references, sorted first by charge, then by precursor mass
    std::vector<precursor> precursors;
    std::vector<unsigned int> ranking;
    std::vector<unsigned int> id_to_rank; //inverse of ranking
    unsigned int id_counter = 0;

public:

    precursor_index(); //Init empty index
    precursor& record_new_precursor(spectrum *spectrum);
    precursor& record_new_precursor(float mz, int charge, std::string peptide);
    bool sort_index();

    int get_size();
    int get_lower_bound(int charge, float mass);
    int get_upper_bound(int charge, float max_mass);
    float get_max_precursor_mass();

    precursor& get_precursor(int i);

};


#endif //SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
