#include <iostream>
#include <numeric>
#include "spectrum.h"
#include "msp_reader.h"
#include "scores.h"
#include "library.h"

using namespace std;

int main() {
    cout << "Welcome, welcome" << endl;
    string s = "string \2s dssd\n";
    cout << "2 test: " << "\2" << s << endl;
    string msp_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_Complete_simulatedSpectra\pyrofur_2019.msp)";
    string mgf_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_SearchFile\pfu_velos.mgf)";
    //vector<spectrum*> library = msp_reader::read_file(R"(C:\Users\ynowatzk\Desktop\data\9MM\simulated_spectra\Brevibacillus+laterosporus.msp)");


    library *search_lib = new library(mgf_file);
    library *lib = new library(msp_file);
    //library *lib = new library();


    spectrum *first = search_lib->spectrum_list[0];
    cout << first->peptide << endl;
    cout << "peak_positions: ";
    for (float p : first->peak_positions) {
        cout << p << " ";
    }
    cout << endl << "inten: ";
    for (float i : first->intensities) {
        cout << i << " ";
    }
    cout << endl << "bins: ";
    for (float i : first->bins) {
        cout << i << " ";
    }
    cout << endl;

    cout << "accumulated " << accumulate(first->bins.begin(), first->bins.end(), 0.0) << endl;

    cout << "# spec: " << lib->spectrum_list.size() << endl;
    cout << "# search specs: " << search_lib->spectrum_list.size() << endl;

    cout << "Target spectrum: " << endl;
    spectrum *target = search_lib->spectrum_list[0];
    cout << "target: " << target->name << " " << target->precursor_mass << " " << target->charge << endl;
    float max_dot = 0;

    spectrum *best_candidate = search_lib->spectrum_list[0];
    for (int i = 0; i < lib->spectrum_list.size(); ++i) {
        spectrum *candidate = lib->spectrum_list[i];
        if (abs(target->precursor_mass - candidate->precursor_mass) > 3 || target->charge!=candidate->charge) {
            continue;
        }
        if (scores::dot_product(candidate->bins, target->bins) > max_dot) {
            max_dot = scores::dot_product(candidate->bins, target->bins);
            best_candidate = candidate;
        }
        if (candidate->peptide == "INKAIEFPIDDLKK") {
            cout << "spectrast candidate:" << candidate->peptide << " " << candidate->precursor_mass << " " << candidate->charge << " dot: " << scores::dot_product(search_lib->spectrum_list[0]->bins, candidate->bins) << endl;
        }
    }
    cout << "my candidate: " << best_candidate->peptide << " " << best_candidate->precursor_mass << " " << best_candidate->charge << " dot: " << max_dot << endl;

    delete lib;
    delete search_lib;

    return 0;
}
