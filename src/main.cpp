#include <iostream>
#include <numeric>
#include "spectrum.h"
#include "msp_reader.h"
#include "scores.h"
#include "spectral_search.h"
#include "library.h"


int main() {
    cout << "Welcome, welcome" << endl;
    string s = "string \2s dssd\n";
    cout << "2 test: " << "\2" << s << endl;
    string msp_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_Complete_simulatedSpectra\pyrofur_2019.msp)";
    string mgf_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_SearchFile\pfu_velos.mgf)";
    //vector<spectrum*> library = msp_reader::read_file(R"(C:\Users\ynowatzk\Desktop\data\9MM\simulated_spectra\Brevibacillus+laterosporus.msp)");


    library *search_lib = new library(mgf_file);

    spectrum *one,*two;
    for (spectrum *s : search_lib->spectrum_list) {
        if (s->name == "824.836730957031_212.9232_2") {
            one = s;
        }
        if (s->name == "508.785949707031_252.8919_2") {
            two = s;
        }
    }


    spectral_search search(search_lib);

    library *lib = new library(msp_file);
    //library *lib = new library();


    spectrum *one_ST, *two_ST;
    for (spectrum *s : lib->spectrum_list) {
        if (s->name == "INKAIEFPIDDLKK/2") {
            one_ST = s;
            cout << "one: st: 0.223 my: " << scores::dot_product(one->bins, one_ST->bins) << endl;
        }
        if (s->name == "KHLEQHPK/2") {
            two_ST = s;
            cout << "two: st: 0.774 my: " << scores::dot_product(two->bins, two_ST->bins) << endl;
        }
    }

    exit(12);

    search.search_target_library(lib);

    vector<match> matches = search.get_results();

    for (int i = 0; i < 10; ++i) {
        cout << matches[i].query_spectrum->name << " " << matches[i].matched_spectrum->peptide << " " << matches[i].dot_product << endl;
    }

    search.save_results_to_file("./out.csv");


    return 0;
}
